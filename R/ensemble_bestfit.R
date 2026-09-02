#' Create prediction and lag-contribution ensembles from best-fit models
#'
#' Builds one or more model-level ensembles from out-of-fold predictions
#' produced by `find_bestfit()` and, optionally, lag-level ensembles from
#' lag-specific contributions. Automatic lag decomposition uses
#' `compute_ecilag()`.
#'
#' Exposure histories supplied through `data` must be chronological, from the
#' earliest to the most recent observation within each group. Their conversion
#' to the retrospective lag representation is handled by `compute_ecilag()`.
#'
#' @param bestfit Data frame returned by `find_bestfit()`. It must contain the
#'   model identifier and the metric selected by `weight_metric`. The current
#'   `find_bestfit()` also stores `family`, `outcome_type`, `rank_metric`, and
#'   `threshold` as attributes, which are inherited here when appropriate.
#' @param predictions Optional data frame of out-of-fold predictions. If `NULL`,
#'   `attr(bestfit, "predictions")` is used. Required columns are `model_col`,
#'   `id_cols`, `observed`, and `predicted`. For binomial outcomes, `predicted`
#'   must contain predicted probabilities in `[0, 1]`, not already-thresholded
#'   class labels.
#' @param lag_data Optional data frame of lag-specific contributions. Required
#'   columns are `model_col`, `lag_group_cols`, and `contribution`.
#' @param fit_list Optional named list of fitted models. If `NULL`,
#'   `attr(bestfit, "fits")` is used for automatic lag decomposition.
#' @param data Optional long-format exposure data used for automatic lag
#'   decomposition. Must contain `group`, `time`, and the requested variables.
#' @param group Character scalar naming the grouping column in `data`.
#' @param var Optional character vector of exposure variables used for automatic
#'   lag decomposition.
#' @param compute_ecilag_args Named list of additional arguments passed to
#'   `compute_ecilag()`. It cannot override `data`, `group`, `fit`, `var`, or
#'   `profile`.
#' @param ensemble_scope Character scalar. One of `"model"`, `"lag"`, or
#'   `"both"`.
#' @param method One or more ensemble methods selected from `"best"`,
#'   `"hard_voting"`, `"unweighted"`, `"weighted"`, and `"stacked"`.
#'   Multiple methods can be evaluated in a single call, for example
#'   `method = c("unweighted", "weighted", "stacked")`. `"stacked"` and
#'   `"hard_voting"` are available only for `ensemble_scope = "model"`.
#'   `"hard_voting"` additionally requires a binomial outcome.
#' @param top_n Positive integer number of top-ranked models included in each
#'   requested ensemble.
#' @param model_ids Optional vector of model identifiers. If supplied, `top_n`
#'   is ignored, but selected models are still ordered by `weight_metric`.
#' @param weight_metric Optional character scalar naming a performance metric in
#'   `bestfit` used to rank/select models and, for `method = "weighted"`, derive
#'   model weights. If `NULL`, the function first inherits
#'   `attr(bestfit, "rank_metric")`; if unavailable, the default is `"CCC"` for
#'   non-binary outcomes and `"ROC_AUC"` for binomial outcomes. Metric direction
#'   is determined automatically; higher-is-better metrics are maximized and
#'   error/loss metrics are minimized.
#' @param threshold Optional numeric probability threshold strictly between 0
#'   and 1. For a binomial outcome, it converts predicted probabilities into
#'   classes 0/1 for classification metrics and converts individual-model
#'   probabilities into votes for `method = "hard_voting"`. It is valid for
#'   probabilities, not predictions already reduced to class labels 0/1. It
#'   does not affect ROC AUC, Brier score, or Log Loss. If `NULL`, the function
#'   inherits `attr(bestfit, "threshold")`; if unavailable, `0.5` is used.
#' @param weight_transform Character scalar. One of `"softmax"`, `"positive"`,
#'   `"rank_inverse"`, or `"uniform"`. Used only by `method = "weighted"`.
#' @param stacking_model Character scalar. `"ridge"` or `"lm"`.
#' @param stack_objective Character scalar. `"regularized"`, `"min_error"`, or
#'   `"constrained"`.
#' @param lambda Non-negative ridge penalty.
#' @param stack_intercept Logical. Include an intercept in stacking.
#' @param model_col Character scalar naming the model identifier column.
#' @param id_cols Character vector identifying prediction rows.
#' @param lag_group_cols Character vector defining each lag-ensemble unit.
#'   Include an epidemic/group column, for example
#'   `c("epi_id", "var", "lag")`, when lag contributions are group-specific.
#' @param lg_strategy Character scalar. One of `"requested_available"`,
#'   `"model"`, or `"strict"`.
#' @param rl_weights Logical. Renormalize model weights within each lag unit.
#' @param verbose Logical. Print informative messages.
#'
#' @return A list containing `ensemble_summary`, `ensemble_predictions`,
#'   `ensemble_by_lag`, `lag_data`, `model_weights`, `selected_models`,
#'   `metric_warnings`, `method`, `ensemble_scope`, `lg_strategy`, `rl_weights`,
#'   `family`, `outcome_type`, `weight_metric`, `metric_direction`, and
#'   `threshold`. Outputs that vary by ensemble method contain an explicit
#'   `method` column, allowing several methods to be compared from one call.
#'
#' @details
#' Model ranking and weighting use the common EpiExposure performance-metric
#' infrastructure also used by `find_bestfit()`. Supported non-binary metrics
#' are CCC, Cb, rho, RMSE, and MAE. Supported binomial metrics are ROC AUC,
#' Brier score, Log Loss, Accuracy, Balanced Accuracy, Sensitivity, Specificity,
#' F1, MCC, and Precision. The direction of each metric is resolved internally,
#' so no higher-is-better flag is required.
#'
#' For binomial outcomes, `best`, `unweighted`, `weighted`, and `stacked`
#' combine or generate predicted probabilities. Classification metrics are then
#' computed after applying `threshold`. `hard_voting` is different: each model
#' probability is thresholded first and the final class is obtained by majority
#' vote. With an even number of models, an exact vote tie is broken using the
#' class predicted by the highest-ranked selected model. The returned
#' `vote_fraction` is the fraction of positive votes and is not treated as a
#' calibrated predicted probability. Consequently, probability-based metrics
#' (ROC AUC, Brier, and Log Loss) are returned as `NA` for hard voting.
#'
#' Lag ensembles combine lag-specific contributions, not raw DLNM coefficients:
#' \deqn{C_{ens,l}=\sum_m \alpha_m C_{m,l}}
#'
#' `lag_data` must contain exactly one row per combination of `model_col` and
#' `lag_group_cols`. If automatic decomposition returns several epidemics,
#' include the grouping column in `lag_group_cols`.
#'
#' When `compute_ecilag_args` requests uncertainty, this function combines the
#' summarized `by_lag` estimates returned by `compute_ecilag()`. It does not
#' combine `by_lag_samples` draw by draw and therefore does not return ensemble
#' uncertainty intervals.
#'
#' Model-level ensembles combine out-of-fold predictions that have already been
#' generated by the underlying models. Therefore, deterministic and
#' uncertainty-related behaviour is inherited from the prediction workflows
#' used to create the input `predictions` object.
#'
#' @export
ensemble_bestfit <- function(
    bestfit,
    predictions = NULL,
    lag_data = NULL,
    fit_list = NULL,
    data = NULL,
    group = "epi_id",
    var = NULL,
    compute_ecilag_args = list(),
    ensemble_scope = c("model", "lag", "both"),
    method = "unweighted",
    top_n = 3,
    model_ids = NULL,
    weight_metric = NULL,
    threshold = NULL,
    weight_transform = c("softmax", "positive", "rank_inverse", "uniform"),
    stacking_model = c("ridge", "lm"),
    stack_objective = c("regularized", "min_error", "constrained"),
    lambda = 1e-6,
    stack_intercept = TRUE,
    model_col = "model_id",
    id_cols = c("group", "fold"),
    lag_group_cols = c("var", "lag"),
    lg_strategy = c("requested_available", "model", "strict"),
    rl_weights = TRUE,
    verbose = TRUE
) {
  ensemble_scope <- match.arg(ensemble_scope)
  weight_transform <- match.arg(weight_transform)
  stacking_model <- match.arg(stacking_model)
  stack_objective <- match.arg(stack_objective)
  lg_strategy <- match.arg(lg_strategy)

  allowed_methods <- c(
    "unweighted",
    "weighted",
    "hard_voting",
    "best",
    "stacked"
  )

  if (!is.character(method) || !length(method) || anyNA(method) ||
      any(!nzchar(method))) {
    stop(
      "`method` must contain one or more non-empty ensemble method names.",
      call. = FALSE
    )
  }

  invalid_methods <- setdiff(method, allowed_methods)
  if (length(invalid_methods)) {
    stop(
      "Unsupported `method`: ", paste(invalid_methods, collapse = ", "),
      ". Available methods are: ", paste(allowed_methods, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # Repeated methods have no additional meaning. Preserve the user's requested
  # order while evaluating each distinct method only once.
  methods <- unique(method)

  response_col <- "observed"
  prediction_col <- "predicted"
  lag_col <- "lag"
  lag_contribution_col <- "contribution"
  lag_weight_col <- "weight"

  valid_name <- function(x) {
    is.character(x) && length(x) >= 1L && !anyNA(x) &&
      all(nzchar(x)) && !anyDuplicated(x)
  }

  valid_flag <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  if (!is.data.frame(bestfit)) stop("`bestfit` must be a data.frame.", call. = FALSE)

  if (!valid_name(model_col) || length(model_col) != 1L) {
    stop("`model_col` must be one non-empty column name.", call. = FALSE)
  }

  if (!valid_name(id_cols)) {
    stop("`id_cols` must contain unique non-empty names.", call. = FALSE)
  }

  if (!valid_name(lag_group_cols)) {
    stop("`lag_group_cols` must contain unique non-empty names.", call. = FALSE)
  }

  if (!is.character(group) || length(group) != 1L ||
      is.na(group) || !nzchar(group)) {
    stop("`group` must be one non-empty column name.", call. = FALSE)
  }

  if (!is.null(var) &&
      (!is.character(var) || !length(var) || anyNA(var) ||
       any(!nzchar(var)) || anyDuplicated(var))) {
    stop(
      "`var` must be NULL or a unique non-empty character vector.",
      call. = FALSE
    )
  }

  if (!model_col %in% names(bestfit)) {
    stop(
      "`model_col` ('", model_col, "') was not found in `bestfit`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(bestfit[[model_col]])) {
    stop("`model_col` must uniquely identify rows in `bestfit`.", call. = FALSE)
  }

  if (!is.numeric(top_n) || length(top_n) != 1L ||
      !is.finite(top_n) || top_n < 1) {
    stop("`top_n` must be a positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0) {
    stop("`lambda` must be a non-negative numeric scalar.", call. = FALSE)
  }

  if (!valid_flag(stack_intercept)) {
    stop("`stack_intercept` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!valid_flag(rl_weights)) {
    stop("`rl_weights` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!valid_flag(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.list(compute_ecilag_args)) {
    stop("`compute_ecilag_args` must be a named list.", call. = FALSE)
  }

  if (length(compute_ecilag_args) &&
      (is.null(names(compute_ecilag_args)) || anyNA(names(compute_ecilag_args)) ||
       any(!nzchar(names(compute_ecilag_args))) ||
       anyDuplicated(names(compute_ecilag_args)))) {
    stop(
      "`compute_ecilag_args` must have unique non-empty names.",
      call. = FALSE
    )
  }

  reserved_args <- c("data", "group", "fit", "var", "profile")
  forbidden_args <- intersect(names(compute_ecilag_args), reserved_args)
  if (length(forbidden_args)) {
    stop(
      "`compute_ecilag_args` cannot override: ",
      paste(forbidden_args, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.null(model_ids) && anyDuplicated(model_ids)) {
    stop("`model_ids` must contain unique identifiers.", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Resolve family/outcome metadata produced by find_bestfit(). If an older
  # bestfit object is supplied, attempt a conservative recovery from stored fits
  # before failing explicitly. Never infer family from the numerical range of y.
  # -------------------------------------------------------------------------
  family_input <- attr(bestfit, "family")

  recover_family_from_fits <- function(fits) {
    if (!is.list(fits) || !length(fits)) return(NULL)

    recovered <- vapply(
      fits,
      function(fit) {
        value <- attr(fit, "epiexposure_family_name")
        if (is.null(value) || !length(value)) NA_character_ else as.character(value[[1]])
      },
      character(1)
    )

    recovered <- unique(recovered[!is.na(recovered) & nzchar(recovered)])
    if (!length(recovered)) return(NULL)
    if (length(recovered) > 1L) {
      stop(
        "Stored fitted models contain inconsistent family metadata: ",
        paste(recovered, collapse = ", "), ".",
        call. = FALSE
      )
    }
    recovered[[1]]
  }

  if (is.null(family_input)) {
    candidate_fit_list <- if (!is.null(fit_list)) fit_list else attr(bestfit, "fits")
    family_input <- recover_family_from_fits(candidate_fit_list)
  }

  if (is.null(family_input)) {
    stop(
      "`bestfit` does not contain family metadata. Use the current ",
      "`find_bestfit()` (which stores attr(bestfit, 'family')) or provide a ",
      "`fit_list` whose fitted models contain `epiexposure_family_name` metadata.",
      call. = FALSE
    )
  }

  family_name <- .resolve_family_name(family_input)
  outcome_type <- .resolve_outcome_type(family_name)

  stored_outcome_type <- attr(bestfit, "outcome_type")
  if (!is.null(stored_outcome_type) &&
      length(stored_outcome_type) == 1L &&
      !is.na(stored_outcome_type) &&
      !identical(as.character(stored_outcome_type), outcome_type)) {
    stop(
      "`bestfit` contains inconsistent family/outcome metadata: family = '",
      family_name, "' implies outcome_type = '", outcome_type,
      "', but attr(bestfit, 'outcome_type') is '", stored_outcome_type, "'.",
      call. = FALSE
    )
  }

  if (identical(outcome_type, "ordinal")) {
    stop(
      "`ensemble_bestfit()` does not yet define performance metrics for ",
      "`family = 'ordinal'`. Ordinal-specific ensemble evaluation must be ",
      "defined before this family can be used here.",
      call. = FALSE
    )
  }

  available_metrics <- .available_metrics(family_name)

  if (is.null(weight_metric)) {
    inherited_metric <- attr(bestfit, "rank_metric")
    if (is.character(inherited_metric) && length(inherited_metric) == 1L &&
        !is.na(inherited_metric) && nzchar(inherited_metric) &&
        inherited_metric %in% available_metrics &&
        inherited_metric %in% names(bestfit)) {
      weight_metric <- inherited_metric
    } else {
      weight_metric <- if (identical(outcome_type, "binary")) "ROC_AUC" else "CCC"
    }
  }

  if (!is.character(weight_metric) || length(weight_metric) != 1L ||
      is.na(weight_metric) || !nzchar(weight_metric)) {
    stop(
      "`weight_metric` must be NULL or one non-empty metric name.",
      call. = FALSE
    )
  }

  if (!weight_metric %in% available_metrics) {
    stop(
      "`weight_metric = ", sQuote(weight_metric), "` is not available for family '",
      family_name, "'. Available metrics are: ",
      paste(available_metrics, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!weight_metric %in% names(bestfit)) {
    stop(
      "`weight_metric` ('", weight_metric, "') was not found in `bestfit`.",
      call. = FALSE
    )
  }

  if (!is.numeric(bestfit[[weight_metric]])) {
    stop("`weight_metric` must identify a numeric column in `bestfit`.", call. = FALSE)
  }

  metric_direction <- .metric_direction(weight_metric)

  if (is.null(threshold)) {
    inherited_threshold <- attr(bestfit, "threshold")
    if (is.numeric(inherited_threshold) && length(inherited_threshold) == 1L &&
        is.finite(inherited_threshold)) {
      threshold <- inherited_threshold
    } else {
      threshold <- 0.5
    }
  }

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0 || threshold >= 1) {
    stop(
      "`threshold` must be one finite numeric probability strictly between 0 and 1.",
      call. = FALSE
    )
  }

  if (ensemble_scope %in% c("lag", "both") && "stacked" %in% methods) {
    stop(
      "`method = 'stacked'` is supported only for `ensemble_scope = 'model'`.",
      call. = FALSE
    )
  }

  if (ensemble_scope %in% c("lag", "both") && "hard_voting" %in% methods) {
    stop(
      "`method = 'hard_voting'` is supported only for `ensemble_scope = 'model'`.",
      call. = FALSE
    )
  }

  if ("hard_voting" %in% methods && !identical(outcome_type, "binary")) {
    stop(
      "`method = 'hard_voting'` is available only for binomial outcomes. ",
      "For non-binary outcomes, use 'best', 'unweighted', 'weighted', or 'stacked'.",
      call. = FALSE
    )
  }

  if (verbose) {
    message(
      "Ensemble performance: family = ", family_name,
      ", outcome_type = ", outcome_type,
      ", weight_metric = ", weight_metric,
      " (", metric_direction, ")."
    )

    if (identical(outcome_type, "binary") &&
        ensemble_scope %in% c("model", "both")) {
      message(
        "Classification threshold = ", format(threshold),
        ". It is applied to predicted probabilities to generate 0/1 classes",
        if ("hard_voting" %in% methods) " and model votes" else "",
        "; it is not intended for predictions already reduced to class labels. ",
        "ROC_AUC, Brier, and LogLoss are unaffected by the threshold."
      )
    }
  }

  project_simplex <- function(v) {
    v <- as.numeric(v)
    if (!length(v)) return(v)
    if (all(!is.finite(v))) return(rep(1 / length(v), length(v)))
    v[!is.finite(v)] <- 0
    u <- sort(v, decreasing = TRUE)
    cssv <- cumsum(u)
    candidates <- which(u + (1 - cssv) / seq_along(u) > 0)
    if (!length(candidates)) return(rep(1 / length(v), length(v)))
    rho <- max(candidates)
    theta <- (cssv[rho] - 1) / rho
    weights <- pmax(v - theta, 0)
    if (sum(weights) <= 0) rep(1 / length(v), length(v)) else weights / sum(weights)
  }

  rank_models <- function(model_table) {
    metric <- model_table[[weight_metric]]

    primary <- if (identical(metric_direction, "maximize")) {
      -metric
    } else {
      metric
    }

    tie_rank <- if ("rank" %in% names(model_table) && is.numeric(model_table$rank)) {
      model_table$rank
    } else {
      seq_len(nrow(model_table))
    }

    ord <- order(
      primary,
      tie_rank,
      as.character(model_table[[model_col]]),
      na.last = NA
    )

    model_table[ord, , drop = FALSE]
  }

  if (!is.null(model_ids)) {
    missing_models <- setdiff(model_ids, bestfit[[model_col]])
    if (length(missing_models)) {
      stop(
        "The following `model_ids` were not found: ",
        paste(missing_models, collapse = ", "), ".",
        call. = FALSE
      )
    }

    selected_models <- bestfit[
      bestfit[[model_col]] %in% model_ids,
      ,
      drop = FALSE
    ]
    selected_models <- rank_models(selected_models)
  } else {
    selected_models <- utils::head(rank_models(bestfit), top_n)
  }

  if (!nrow(selected_models)) {
    stop("No models with finite `weight_metric` values were selected for the ensemble.", call. = FALSE)
  }

  if (all(!is.finite(selected_models[[weight_metric]]))) {
    stop(
      "No finite `weight_metric` values are available for selected models.",
      call. = FALSE
    )
  }

  selected_ids <- selected_models[[model_col]]
  selected_ids_chr <- as.character(selected_ids)

  if (verbose) {
    message(
      "Selected ", length(selected_ids), " model(s): ",
      paste(selected_ids, collapse = ", "),
      ". Methods: ", paste(methods, collapse = ", "), "."
    )
  }

  compute_model_weights <- function(model_table, ensemble_method) {
    n_models <- nrow(model_table)
    ids <- as.character(model_table[[model_col]])

    if (ensemble_method == "unweighted") {
      weights <- rep(1 / n_models, n_models)
    } else if (ensemble_method == "best") {
      weights <- c(1, rep(0, n_models - 1L))
    } else if (ensemble_method == "weighted") {
      metric <- model_table[[weight_metric]]
      finite <- is.finite(metric)

      if (!any(finite)) {
        stop("No finite values are available for weighting.", call. = FALSE)
      }

      if (weight_transform == "uniform") {
        weights <- rep(1 / n_models, n_models)
      } else if (weight_transform == "rank_inverse") {
        weights <- 1 / seq_len(n_models)
        weights <- weights / sum(weights)
      } else if (weight_transform == "positive") {
        if (identical(metric_direction, "maximize")) {
          utility <- pmax(metric, 0)
        } else {
          # All currently supported minimize metrics (RMSE, MAE, Brier,
          # LogLoss) are non-negative losses. Their reciprocal is therefore a
          # natural positive utility in which smaller loss receives more weight.
          utility <- 1 / pmax(metric, .Machine$double.eps)
        }
        utility[!finite | !is.finite(utility)] <- 0
        weights <- if (sum(utility) > 0) {
          utility / sum(utility)
        } else {
          rep(1 / n_models, n_models)
        }
      } else {
        score <- if (identical(metric_direction, "maximize")) metric else -metric
        score[!finite] <- -Inf
        shifted <- score - max(score[finite])
        exp_score <- exp(shifted)
        exp_score[!is.finite(exp_score)] <- 0
        weights <- exp_score / sum(exp_score)
      }
    } else {
      stop(
        "Internal weighting is not defined for method '", ensemble_method, "'.",
        call. = FALSE
      )
    }

    names(weights) <- ids
    weights
  }

  build_prediction_matrix <- function(prediction_data) {
    required <- c(model_col, id_cols, response_col, prediction_col)
    missing <- setdiff(required, names(prediction_data))
    if (length(missing)) {
      stop(
        "Missing prediction columns: ", paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }

    prediction_data <- prediction_data[
      as.character(prediction_data[[model_col]]) %in% selected_ids_chr,
      ,
      drop = FALSE
    ]

    if (!nrow(prediction_data)) {
      stop("No predictions were found for selected models.", call. = FALSE)
    }

    if (!is.numeric(prediction_data[[response_col]]) ||
        !is.numeric(prediction_data[[prediction_col]])) {
      stop("`observed` and `predicted` must be numeric.", call. = FALSE)
    }

    key_cols <- c(model_col, id_cols)
    if (anyDuplicated(prediction_data[key_cols])) {
      stop(
        "`predictions` contains duplicate rows for a model and prediction identifier.",
        call. = FALSE
      )
    }

    make_key <- function(d, columns) {
      if (length(columns) == 1L) return(as.character(d[[columns]]))
      apply(d[, columns, drop = FALSE], 1L, paste, collapse = "||")
    }

    prediction_data$.obs_key <- make_key(prediction_data, id_cols)
    observed_sets <- split(prediction_data[[response_col]], prediction_data$.obs_key)
    inconsistent <- vapply(
      observed_sets,
      function(z) length(unique(z[is.finite(z)])) > 1L,
      logical(1)
    )

    if (any(inconsistent)) {
      stop(
        "Different observed values were found for the same prediction identifier.",
        call. = FALSE
      )
    }

    obs_df <- prediction_data[, c(".obs_key", id_cols, response_col), drop = FALSE]
    obs_df <- obs_df[!duplicated(obs_df$.obs_key), , drop = FALSE]
    obs_keys <- obs_df$.obs_key

    pred_mat <- matrix(
      NA_real_,
      nrow = length(obs_keys),
      ncol = length(selected_ids),
      dimnames = list(obs_keys, paste0("model_", selected_ids_chr))
    )

    for (j in seq_along(selected_ids_chr)) {
      temporary <- prediction_data[
        as.character(prediction_data[[model_col]]) == selected_ids_chr[j],
        ,
        drop = FALSE
      ]
      row_index <- match(temporary$.obs_key, obs_keys)
      pred_mat[row_index, j] <- temporary[[prediction_col]]
    }

    complete <- stats::complete.cases(pred_mat) & is.finite(obs_df[[response_col]])
    if (sum(complete) < 2L) {
      stop(
        "Fewer than two complete OOF observations are shared by selected models.",
        call. = FALSE
      )
    }

    pred_mat <- pred_mat[complete, , drop = FALSE]
    obs_df <- obs_df[complete, , drop = FALSE]
    observed <- obs_df[[response_col]]

    if (identical(outcome_type, "binary")) {
      if (!all(unique(observed) %in% c(0, 1))) {
        stop(
          "Binomial ensemble performance requires observed outcomes coded as 0/1.",
          call. = FALSE
        )
      }

      if (any(pred_mat < 0 | pred_mat > 1)) {
        stop(
          "Binomial ensemble predictions must be probabilities between 0 and 1. ",
          "The supplied OOF predictions contain values outside this range.",
          call. = FALSE
        )
      }

      if (all(pred_mat %in% c(0, 1))) {
        warning(
          "Selected binomial-model predictions contain only 0/1 values. ",
          "`threshold` is intended for predicted probabilities; probability-based ",
          "metrics will not retain probability-resolution information.",
          call. = FALSE
        )
      }
    }

    list(
      pred_mat = pred_mat,
      obs_df = obs_df,
      obs = observed
    )
  }

  compute_stacking <- function(pred_mat, observed) {
    fit_stack <- function(X, y) {
      design <- if (stack_intercept) cbind("(Intercept)" = 1, X) else X

      if (stacking_model == "lm" || stack_objective == "min_error") {
        coefficients <- stats::lm.fit(design, y)$coefficients
        coefficients[!is.finite(coefficients)] <- 0
        return(coefficients)
      }

      penalty <- diag(lambda, ncol(design))
      if (stack_intercept) penalty[1, 1] <- 0

      tryCatch(
        as.numeric(
          solve(
            crossprod(design) + penalty,
            crossprod(design, y)
          )
        ),
        error = function(e) rep(0, ncol(design))
      )
    }

    leave_one_out <- numeric(nrow(pred_mat))

    for (i in seq_len(nrow(pred_mat))) {
      train <- setdiff(seq_len(nrow(pred_mat)), i)
      coefficients <- fit_stack(
        pred_mat[train, , drop = FALSE],
        observed[train]
      )
      row_design <- if (stack_intercept) c(1, pred_mat[i, ]) else pred_mat[i, ]
      leave_one_out[i] <- sum(row_design * coefficients)
    }

    final <- fit_stack(pred_mat, observed)
    intercept <- if (stack_intercept) final[1] else 0
    weights <- if (stack_intercept) final[-1] else final

    if (stack_objective == "constrained") {
      weights <- project_simplex(weights)
      intercept <- 0
      leave_one_out <- as.numeric(pred_mat %*% weights)
    }

    names(weights) <- selected_ids_chr

    list(
      ensemble_pred = leave_one_out,
      weights = weights,
      intercept = intercept
    )
  }

  predictions_needed <- ensemble_scope %in% c("model", "both")

  if (predictions_needed) {
    if (is.null(predictions)) predictions <- attr(bestfit, "predictions")
    if (!is.data.frame(predictions)) {
      stop(
        "`predictions` must be supplied or stored in `bestfit`.",
        call. = FALSE
      )
    }
  }

  ensemble_summary <- NULL
  ensemble_predictions <- NULL
  ensemble_by_lag <- NULL
  prediction_weights_by_method <- list()
  lag_weights_by_method <- list()
  stack_intercepts <- setNames(rep(NA_real_, length(methods)), methods)
  metric_warning_rows <- list()
  metric_warning_index <- 1L

  collect_performance_metrics <- function(
    observed,
    predicted,
    method_name,
    hard_voting = FALSE
  ) {
    captured <- character(0)

    metrics <- withCallingHandlers(
      .compute_performance_metrics(
        observed = observed,
        predicted = predicted,
        family = family_name,
        threshold = threshold
      ),
      warning = function(w) {
        message_text <- conditionMessage(w)

        # For hard voting, 0/1 class predictions are expected by construction.
        # The generic probability-resolution warning is therefore not useful.
        if (!(hard_voting && grepl(
          "threshold.*predicted probabilities.*0/1",
          message_text,
          ignore.case = TRUE
        ))) {
          captured <<- c(captured, message_text)
        }

        invokeRestart("muffleWarning")
      }
    )

    if (hard_voting) {
      # Hard voting returns classes, not calibrated probabilities. Classification
      # metrics from the shared performance engine remain valid; probability
      # metrics are deliberately suppressed.
      metrics$ROC_AUC <- NA_real_
      metrics$Brier <- NA_real_
      metrics$LogLoss <- NA_real_
    }

    if (length(captured)) {
      captured <- unique(captured)
      for (warning_text in captured) {
        metric_warning_rows[[metric_warning_index]] <<- data.frame(
          method = method_name,
          warning = warning_text,
          stringsAsFactors = FALSE
        )
        metric_warning_index <<- metric_warning_index + 1L
      }
    }

    metrics
  }

  if (predictions_needed) {
    prediction_matrix <- build_prediction_matrix(predictions)

    summary_list <- vector("list", length(methods))
    prediction_list <- vector("list", length(methods))

    for (method_index in seq_along(methods)) {
      method_name <- methods[[method_index]]
      prediction_weights <- NULL
      stack_intercept_value <- NA_real_
      vote_fraction <- rep(NA_real_, nrow(prediction_matrix$pred_mat))

      if (method_name == "stacked") {
        stacked <- compute_stacking(
          prediction_matrix$pred_mat,
          prediction_matrix$obs
        )
        ensemble_prediction <- stacked$ensemble_pred
        prediction_weights <- stacked$weights
        stack_intercept_value <- stacked$intercept
        stack_intercepts[[method_name]] <- stack_intercept_value

        if (identical(outcome_type, "binary")) {
          tolerance <- 1e-12
          if (any(ensemble_prediction < -tolerance | ensemble_prediction > 1 + tolerance)) {
            stop(
              "`method = 'stacked'` produced values outside [0, 1] for a binomial ",
              "outcome. These cannot be interpreted as probabilities. Consider ",
              "`stack_objective = 'constrained'` or another ensemble method.",
              call. = FALSE
            )
          }
          ensemble_prediction <- pmin(pmax(ensemble_prediction, 0), 1)
        }
      } else if (method_name == "hard_voting") {
        class_matrix <- prediction_matrix$pred_mat >= threshold
        vote_fraction <- rowMeans(class_matrix)
        highest_ranked_class <- as.numeric(class_matrix[, 1L])

        ensemble_prediction <- ifelse(
          vote_fraction > 0.5,
          1,
          ifelse(
            vote_fraction < 0.5,
            0,
            highest_ranked_class
          )
        )

        prediction_weights <- setNames(
          rep(NA_real_, length(selected_ids_chr)),
          selected_ids_chr
        )
      } else {
        prediction_weights <- compute_model_weights(
          selected_models,
          method_name
        )
        ensemble_prediction <- as.numeric(
          prediction_matrix$pred_mat %*% as.numeric(prediction_weights)
        )
      }

      prediction_weights_by_method[[method_name]] <- prediction_weights

      current_predictions <- prediction_matrix$obs_df[, id_cols, drop = FALSE]
      current_predictions$method <- method_name
      current_predictions$observed <- prediction_matrix$obs
      current_predictions$predicted_ensemble <- ensemble_prediction

      if (identical(outcome_type, "binary")) {
        if (method_name == "hard_voting") {
          current_predictions$predicted_probability <- NA_real_
          current_predictions$predicted_class <- as.integer(ensemble_prediction)
          current_predictions$vote_fraction <- vote_fraction
        } else {
          current_predictions$predicted_probability <- ensemble_prediction
          current_predictions$predicted_class <- as.integer(
            ensemble_prediction >= threshold
          )
          current_predictions$vote_fraction <- NA_real_
        }
      }

      current_predictions <- cbind(
        current_predictions,
        as.data.frame(
          prediction_matrix$pred_mat,
          stringsAsFactors = FALSE
        )
      )

      prediction_list[[method_index]] <- current_predictions

      metrics <- collect_performance_metrics(
        observed = prediction_matrix$obs,
        predicted = ensemble_prediction,
        method_name = method_name,
        hard_voting = identical(method_name, "hard_voting")
      )

      summary_list[[method_index]] <- data.frame(
        method = method_name,
        top_n = length(selected_ids),
        family = family_name,
        outcome_type = outcome_type,
        weight_metric = weight_metric,
        metric_direction = metric_direction,
        threshold = if (identical(outcome_type, "binary")) threshold else NA_real_,
        stacking_model = if (method_name == "stacked") stacking_model else NA_character_,
        stack_objective = if (method_name == "stacked") stack_objective else NA_character_,
        stack_intercept = if (method_name == "stacked") stack_intercept_value else NA_real_,
        weight_transform = if (method_name == "weighted") weight_transform else NA_character_,
        metrics,
        stringsAsFactors = FALSE
      )
    }

    ensemble_summary <- do.call(rbind, summary_list)
    rownames(ensemble_summary) <- NULL

    ensemble_predictions <- do.call(rbind, prediction_list)
    rownames(ensemble_predictions) <- NULL
  }

  # -------------------------------------------------------------------------
  # Lag decomposition is generated only once for the selected models, even when
  # several compatible ensemble methods are requested. Method-specific weights
  # are applied afterwards.
  # -------------------------------------------------------------------------
  if (ensemble_scope %in% c("lag", "both") && is.null(lag_data)) {
    if (is.null(fit_list)) fit_list <- attr(bestfit, "fits")

    if (!is.list(fit_list) || !length(fit_list)) {
      stop(
        "No fitted models were found; provide `fit_list` or precomputed `lag_data`.",
        call. = FALSE
      )
    }

    if (is.null(names(fit_list)) || anyNA(names(fit_list)) ||
        any(!nzchar(names(fit_list))) || anyDuplicated(names(fit_list))) {
      stop(
        "`fit_list` must have unique non-empty model-ID names.",
        call. = FALSE
      )
    }

    if (!is.data.frame(data)) {
      stop("`data` is required for automatic lag decomposition.", call. = FALSE)
    }

    required_data <- c(group, "time")
    missing_data <- setdiff(required_data, names(data))
    if (length(missing_data)) {
      stop(
        "Missing columns in `data`: ", paste(missing_data, collapse = ", "), ".",
        call. = FALSE
      )
    }

    if (!is.numeric(data$time) || any(!is.finite(data$time))) {
      stop("`time` must contain finite numeric values.", call. = FALSE)
    }

    lag_list <- list()
    lag_index <- 1L

    for (model_id in selected_ids_chr) {
      fitted <- fit_list[[model_id]]
      if (is.null(fitted)) {
        stop("No fitted model found for model_id = ", model_id, ".", call. = FALSE)
      }

      fitted_vars <- attr(fitted, "epiexposure_vars")
      if (is.null(fitted_vars) || !length(fitted_vars)) {
        stop(
          "Model ", model_id, " lacks `epiexposure_vars` metadata.",
          call. = FALSE
        )
      }

      if (lg_strategy == "model") {
        selected_vars <- fitted_vars
      } else if (lg_strategy == "requested_available") {
        selected_vars <- if (is.null(var)) fitted_vars else intersect(var, fitted_vars)

        if (!length(selected_vars)) {
          if (verbose) {
            warning(
              "Skipping model ", model_id,
              ": no requested variables are available.",
              call. = FALSE
            )
          }
          next
        }
      } else {
        selected_vars <- if (is.null(var)) fitted_vars else var
        missing_vars <- setdiff(selected_vars, fitted_vars)
        if (length(missing_vars)) {
          stop(
            "Model ", model_id, " lacks variables: ",
            paste(missing_vars, collapse = ", "), ".",
            call. = FALSE
          )
        }
      }

      call_args <- c(
        list(
          data = data,
          group = group,
          fit = fitted,
          var = selected_vars
        ),
        compute_ecilag_args
      )

      decomposition <- tryCatch(
        do.call(compute_ecilag, call_args),
        error = function(e) {
          stop(
            "Lag decomposition failed for model ", model_id,
            ": ", conditionMessage(e),
            call. = FALSE
          )
        }
      )

      if (!is.list(decomposition) || !is.data.frame(decomposition$by_lag)) {
        stop(
          "`compute_ecilag()` returned no valid `by_lag` for model ",
          model_id, ".",
          call. = FALSE
        )
      }

      temporary <- decomposition$by_lag
      temporary[[model_col]] <- type.convert(model_id, as.is = TRUE)
      lag_list[[lag_index]] <- temporary
      lag_index <- lag_index + 1L
    }

    if (!length(lag_list)) {
      stop("No lag decompositions were generated.", call. = FALSE)
    }

    lag_data <- do.call(rbind, lag_list)
    rownames(lag_data) <- NULL

    if (verbose) message("Lag data generated from stored fitted models.")
  }

  if (ensemble_scope %in% c("lag", "both")) {
    if (!is.data.frame(lag_data)) {
      stop("`lag_data` must be a data.frame.", call. = FALSE)
    }

    required_lag <- unique(
      c(model_col, lag_group_cols, lag_col, lag_contribution_col)
    )
    missing_lag <- setdiff(required_lag, names(lag_data))

    if (length(missing_lag)) {
      stop(
        "Missing columns in `lag_data`: ",
        paste(missing_lag, collapse = ", "), ".",
        call. = FALSE
      )
    }

    if (!is.numeric(lag_data[[lag_contribution_col]]) ||
        any(!is.finite(lag_data[[lag_contribution_col]]))) {
      stop("Lag contributions must be finite numeric values.", call. = FALSE)
    }

    lag_use <- lag_data[
      as.character(lag_data[[model_col]]) %in% selected_ids_chr,
      ,
      drop = FALSE
    ]

    if (!nrow(lag_use)) {
      stop("No lag contributions were found for selected models.", call. = FALSE)
    }

    lag_key_cols <- unique(c(model_col, lag_group_cols))
    if (anyDuplicated(lag_use[lag_key_cols])) {
      stop(
        "`lag_data` has multiple rows for the same model and lag unit. ",
        "Include the relevant grouping column in `lag_group_cols`, for example ",
        "c('epi_id', 'var', 'lag'), or aggregate rows before calling this function.",
        call. = FALSE
      )
    }

    make_group_key <- function(d, columns) {
      if (length(columns) == 1L) return(as.character(d[[columns]]))
      apply(d[, columns, drop = FALSE], 1L, paste, collapse = "||")
    }

    lag_method_rows <- vector("list", length(methods))

    for (method_index in seq_along(methods)) {
      method_name <- methods[[method_index]]
      lag_weights <- compute_model_weights(selected_models, method_name)
      lag_weights_by_method[[method_name]] <- lag_weights

      current_lag_use <- lag_use
      current_lag_use$.model_weight <- as.numeric(
        lag_weights[
          match(
            as.character(current_lag_use[[model_col]]),
            names(lag_weights)
          )
        ]
      )

      if (any(!is.finite(current_lag_use$.model_weight))) {
        stop("Could not match model weights to `lag_data`.", call. = FALSE)
      }

      group_key <- make_group_key(current_lag_use, lag_group_cols)

      ensemble_rows <- lapply(
        unique(group_key),
        function(current_key) {
          temporary <- current_lag_use[group_key == current_key, , drop = FALSE]
          weights <- temporary$.model_weight

          if (rl_weights) {
            total_weight <- sum(weights)
            if (total_weight > 0) weights <- weights / total_weight
          }

          output <- temporary[1, lag_group_cols, drop = FALSE]
          output$method <- method_name
          output$ensemble_contribution <- sum(
            temporary[[lag_contribution_col]] * weights
          )

          if (lag_weight_col %in% names(temporary)) {
            output$ensemble_weight <- sum(
              temporary[[lag_weight_col]] * weights
            )
          }

          output$n_models <- nrow(temporary)
          output$models_used <- paste(
            sort(unique(temporary[[model_col]])),
            collapse = ", "
          )
          output
        }
      )

      method_lag <- do.call(rbind, ensemble_rows)
      rownames(method_lag) <- NULL

      percent_groups <- setdiff(lag_group_cols, lag_col)

      if (!length(percent_groups)) {
        denominator <- sum(abs(method_lag$ensemble_contribution))
        method_lag$ensemble_percent_contribution <- if (denominator > 0) {
          100 * abs(method_lag$ensemble_contribution) / denominator
        } else {
          NA_real_
        }
      } else {
        percent_key <- make_group_key(method_lag, percent_groups)
        method_lag$ensemble_percent_contribution <- NA_real_

        for (current_key in unique(percent_key)) {
          rows <- which(percent_key == current_key)
          denominator <- sum(abs(method_lag$ensemble_contribution[rows]))
          method_lag$ensemble_percent_contribution[rows] <- if (denominator > 0) {
            100 * abs(method_lag$ensemble_contribution[rows]) / denominator
          } else {
            NA_real_
          }
        }
      }

      lag_method_rows[[method_index]] <- method_lag
    }

    ensemble_by_lag <- do.call(rbind, lag_method_rows)
    rownames(ensemble_by_lag) <- NULL
  }

  # -------------------------------------------------------------------------
  # Model weights are method-specific. Keeping a method column makes a single
  # multi-method call directly comparable without separate function runs.
  # -------------------------------------------------------------------------
  model_weight_rows <- vector("list", length(methods))

  for (method_index in seq_along(methods)) {
    method_name <- methods[[method_index]]
    prediction_weights <- prediction_weights_by_method[[method_name]]
    lag_weights <- lag_weights_by_method[[method_name]]

    current <- data.frame(
      method = method_name,
      selected_ids,
      prediction_weight = if (!is.null(prediction_weights)) {
        as.numeric(
          prediction_weights[
            match(selected_ids_chr, names(prediction_weights))
          ]
        )
      } else {
        NA_real_
      },
      lag_weight = if (!is.null(lag_weights)) {
        as.numeric(
          lag_weights[
            match(selected_ids_chr, names(lag_weights))
          ]
        )
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )

    names(current)[2] <- model_col
    model_weight_rows[[method_index]] <- current
  }

  model_weights <- do.call(rbind, model_weight_rows)
  rownames(model_weights) <- NULL

  # Preserve the legacy single-stacked-method attribute when it is
  # unambiguous, while the summary table remains the authoritative location for
  # multi-method stack intercepts.
  if (length(methods) == 1L && identical(methods, "stacked")) {
    attr(model_weights, "stack_intercept") <- stack_intercepts[["stacked"]]
  }

  metric_warnings <- if (length(metric_warning_rows)) {
    output <- do.call(rbind, metric_warning_rows)
    rownames(output) <- NULL
    unique(output)
  } else {
    data.frame(
      method = character(0),
      warning = character(0),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(metric_warnings)) {
    warning(
      "One or more ensemble performance metrics were undefined or required ",
      "special handling. See `result$metric_warnings` for details.",
      call. = FALSE
    )
  }

  out <- list(
    ensemble_summary = ensemble_summary,
    ensemble_predictions = ensemble_predictions,
    ensemble_by_lag = ensemble_by_lag,
    lag_data = lag_data,
    model_weights = model_weights,
    selected_models = selected_models,
    metric_warnings = metric_warnings,
    method = methods,
    ensemble_scope = ensemble_scope,
    lg_strategy = lg_strategy,
    rl_weights = rl_weights,
    family = family_name,
    outcome_type = outcome_type,
    weight_metric = weight_metric,
    metric_direction = metric_direction,
    threshold = threshold
  )

  out
}
