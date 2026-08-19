#' Create prediction and lag-contribution ensembles from best-fit models
#'
#' Builds model-level ensembles from out-of-fold predictions produced by
#' `find_bestfit()` and, optionally, lag-level ensembles from lag-specific
#' contributions. Automatic lag decomposition uses `compute_ecilag()`.
#'
#' Exposure histories supplied through `data` must be chronological, from the
#' earliest to the most recent observation within each group. Their conversion
#' to the retrospective lag representation is handled by `compute_ecilag()`.
#'
#' @param bestfit Data frame returned by `find_bestfit()`. It must contain the
#'   model identifier and the metric selected by `weight_metric`.
#' @param predictions Optional data frame of out-of-fold predictions. If `NULL`,
#'   `attr(bestfit, "predictions")` is used. Required columns are `model_col`,
#'   `id_cols`, `observed`, and `predicted`.
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
#' @param ensemble_scope Character. `"model"`, `"lag"`, or `"both"`.
#' @param method Character. `"best"`, `"hard_voting"`, `"unweighted"`,
#'   `"weighted"`, or `"stacked"`. Stacking is available only for model-level
#'   ensembles.
#' @param top_n Positive integer number of top-ranked models.
#' @param model_ids Optional vector of model identifiers. If supplied, `top_n`
#'   is ignored, but selected models are still ordered by `weight_metric`.
#' @param weight_metric Numeric column in `bestfit` used for ranking and weights.
#' @param weight_transform Character. `"softmax"`, `"positive"`,
#'   `"rank_inverse"`, or `"uniform"`.
#' @param stacking_model Character. `"ridge"` or `"lm"`.
#' @param stack_objective Character. `"regularized"`, `"min_error"`, or
#'   `"constrained"`.
#' @param lambda Non-negative ridge penalty.
#' @param stack_intercept Logical. Include an intercept in stacking.
#' @param model_col Character scalar naming the model identifier column.
#' @param id_cols Character vector identifying prediction rows.
#' @param lag_group_cols Character vector defining each lag-ensemble unit.
#'   Include an epidemic/group column, for example
#'   `c("epi_id", "var", "lag")`, when lag contributions are group-specific.
#' @param lg_strategy Character. `"requested_available"`, `"model"`, or
#'   `"strict"`.
#' @param rl_weights Logical. Renormalize model weights within each lag unit.
#' @param hb Logical. If `TRUE`, higher `weight_metric` values are better.
#' @param verbose Logical. Print informative messages.
#'
#' @return A list containing `ensemble_summary`, `ensemble_predictions`,
#'   `ensemble_by_lag`, `lag_data`, `model_weights`, `selected_models`,
#'   `method`, `ensemble_scope`, `lg_strategy`, and `rl_weights`.
#'
#' @details
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
    method = c("unweighted", "weighted", "hard_voting", "best", "stacked"),
    top_n = 3,
    model_ids = NULL,
    weight_metric = "CCC",
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
    hb = TRUE,
    verbose = TRUE
) {
  ensemble_scope <- match.arg(ensemble_scope)
  method <- match.arg(method)
  weight_transform <- match.arg(weight_transform)
  stacking_model <- match.arg(stacking_model)
  stack_objective <- match.arg(stack_objective)
  lg_strategy <- match.arg(lg_strategy)

  response_col <- "observed"
  prediction_col <- "predicted"
  lag_col <- "lag"
  lag_contribution_col <- "contribution"
  lag_weight_col <- "weight"

  valid_name <- function(x) {
    is.character(x) && length(x) >= 1L && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
  }
  valid_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)

  if (!is.data.frame(bestfit)) stop("`bestfit` must be a data.frame.")
  if (!valid_name(model_col) || length(model_col) != 1L) {
    stop("`model_col` must be one non-empty column name.")
  }
  if (!valid_name(id_cols)) stop("`id_cols` must contain unique non-empty names.")
  if (!valid_name(lag_group_cols)) {
    stop("`lag_group_cols` must contain unique non-empty names.")
  }
  if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(group)) {
    stop("`group` must be one non-empty column name.")
  }
  if (!is.null(var) && (!is.character(var) || !length(var) || anyNA(var) ||
                        any(!nzchar(var)) || anyDuplicated(var))) {
    stop("`var` must be NULL or a unique non-empty character vector.")
  }
  if (!model_col %in% names(bestfit)) {
    stop("`model_col` ('", model_col, "') was not found in `bestfit`.")
  }
  if (!weight_metric %in% names(bestfit)) {
    stop("`weight_metric` ('", weight_metric, "') was not found in `bestfit`.")
  }
  if (!is.numeric(bestfit[[weight_metric]])) {
    stop("`weight_metric` must identify a numeric column in `bestfit`.")
  }
  if (anyDuplicated(bestfit[[model_col]])) {
    stop("`model_col` must uniquely identify rows in `bestfit`.")
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
    stop("`top_n` must be a positive integer.")
  }
  top_n <- as.integer(top_n)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("`lambda` must be a non-negative numeric scalar.")
  }
  if (!valid_flag(stack_intercept)) stop("`stack_intercept` must be TRUE or FALSE.")
  if (!valid_flag(rl_weights)) stop("`rl_weights` must be TRUE or FALSE.")
  if (!valid_flag(hb)) stop("`hb` must be TRUE or FALSE.")
  if (!valid_flag(verbose)) stop("`verbose` must be TRUE or FALSE.")
  if (!is.list(compute_ecilag_args)) stop("`compute_ecilag_args` must be a named list.")
  if (length(compute_ecilag_args) &&
      (is.null(names(compute_ecilag_args)) || anyNA(names(compute_ecilag_args)) ||
       any(!nzchar(names(compute_ecilag_args))) || anyDuplicated(names(compute_ecilag_args)))) {
    stop("`compute_ecilag_args` must have unique non-empty names.")
  }
  reserved_args <- c("data", "group", "fit", "var", "profile")
  forbidden_args <- intersect(names(compute_ecilag_args), reserved_args)
  if (length(forbidden_args)) {
    stop("`compute_ecilag_args` cannot override: ",
         paste(forbidden_args, collapse = ", "), ".")
  }
  if (!is.null(model_ids) && anyDuplicated(model_ids)) {
    stop("`model_ids` must contain unique identifiers.")
  }
  if (ensemble_scope %in% c("lag", "both") && method == "stacked") {
    stop("`method = 'stacked'` is supported only for `ensemble_scope = 'model'`.")
  }

  ccc_lins <- function(obs, pred) {
    ok <- is.finite(obs) & is.finite(pred)
    obs <- obs[ok]
    pred <- pred[ok]
    if (length(obs) < 2L) {
      return(data.frame(CCC = NA_real_, Cb = NA_real_, rho = NA_real_,
                        RMSE = NA_real_, MAE = NA_real_))
    }
    mx <- mean(obs)
    my <- mean(pred)
    vx <- stats::var(obs)
    vy <- stats::var(pred)
    sxy <- stats::cov(obs, pred)
    rho <- suppressWarnings(stats::cor(obs, pred))
    ccc <- (2 * sxy) / (vx + vy + (mx - my)^2)
    cb <- if (is.finite(rho) && abs(rho) > .Machine$double.eps) ccc / rho else NA_real_
    data.frame(
      CCC = as.numeric(ccc),
      Cb = as.numeric(cb),
      rho = as.numeric(rho),
      RMSE = sqrt(mean((obs - pred)^2)),
      MAE = mean(abs(obs - pred))
    )
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
    score <- if (hb) metric else -metric
    model_table[order(-score, na.last = NA), , drop = FALSE]
  }

  if (!is.null(model_ids)) {
    missing_models <- setdiff(model_ids, bestfit[[model_col]])
    if (length(missing_models)) {
      stop("The following `model_ids` were not found: ",
           paste(missing_models, collapse = ", "), ".")
    }
    selected_models <- bestfit[bestfit[[model_col]] %in% model_ids, , drop = FALSE]
    selected_models <- rank_models(selected_models)
  } else {
    selected_models <- utils::head(rank_models(bestfit), top_n)
  }

  if (!nrow(selected_models)) stop("No models were selected for the ensemble.")
  if (all(!is.finite(selected_models[[weight_metric]]))) {
    stop("No finite `weight_metric` values are available for selected models.")
  }
  selected_ids <- selected_models[[model_col]]
  selected_ids_chr <- as.character(selected_ids)

  if (verbose) {
    message("Selected ", length(selected_ids), " model(s): ",
            paste(selected_ids, collapse = ", "))
  }

  compute_model_weights <- function(model_table, ensemble_method) {
    n_models <- nrow(model_table)
    ids <- as.character(model_table[[model_col]])
    if (ensemble_method == "unweighted") {
      weights <- rep(1 / n_models, n_models)
    } else if (ensemble_method %in% c("best", "hard_voting")) {
      weights <- c(1, rep(0, n_models - 1L))
    } else {
      metric <- model_table[[weight_metric]]
      score <- if (hb) metric else -metric
      finite <- is.finite(score)
      if (!any(finite)) stop("No finite values are available for weighting.")
      if (weight_transform == "uniform") {
        weights <- rep(1 / n_models, n_models)
      } else if (weight_transform == "rank_inverse") {
        weights <- 1 / seq_len(n_models)
        weights <- weights / sum(weights)
      } else if (weight_transform == "positive") {
        positive <- pmax(score, 0)
        positive[!is.finite(positive)] <- 0
        weights <- if (sum(positive) > 0) positive / sum(positive) else rep(1 / n_models, n_models)
      } else {
        safe_score <- score
        safe_score[!finite] <- -Inf
        shifted <- safe_score - max(safe_score[finite])
        exp_score <- exp(shifted)
        exp_score[!is.finite(exp_score)] <- 0
        weights <- exp_score / sum(exp_score)
      }
    }
    names(weights) <- ids
    weights
  }

  build_prediction_matrix <- function(prediction_data) {
    required <- c(model_col, id_cols, response_col, prediction_col)
    missing <- setdiff(required, names(prediction_data))
    if (length(missing)) {
      stop("Missing prediction columns: ", paste(missing, collapse = ", "), ".")
    }
    prediction_data <- prediction_data[
      as.character(prediction_data[[model_col]]) %in% selected_ids_chr,
      ,
      drop = FALSE
    ]
    if (!nrow(prediction_data)) stop("No predictions were found for selected models.")
    if (!is.numeric(prediction_data[[response_col]]) ||
        !is.numeric(prediction_data[[prediction_col]])) {
      stop("`observed` and `predicted` must be numeric.")
    }
    key_cols <- c(model_col, id_cols)
    if (anyDuplicated(prediction_data[key_cols])) {
      stop("`predictions` contains duplicate rows for a model and prediction identifier.")
    }
    make_key <- function(d, columns) {
      if (length(columns) == 1L) return(as.character(d[[columns]]))
      apply(d[, columns, drop = FALSE], 1L, paste, collapse = "||")
    }
    prediction_data$.obs_key <- make_key(prediction_data, id_cols)
    observed_sets <- split(prediction_data[[response_col]], prediction_data$.obs_key)
    inconsistent <- vapply(observed_sets, function(z) {
      length(unique(z[is.finite(z)])) > 1L
    }, logical(1))
    if (any(inconsistent)) {
      stop("Different observed values were found for the same prediction identifier.")
    }
    obs_df <- prediction_data[, c(".obs_key", id_cols, response_col), drop = FALSE]
    obs_df <- obs_df[!duplicated(obs_df$.obs_key), , drop = FALSE]
    obs_keys <- obs_df$.obs_key
    pred_mat <- matrix(
      NA_real_, nrow = length(obs_keys), ncol = length(selected_ids),
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
      stop("Fewer than two complete OOF observations are shared by selected models.")
    }
    list(
      pred_mat = pred_mat[complete, , drop = FALSE],
      obs_df = obs_df[complete, , drop = FALSE],
      obs = obs_df[[response_col]][complete]
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
        as.numeric(solve(crossprod(design) + penalty, crossprod(design, y))),
        error = function(e) rep(0, ncol(design))
      )
    }
    leave_one_out <- numeric(nrow(pred_mat))
    for (i in seq_len(nrow(pred_mat))) {
      train <- setdiff(seq_len(nrow(pred_mat)), i)
      coefficients <- fit_stack(pred_mat[train, , drop = FALSE], observed[train])
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
    list(ensemble_pred = leave_one_out, weights = weights, intercept = intercept)
  }

  predictions_needed <- ensemble_scope %in% c("model", "both")
  if (predictions_needed) {
    if (is.null(predictions)) predictions <- attr(bestfit, "predictions")
    if (!is.data.frame(predictions)) {
      stop("`predictions` must be supplied or stored in `bestfit`.")
    }
  }

  ensemble_summary <- NULL
  ensemble_predictions <- NULL
  prediction_weights <- NULL
  lag_weights <- NULL
  stack_intercept_value <- NULL

  if (predictions_needed) {
    prediction_matrix <- build_prediction_matrix(predictions)
    if (method == "stacked") {
      stacked <- compute_stacking(prediction_matrix$pred_mat, prediction_matrix$obs)
      ensemble_prediction <- stacked$ensemble_pred
      prediction_weights <- stacked$weights
      stack_intercept_value <- stacked$intercept
    } else {
      prediction_weights <- compute_model_weights(selected_models, method)
      ensemble_prediction <- as.numeric(
        prediction_matrix$pred_mat %*% as.numeric(prediction_weights)
      )
    }
    ensemble_predictions <- prediction_matrix$obs_df[, id_cols, drop = FALSE]
    ensemble_predictions$observed <- prediction_matrix$obs
    ensemble_predictions$predicted_ensemble <- ensemble_prediction
    ensemble_predictions <- cbind(
      ensemble_predictions,
      as.data.frame(prediction_matrix$pred_mat, stringsAsFactors = FALSE)
    )
    metrics <- ccc_lins(ensemble_predictions$observed,
                        ensemble_predictions$predicted_ensemble)
    ensemble_summary <- data.frame(
      method = method,
      top_n = length(selected_ids),
      stacking_model = if (method == "stacked") stacking_model else NA_character_,
      stack_objective = if (method == "stacked") stack_objective else NA_character_,
      weight_metric = if (method == "weighted") weight_metric else NA_character_,
      weight_transform = if (method == "weighted") weight_transform else NA_character_,
      metrics,
      stringsAsFactors = FALSE
    )
  }

  if (ensemble_scope %in% c("lag", "both")) {
    lag_weights <- compute_model_weights(selected_models, method)
  }

  if (ensemble_scope %in% c("lag", "both") && is.null(lag_data)) {
    if (is.null(fit_list)) fit_list <- attr(bestfit, "fits")
    if (!is.list(fit_list) || !length(fit_list)) {
      stop("No fitted models were found; provide `fit_list` or precomputed `lag_data`.")
    }
    if (is.null(names(fit_list)) || anyNA(names(fit_list)) ||
        any(!nzchar(names(fit_list))) || anyDuplicated(names(fit_list))) {
      stop("`fit_list` must have unique non-empty model-ID names.")
    }
    if (!is.data.frame(data)) {
      stop("`data` is required for automatic lag decomposition.")
    }
    required_data <- c(group, "time")
    missing_data <- setdiff(required_data, names(data))
    if (length(missing_data)) {
      stop("Missing columns in `data`: ", paste(missing_data, collapse = ", "), ".")
    }
    if (!is.numeric(data$time) || any(!is.finite(data$time))) {
      stop("`time` must contain finite numeric values.")
    }

    lag_list <- list()
    lag_index <- 1L
    for (model_id in selected_ids_chr) {
      fitted <- fit_list[[model_id]]
      if (is.null(fitted)) stop("No fitted model found for model_id = ", model_id, ".")
      fitted_vars <- attr(fitted, "epiexposure_vars")
      if (is.null(fitted_vars) || !length(fitted_vars)) {
        stop("Model ", model_id, " lacks `epiexposure_vars` metadata.")
      }
      if (lg_strategy == "model") {
        selected_vars <- fitted_vars
      } else if (lg_strategy == "requested_available") {
        selected_vars <- if (is.null(var)) fitted_vars else intersect(var, fitted_vars)
        if (!length(selected_vars)) {
          if (verbose) warning("Skipping model ", model_id,
                               ": no requested variables are available.", call. = FALSE)
          next
        }
      } else {
        selected_vars <- if (is.null(var)) fitted_vars else var
        missing_vars <- setdiff(selected_vars, fitted_vars)
        if (length(missing_vars)) {
          stop("Model ", model_id, " lacks variables: ",
               paste(missing_vars, collapse = ", "), ".")
        }
      }
      call_args <- c(
        list(data = data, group = group, fit = fitted, var = selected_vars),
        compute_ecilag_args
      )
      decomposition <- tryCatch(
        do.call(compute_ecilag, call_args),
        error = function(e) stop("Lag decomposition failed for model ", model_id,
                                 ": ", conditionMessage(e), call. = FALSE)
      )
      if (!is.list(decomposition) || !is.data.frame(decomposition$by_lag)) {
        stop("`compute_ecilag()` returned no valid `by_lag` for model ", model_id, ".")
      }
      temporary <- decomposition$by_lag
      temporary[[model_col]] <- type.convert(model_id, as.is = TRUE)
      lag_list[[lag_index]] <- temporary
      lag_index <- lag_index + 1L
    }
    if (!length(lag_list)) stop("No lag decompositions were generated.")
    lag_data <- do.call(rbind, lag_list)
    rownames(lag_data) <- NULL
    if (verbose) message("Lag data generated from stored fitted models.")
  }

  ensemble_by_lag <- NULL
  if (ensemble_scope %in% c("lag", "both")) {
    if (!is.data.frame(lag_data)) stop("`lag_data` must be a data.frame.")
    required_lag <- unique(c(model_col, lag_group_cols, lag_col, lag_contribution_col))
    missing_lag <- setdiff(required_lag, names(lag_data))
    if (length(missing_lag)) {
      stop("Missing columns in `lag_data`: ", paste(missing_lag, collapse = ", "), ".")
    }
    if (!is.numeric(lag_data[[lag_contribution_col]]) ||
        any(!is.finite(lag_data[[lag_contribution_col]]))) {
      stop("Lag contributions must be finite numeric values.")
    }
    lag_use <- lag_data[
      as.character(lag_data[[model_col]]) %in% selected_ids_chr,
      ,
      drop = FALSE
    ]
    if (!nrow(lag_use)) stop("No lag contributions were found for selected models.")
    lag_key_cols <- unique(c(model_col, lag_group_cols))
    if (anyDuplicated(lag_use[lag_key_cols])) {
      stop(
        "`lag_data` has multiple rows for the same model and lag unit. ",
        "Include the relevant grouping column in `lag_group_cols`, for example ",
        "c('epi_id', 'var', 'lag'), or aggregate rows before calling this function."
      )
    }
    lag_use$.model_weight <- as.numeric(
      lag_weights[match(as.character(lag_use[[model_col]]), names(lag_weights))]
    )
    if (any(!is.finite(lag_use$.model_weight))) {
      stop("Could not match model weights to `lag_data`.")
    }
    make_group_key <- function(d, columns) {
      if (length(columns) == 1L) return(as.character(d[[columns]]))
      apply(d[, columns, drop = FALSE], 1L, paste, collapse = "||")
    }
    group_key <- make_group_key(lag_use, lag_group_cols)
    ensemble_rows <- lapply(unique(group_key), function(current_key) {
      temporary <- lag_use[group_key == current_key, , drop = FALSE]
      weights <- temporary$.model_weight
      if (rl_weights) {
        total_weight <- sum(weights)
        if (total_weight > 0) weights <- weights / total_weight
      }
      output <- temporary[1, lag_group_cols, drop = FALSE]
      output$ensemble_contribution <- sum(
        temporary[[lag_contribution_col]] * weights
      )
      if (lag_weight_col %in% names(temporary)) {
        output$ensemble_weight <- sum(temporary[[lag_weight_col]] * weights)
      }
      output$n_models <- nrow(temporary)
      output$models_used <- paste(sort(unique(temporary[[model_col]])), collapse = ", ")
      output
    })
    ensemble_by_lag <- do.call(rbind, ensemble_rows)
    rownames(ensemble_by_lag) <- NULL

    percent_groups <- setdiff(lag_group_cols, lag_col)
    if (!length(percent_groups)) {
      denominator <- sum(abs(ensemble_by_lag$ensemble_contribution))
      ensemble_by_lag$ensemble_percent_contribution <- if (denominator > 0) {
        100 * abs(ensemble_by_lag$ensemble_contribution) / denominator
      } else NA_real_
    } else {
      percent_key <- make_group_key(ensemble_by_lag, percent_groups)
      ensemble_by_lag$ensemble_percent_contribution <- NA_real_
      for (current_key in unique(percent_key)) {
        rows <- which(percent_key == current_key)
        denominator <- sum(abs(ensemble_by_lag$ensemble_contribution[rows]))
        ensemble_by_lag$ensemble_percent_contribution[rows] <- if (denominator > 0) {
          100 * abs(ensemble_by_lag$ensemble_contribution[rows]) / denominator
        } else NA_real_
      }
    }
  }

  model_weights <- data.frame(
    selected_ids,
    prediction_weight = if (!is.null(prediction_weights)) {
      as.numeric(prediction_weights[match(selected_ids_chr, names(prediction_weights))])
    } else NA_real_,
    lag_weight = if (!is.null(lag_weights)) {
      as.numeric(lag_weights[match(selected_ids_chr, names(lag_weights))])
    } else NA_real_,
    stringsAsFactors = FALSE
  )
  names(model_weights)[1] <- model_col
  if (!is.null(stack_intercept_value)) {
    attr(model_weights, "stack_intercept") <- stack_intercept_value
  }

  out <- list(
    ensemble_summary = ensemble_summary,
    ensemble_predictions = ensemble_predictions,
    ensemble_by_lag = ensemble_by_lag,
    lag_data = lag_data,
    model_weights = model_weights,
    selected_models = selected_models,
    method = method,
    ensemble_scope = ensemble_scope,
    lg_strategy = lg_strategy,
    rl_weights = rl_weights
  )

  out
}
