#' Find the best DLNM model structure using LOOCV
#'
#' Tests combinations of exposure variables and basis dimensions using
#' leave-one-group-out cross-validation and family-aware predictive performance
#' metrics. Candidate models are ranked by a user-selected performance metric,
#' with metric direction (maximize or minimize) determined automatically.
#'
#' Input observations must be chronological within each group, from the earliest
#' to the most recent observation. The function orders rows internally by
#' `group` and `time`. The DLNM lag representation is constructed internally by
#' `define_exposure()` and `build_design()`.
#'
#' @param data Long-format data frame containing the response, grouping column,
#'   chronological time column, and candidate exposure variables.
#' @param response Character scalar naming the response column in `data`. The
#'   response may be repeated across rows but must be constant within each group.
#' @param group Character scalar naming the grouping column used for LOOCV.
#' @param time Character scalar naming the chronological time column.
#' @param vars Unique character vector of candidate exposure variables.
#' @param lag_max Maximum lag. A numeric vector such as `c(0, 85)` is normalized
#'   to its finite maximum.
#' @param df_var_grid Positive finite candidate dimensions for the exposure basis.
#' @param df_lag_grid Positive finite candidate dimensions for the lag basis.
#' @param min_vars Positive integer minimum number of exposure variables.
#' @param max_vars Positive integer maximum number of exposure variables, or
#'   `NULL` to use `length(vars)`.
#' @param var_sets Optional list of exact variable combinations. When supplied,
#'   `min_vars` and `max_vars` are ignored.
#' @param fun_var Character exposure-basis function passed to `define_exposure()`.
#' @param fun_lag Character lag-basis function passed to `define_exposure()`.
#' @param model_engine Character model engine supported by `fit_epidlnm()`:
#'   `"glm"`, `"glmmTMB"`, `"gam"`, `"gamm"`, `"gls"`, `"spamm"`,
#'   `"brms"`, `"inla"`, or `"bdlnm"`.
#' @param family Character or family object passed to `prepare_response()` and
#'   `fit_epidlnm()`. The same family is used to determine the appropriate
#'   performance metrics. Supported canonical families are `"beta"`,
#'   `"binomial"`, `"poisson"`, `"gamma"`, `"gaussian"`,
#'   `"negative_binomial"`, and `"ordinal"`; aliases are normalized by the
#'   shared internal family resolver.
#' @param random_effect Optional character scalar naming a random-effect column.
#' @param min_success Minimum number of successful LOOCV folds required.
#' @param rank_metric Optional character scalar selecting the primary metric used
#'   to rank candidate models. If `NULL`, the default is `"CCC"` for non-binary
#'   outcomes and `"ROC_AUC"` for binary outcomes. For non-binary outcomes the
#'   available metrics are `"CCC"`, `"Cb"`, `"rho"`, `"RMSE"`, and `"MAE"`.
#'   For binary outcomes the available metrics are `"ROC_AUC"`, `"Brier"`,
#'   `"LogLoss"`, `"Accuracy"`, `"Balanced_Accuracy"`, `"Sensitivity"`,
#'   `"Specificity"`, `"F1"`, `"MCC"`, and `"Precision"`. The function
#'   automatically maximizes or minimizes the selected metric as appropriate.
#' @param threshold Numeric probability threshold strictly between 0 and 1 used
#'   only for `family = "binomial"` to convert predicted probabilities into
#'   classes 0/1 for classification metrics. The default is `0.5`. The threshold
#'   applies to predicted probabilities, not predictions that have already been
#'   reduced to class labels 0/1. It affects `Accuracy`, `Balanced_Accuracy`,
#'   `Sensitivity`, `Specificity`, `F1`, `MCC`, and `Precision`, but does not
#'   affect `ROC_AUC`, `Brier`, or `LogLoss`. It is ignored for non-binary
#'   outcomes.
#' @param top_n Positive integer or `Inf`; number of ranked candidates returned.
#' @param keep_fits Logical. Refit successful retained candidates on all data and
#'   store them in `attr(result, "fits")`.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments passed to `fit_epidlnm()`.
#'
#' @return A data frame ranked by `rank_metric`. Non-binary families report
#'   `CCC`, `Cb`, `rho`, `RMSE`, and `MAE`. Binary families report `ROC_AUC`,
#'   `Brier`, `LogLoss`, `Accuracy`, `Balanced_Accuracy`, `Sensitivity`,
#'   `Specificity`, `F1`, `MCC`, and `Precision`.
#'
#'   In addition to the original `"predictions"` and `"failures"` attributes,
#'   the returned object stores `"predictions_by_model"` (a named list of LOOCV
#'   prediction data frames), `"warnings"` (warning details and classes by
#'   candidate/fold/stage), `"warning_summary"` (one warning summary per
#'   candidate, including counts for knot, convergence, Hessian, other, and
#'   performance-metric warnings), `"timing"` (overall run timing), and
#'   `"candidate_times"` (elapsed time per candidate). If `keep_fits = TRUE`,
#'   fitted full-data models are stored in `"fits"`.
#'
#'   The result also stores standardized metadata in attributes `"family"`,
#'   `"outcome_type"`, `"rank_metric"`, and `"threshold"` for downstream use,
#'   including by `ensemble_bestfit()`.
#'
#' @details
#' For each fold, basis templates are estimated from the training groups only
#' and then applied to the held-out group. This prevents the held-out exposure
#' history from determining the training basis specification. The LOOCV fold
#' construction, basis estimation, model fitting, and held-out prediction steps
#' are not altered by the performance-metric layer.
#'
#' Performance metrics are computed only after all successful out-of-fold
#' predictions for a candidate have been collected. Both `find_bestfit()` and
#' `ensemble_bestfit()` use the shared internal `.compute_performance_metrics()`
#' implementation so the same metric definitions are used throughout the package.
#'
#' For non-binary outcomes (`beta`, `poisson`, `gamma`, `gaussian`, and
#' `negative_binomial`), performance is summarized by CCC, Cb, Pearson
#' correlation, RMSE, and MAE. For binary outcomes (`binomial`), predicted
#' probabilities are retained and evaluated using ROC AUC, Brier score, and Log
#' Loss; the configured `threshold` is then applied to those probabilities to
#' obtain 0/1 classes for Accuracy, Balanced Accuracy, Sensitivity, Specificity,
#' F1, MCC, and Precision.
#'
#' If a binary outcome contains only one observed class, metrics that require
#' both classes are returned as `NA` and a diagnostic warning is recorded. A
#' threshold is meaningful for predicted probabilities; if binary predictions
#' contain only 0/1 values, a warning is recorded because probability-resolution
#' information has already been lost.
#'
#' The `ordinal` family is recognized and normalized consistently with
#' `fit_epidlnm()`, but ordinal-specific performance metrics have not yet been
#' defined. `find_bestfit()` therefore stops explicitly for ordinal outcomes
#' rather than silently applying continuous or binary metrics.
#'
#' Candidate ranking uses `rank_metric` as the primary criterion. Remaining
#' family-appropriate metrics are used only as deterministic tie-breakers in
#' their canonical order, each with its own automatically recognized direction.
#'
#' Direct engine-specific prediction is attempted first. If unavailable, a
#' fixed-effect fallback uses the official cross-basis column order stored in
#' `epiexposure_cb_cols`, or `epiexposure_data_template` as a fallback. Missing
#' model or design columns produce an explicit error rather than being omitted.
#'
#' The function ranks point predictions and does not propagate coefficient or
#' posterior-draw uncertainty.
#'
#' Warnings raised while fitting or predicting individual LOOCV folds are
#' captured internally and muffled. They therefore do not flood the console.
#' The number of affected folds and warning events is recorded for each
#' candidate. Fit/prediction warnings are additionally classified as `"knot"`,
#' `"convergence"`, `"hessian"`, or `"other"`, without changing model fitting,
#' eligibility, or ranking. Performance-metric warnings are stored separately
#' at the candidate level and do not count as failed LOOCV folds. When
#' `verbose = TRUE`, at most one warning summary message is emitted per affected
#' candidate after candidate evaluation finishes.
#'
#' Timing is diagnostic only and does not alter model fitting. Total elapsed
#' time and per-candidate elapsed times are stored as attributes of the result.
#'
#' @export
find_bestfit <- function(
    data,
    response = "y",
    group = "epi_id",
    time = "time",
    vars,
    lag_max,
    df_var_grid = c(3, 4, 5),
    df_lag_grid = c(3, 4, 5),
    min_vars = 1,
    max_vars = NULL,
    var_sets = NULL,
    fun_var = "ns",
    fun_lag = "ns",
    model_engine = "glmmTMB",
    family = "beta",
    random_effect = NULL,
    min_success = 2,
    rank_metric = NULL,
    threshold = 0.5,
    top_n = Inf,
    keep_fits = FALSE,
    verbose = TRUE,
    ...
) {
  model_engine <- match.arg(
    model_engine,
    choices = c("glm", "glmmTMB", "gam", "gamm", "gls", "spamm", "brms", "inla", "bdlnm")
  )

  # Diagnostic timer only. This is measured in the calling R process and does
  # not alter candidate fitting or LOOCV execution.
  function_start_time <- proc.time()[["elapsed"]]

  # Force `...` once in the calling R session before any future is created.
  # This avoids exporting unresolved promises to multisession workers.
  fit_dots <- list(...)

  valid_scalar_name <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  valid_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!valid_scalar_name(response)) stop("`response` must be one non-empty column name.")
  if (!valid_scalar_name(group)) stop("`group` must be one non-empty column name.")
  if (!valid_scalar_name(time)) stop("`time` must be one non-empty column name.")
  required_basic <- c(response, group, time)
  missing_basic <- setdiff(required_basic, names(data))
  if (length(missing_basic)) {
    stop("Required columns missing from `data`: ", paste(missing_basic, collapse = ", "), ".")
  }
  if (!is.character(vars) || !length(vars) || anyNA(vars) ||
      any(!nzchar(vars)) || anyDuplicated(vars)) {
    stop("`vars` must contain unique non-empty variable names.")
  }
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop("Variables missing from `data`: ", paste(missing_vars, collapse = ", "), ".")
  }
  if (anyNA(data[[group]])) stop("The grouping column cannot contain missing values.")
  if (!is.numeric(data[[time]]) || any(!is.finite(data[[time]]))) {
    stop("The chronological time column must contain finite numeric values.")
  }
  if (!is.numeric(data[[response]]) || any(!is.finite(data[[response]]))) {
    stop("The response column must contain finite numeric values.")
  }
  for (variable in vars) {
    if (!is.numeric(data[[variable]]) || any(!is.finite(data[[variable]]))) {
      stop("Exposure variable '", variable, "' must contain finite numeric values.")
    }
  }
  if (!is.numeric(lag_max) || !length(lag_max) || any(!is.finite(lag_max)) || max(lag_max) < 0) {
    stop("`lag_max` must contain finite non-negative values.")
  }
  lag_max <- as.integer(max(lag_max))

  validate_df_grid <- function(x, argument) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x <= 0)) {
      stop("`", argument, "` must contain positive finite values.")
    }
    sort(unique(as.integer(x)))
  }
  df_var_grid <- validate_df_grid(df_var_grid, "df_var_grid")
  df_lag_grid <- validate_df_grid(df_lag_grid, "df_lag_grid")

  if (!is.numeric(min_vars) || length(min_vars) != 1L ||
      !is.finite(min_vars) || min_vars < 1) {
    stop("`min_vars` must be a positive integer.")
  }
  min_vars <- as.integer(min_vars)
  if (is.null(max_vars)) max_vars <- length(vars)
  if (!is.numeric(max_vars) || length(max_vars) != 1L || !is.finite(max_vars) ||
      max_vars < min_vars || max_vars > length(vars)) {
    stop("`max_vars` must be between `min_vars` and `length(vars)`.")
  }
  max_vars <- as.integer(max_vars)

  if (!is.null(random_effect)) {
    if (!valid_scalar_name(random_effect)) {
      stop("`random_effect` must be NULL or one non-empty column name.")
    }
    if (!random_effect %in% names(data)) {
      stop("`random_effect` ('", random_effect, "') was not found in `data`.")
    }
    if (anyNA(data[[random_effect]])) stop("`random_effect` cannot contain missing values.")
  }
  if (!is.numeric(min_success) || length(min_success) != 1L ||
      !is.finite(min_success) || min_success < 2) {
    stop("`min_success` must be an integer greater than or equal to 2.")
  }
  min_success <- as.integer(min_success)
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n <= 0) {
    stop("`top_n` must be a positive integer or `Inf`.")
  }
  if (is.finite(top_n)) top_n <- as.integer(top_n)
  if (!valid_flag(keep_fits)) stop("`keep_fits` must be TRUE or FALSE.")
  if (!valid_flag(verbose)) stop("`verbose` must be TRUE or FALSE.")

  data_long <- data
  data_long$epi_id <- data_long[[group]]
  data_long$time <- data_long[[time]]
  data_long$y <- data_long[[response]]
  data_long <- data_long[order(data_long$epi_id, data_long$time), , drop = FALSE]
  rownames(data_long) <- NULL

  response_check <- unique(data_long[, c("epi_id", "y"), drop = FALSE])
  if (anyDuplicated(response_check$epi_id)) {
    stop("Multiple distinct response values were found within at least one group.")
  }

  check_temporal_coverage <- function(input_data, variables, maximum_lag) {
    temporal_summary <- input_data |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        n_time = dplyr::n_distinct(time),
        has_missing_time = any(is.na(time)),
        .groups = "drop"
      )
    if (any(temporal_summary$has_missing_time)) {
      stop("Missing chronological time values were detected within groups.")
    }
    duplicate_groups <- temporal_summary[temporal_summary$n_rows != temporal_summary$n_time, , drop = FALSE]
    if (nrow(duplicate_groups)) {
      stop("Duplicated time values were detected within groups. Example group(s): ",
           paste(utils::head(duplicate_groups$epi_id, 5), collapse = ", "), ".")
    }
    insufficient <- temporal_summary[temporal_summary$n_time < maximum_lag + 1L, , drop = FALSE]
    if (nrow(insufficient)) {
      stop("Some groups have insufficient temporal coverage. Required observations: ",
           maximum_lag + 1L, ". Example group(s): ",
           paste(utils::head(insufficient$epi_id, 5), collapse = ", "), ".")
    }
    for (variable in variables) {
      finite_summary <- input_data |>
        dplyr::group_by(epi_id) |>
        dplyr::summarise(
          n_finite = sum(is.finite(.data[[variable]])),
          .groups = "drop"
        )
      bad <- finite_summary[finite_summary$n_finite < maximum_lag + 1L, , drop = FALSE]
      if (nrow(bad)) {
        stop("Insufficient finite observations for variable '", variable, "'.")
      }
    }
    invisible(TRUE)
  }

  get_linkinv_from_family <- function(family_object) {
    if (inherits(family_object, "family") && !is.null(family_object$linkinv)) {
      return(family_object$linkinv)
    }

    canonical_family <- .resolve_family_name(family_object)

    if (canonical_family %in% c("beta", "binomial")) return(stats::plogis)
    if (canonical_family %in% c("poisson", "gamma", "negative_binomial")) return(exp)
    if (canonical_family == "gaussian") return(identity)

    stop(
      "Could not determine the inverse-link function from `family` for prediction fallback.",
      call. = FALSE
    )
  }

  is_gamm_object <- function(model) {
    is.list(model) && !is.null(model$gam) && inherits(model$gam, "gam")
  }

  extract_fixed_coef <- function(model) {
    if (inherits(model, "glmmTMB")) return(glmmTMB::fixef(model)$cond)
    if (inherits(model, "merMod")) return(lme4::fixef(model))
    if (is_gamm_object(model)) return(stats::coef(model$gam))
    if (inherits(model, "lme")) return(nlme::fixef(model))
    if (inherits(model, "gls") || inherits(model, "gam") || inherits(model, "glm")) {
      return(stats::coef(model))
    }
    if (inherits(model, "HLfit")) return(spaMM::fixef(model))
    if (inherits(model, "brmsfit")) {
      fixed <- brms::fixef(model)
      beta <- fixed[, "Estimate"]
      names(beta) <- rownames(fixed)
      return(beta)
    }
    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      if (is.null(names(beta))) names(beta) <- rownames(model$summary.fixed)
      return(beta)
    }
    if (inherits(model, "bdlnm") && !is.null(model$coefficients.summary)) {
      return(model$coefficients.summary[, "mean"])
    }
    beta <- tryCatch(stats::coef(model), error = function(e) NULL)
    if (is.null(beta) || !is.numeric(beta)) {
      stop("Could not extract fixed-effect coefficients for prediction fallback.")
    }
    beta
  }

  predict_response_engine <- function(fitted_model, newdata, family_object) {
    direct <- tryCatch({
      if (inherits(fitted_model, "glmmTMB")) {
        stats::predict(fitted_model, newdata = newdata, type = "response", allow.new.levels = TRUE)
      } else if (inherits(fitted_model, "merMod")) {
        stats::predict(fitted_model, newdata = newdata, type = "response",
                       re.form = NA, allow.new.levels = TRUE)
      } else if (inherits(fitted_model, "brmsfit")) {
        if (!requireNamespace("brms", quietly = TRUE)) stop("Package 'brms' is required.")
        posterior_prediction <- brms::posterior_epred(
          fitted_model,
          newdata = newdata,
          re_formula = NA,
          allow_new_levels = TRUE
        )

        apply(
          posterior_prediction,
          2,
          stats::median,
          na.rm = TRUE
        )
      } else if (is_gamm_object(fitted_model)) {
        stats::predict(fitted_model$gam, newdata = newdata, type = "response")
      } else if (inherits(fitted_model, "gam") || inherits(fitted_model, "glm")) {
        stats::predict(fitted_model, newdata = newdata, type = "response")
      } else if (inherits(fitted_model, "gls")) {
        stats::predict(fitted_model, newdata = newdata)
      } else if (inherits(fitted_model, "lme")) {
        stats::predict(fitted_model, newdata = newdata, level = 0)
      } else if (inherits(fitted_model, "HLfit")) {
        stats::predict(fitted_model, newdata = newdata, type = "response")
      } else {
        stop("No direct prediction method available.")
      }
    }, error = function(e) NULL)

    if (!is.null(direct) && length(direct) == nrow(newdata) &&
        all(is.finite(as.numeric(direct)))) {
      return(as.numeric(direct))
    }

    beta <- extract_fixed_coef(fitted_model)
    coefficient_names <- names(beta)
    if (is.null(coefficient_names)) {
      stop("Prediction fallback failed because coefficients have no names.")
    }
    official_columns <- attr(fitted_model, "epiexposure_cb_cols")
    if (is.null(official_columns)) {
      template <- attr(fitted_model, "epiexposure_data_template")
      if (!is.null(template)) official_columns <- grep("^cb_", names(template), value = TRUE)
    }
    if (is.null(official_columns) || !length(official_columns)) {
      official_columns <- grep("^cb_", coefficient_names, value = TRUE)
    }
    missing_design <- setdiff(official_columns, names(newdata))
    missing_coefficients <- setdiff(official_columns, coefficient_names)
    if (length(missing_design)) {
      stop("Prediction fallback failed. Cross-basis columns missing from newdata: ",
           paste(missing_design, collapse = ", "), ".")
    }
    if (length(missing_coefficients)) {
      stop("Prediction fallback failed. Cross-basis coefficients missing from model: ",
           paste(missing_coefficients, collapse = ", "), ".")
    }
    eta <- rep(0, nrow(newdata))
    if ("(Intercept)" %in% coefficient_names) eta <- eta + as.numeric(beta["(Intercept)"])
    if (length(official_columns)) {
      design_matrix <- as.matrix(newdata[, official_columns, drop = FALSE])
      eta <- eta + as.numeric(design_matrix %*% as.numeric(beta[official_columns]))
    }
    linkinv <- get_linkinv_from_family(family_object)
    as.numeric(linkinv(eta))
  }

  # Resolve the same family used for model fitting into the canonical package
  # name used by the shared performance-metric infrastructure.
  family_name <- .resolve_family_name(family)
  outcome_type <- .resolve_outcome_type(family_name)

  if (identical(outcome_type, "binary") &&
      !all(response_check$y %in% c(0, 1))) {
    stop(
      "`family = 'binomial'` requires the group-level response used by ",
      "`find_bestfit()` to be coded as 0/1 so probability and classification ",
      "performance metrics can be computed consistently.",
      call. = FALSE
    )
  }

  if (identical(outcome_type, "ordinal")) {
    stop(
      "`find_bestfit()` does not yet define performance metrics for ",
      "`family = 'ordinal'`. Ordinal outcomes are recognized by EpiExposure, ",
      "but ordinal-specific ranking metrics must be defined before LOOCV model ",
      "selection can be performed.",
      call. = FALSE
    )
  }

  available_metrics <- .available_metrics(family_name)

  if (is.null(rank_metric)) {
    rank_metric <- if (identical(outcome_type, "binary")) "ROC_AUC" else "CCC"
  } else {
    if (!is.character(rank_metric) || length(rank_metric) != 1L ||
        is.na(rank_metric) || !nzchar(rank_metric)) {
      stop("`rank_metric` must be NULL or one non-empty metric name.", call. = FALSE)
    }
  }

  if (!rank_metric %in% available_metrics) {
    stop(
      "`rank_metric = ", sQuote(rank_metric), "` is not available for family '",
      family_name, "'. Available metrics are: ",
      paste(available_metrics, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0 || threshold >= 1) {
    stop(
      "`threshold` must be one finite numeric probability strictly between 0 and 1.",
      call. = FALSE
    )
  }

  rank_direction <- .metric_direction(rank_metric)

  if (verbose) {
    message(
      "Performance ranking: family = ", family_name,
      ", outcome_type = ", outcome_type,
      ", rank_metric = ", rank_metric,
      " (", rank_direction, ")."
    )

    if (identical(outcome_type, "binary")) {
      message(
        "Classification threshold = ", format(threshold),
        ". It is applied to predicted probabilities to generate 0/1 classes, ",
        "not to predictions already reduced to class labels. ",
        "ROC_AUC, Brier, and LogLoss are unaffected by the threshold."
      )
    }
  }

  build_candidate <- function(input_data, variables, exposure_df, lag_df) {
    templates <- define_exposure(
      data = input_data,
      vars = variables,
      lag_max = lag_max,
      df_var = exposure_df,
      df_lag = lag_df,
      fun_var = fun_var,
      fun_lag = fun_lag
    )
    design <- build_design(
      data = input_data,
      cb_templates = templates,
      lag_max = lag_max,
      include_response = TRUE
    )
    if (!is.null(random_effect)) {
      metadata <- unique(input_data[, c("epi_id", random_effect), drop = FALSE])
      if (anyDuplicated(metadata$epi_id)) {
        stop("`random_effect` must be unique within each epidemic.")
      }
      design <- merge(design, metadata, by = "epi_id", all.x = TRUE, sort = FALSE)
    }
    prepared <- prepare_response(
      data = design,
      y_var = "y",
      family = family_name
    )
    list(templates = templates, design = design, prepared = prepared)
  }

  fit_candidate <- function(candidate) {
    fit_call <- c(
      list(
        data = candidate$prepared,
        model_engine = model_engine,
        family = family,
        random_effect = random_effect,
        epiexposure_spec = attr(candidate$templates, "spec"),
        basis_objects = candidate$templates
      ),
      fit_dots
    )
    do.call(fit_epidlnm, fit_call)
  }

  if (is.null(var_sets)) {
    sizes <- seq.int(min_vars, max_vars)
    var_sets <- unlist(lapply(sizes, function(size) {
      utils::combn(vars, size, simplify = FALSE)
    }), recursive = FALSE)
  } else {
    if (!is.list(var_sets) || !length(var_sets)) {
      stop("`var_sets` must be NULL or a non-empty list.")
    }
    var_sets <- lapply(var_sets, function(set) {
      if (!is.character(set) || !length(set) || anyNA(set) ||
          any(!nzchar(set)) || anyDuplicated(set)) {
        stop("Every element of `var_sets` must contain unique non-empty variable names.")
      }
      missing <- setdiff(set, vars)
      if (length(missing)) {
        stop("Variables in `var_sets` not listed in `vars`: ",
             paste(missing, collapse = ", "), ".")
      }
      sort(set)
    })
    keys <- vapply(var_sets, paste, collapse = "||", character(1))
    var_sets <- var_sets[!duplicated(keys)]
  }

  check_temporal_coverage(data_long, vars, lag_max)
  fold_ids <- unique(data_long$epi_id)
  if (length(fold_ids) < 2L) stop("At least two groups are required for LOOCV.")
  if (min_success > length(fold_ids)) {
    stop("`min_success` cannot exceed the number of LOOCV groups.")
  }

  n_folds <- length(fold_ids)

  # Pre-compute only the row indices of each held-out epidemic.  The exposure
  # bases, lag structures, templates and model matrices are still rebuilt from
  # the training data inside every LOOCV fold.
  fold_code <- match(data_long$epi_id, fold_ids)
  fold_rows <- split(seq_len(nrow(data_long)), fold_code)
  fold_rows <- unname(fold_rows[as.character(seq_len(n_folds))])
  test_data_list <- lapply(
    fold_rows,
    function(idx) data_long[idx, , drop = FALSE]
  )

  # Classify warning messages for diagnostics only. Classification never
  # changes whether a fold is successful, candidate ranking, or model fitting.
  warning_type_levels <- c("knot", "convergence", "hessian", "other")

  classify_warning <- function(message) {
    message_lower <- tolower(as.character(message)[1])

    # Check Hessian warnings before the broader convergence class because some
    # glmmTMB messages include both "convergence problem" and "Hessian".
    if (grepl("hessian", message_lower, fixed = TRUE) ||
        grepl("positive-definite", message_lower, fixed = TRUE) ||
        grepl("positive definite", message_lower, fixed = TRUE)) {
      return("hessian")
    }

    if (grepl("shoving 'interior' knots matching boundary knots to inside",
              message_lower, fixed = TRUE) ||
        (grepl("knot", message_lower, fixed = TRUE) &&
         grepl("boundary", message_lower, fixed = TRUE))) {
      return("knot")
    }

    if (grepl("converg", message_lower) ||
        grepl("function evaluation limit reached", message_lower, fixed = TRUE) ||
        grepl("iteration limit", message_lower, fixed = TRUE) ||
        grepl("maximum number of iterations", message_lower, fixed = TRUE)) {
      return("convergence")
    }

    "other"
  }

  # Reproduce exactly the order of the original nested loops:
  # df_var (outer) -> df_lag -> variable set (inner).
  candidate_grid <- expand.grid(
    var_set_id = seq_along(var_sets),
    df_lag = df_lag_grid,
    df_var = df_var_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  candidate_grid <- candidate_grid[, c("df_var", "df_lag", "var_set_id"), drop = FALSE]
  candidate_grid$model_id <- seq_len(nrow(candidate_grid))
  total_candidates <- nrow(candidate_grid)

  evaluate_candidate <- function(candidate_index) {
    candidate_start_time <- proc.time()[["elapsed"]]

    candidate_row <- candidate_grid[candidate_index, , drop = FALSE]
    model_id <- candidate_row$model_id[[1]]
    exposure_df <- candidate_row$df_var[[1]]
    lag_df <- candidate_row$df_lag[[1]]
    variables <- var_sets[[candidate_row$var_set_id[[1]]]]
    variables_label <- paste(variables, collapse = " + ")

    if (verbose) {
      message(
        "[", candidate_index, "/", total_candidates, "] df_var = ",
        exposure_df, ", df_lag = ", lag_df,
        ", vars = ", variables_label
      )
    }

    observed <- rep(NA_real_, n_folds)
    predicted <- rep(NA_real_, n_folds)
    success <- rep(FALSE, n_folds)

    # Warning bookkeeping is local to this candidate/future. Warnings are
    # muffled inside each fold and summarized only after all candidates finish.
    warning_fold <- rep(FALSE, n_folds)
    warning_events <- integer(n_folds)
    candidate_warning_messages <- vector("list", n_folds)
    warning_type_fold <- matrix(
      FALSE,
      nrow = n_folds,
      ncol = length(warning_type_levels),
      dimnames = list(NULL, warning_type_levels)
    )
    warning_type_events <- matrix(
      0L,
      nrow = n_folds,
      ncol = length(warning_type_levels),
      dimnames = list(NULL, warning_type_levels)
    )

    # Slots are indexed by fold so their order is deterministic even when
    # candidates themselves are evaluated in parallel.
    candidate_predictions <- vector("list", n_folds)
    candidate_failures <- vector("list", n_folds + as.integer(keep_fits))
    # One additional slot is reserved for candidate-level performance-metric
    # warnings; a second additional slot is used for the full-data refit when
    # keep_fits = TRUE.
    candidate_warnings <- vector("list", n_folds + 1L + as.integer(keep_fits))

    # IMPORTANT: this loop is deliberately sequential. One future owns one
    # candidate and executes the complete LOOCV workflow for that candidate.
    for (fold_index in seq_len(n_folds)) {
      test_group <- fold_ids[fold_index]
      test_idx <- fold_rows[[fold_index]]

      # Negative integer subsetting avoids rescanning epi_id with != at every
      # candidate/fold combination. No training information is reused.
      train_data <- data_long[-test_idx, , drop = FALSE]
      test_data <- test_data_list[[fold_index]]

      fold_warning_messages <- character(0)

      fold_result <- tryCatch(
        withCallingHandlers(
          {
            # Templates are estimated ONLY on the training epidemics.
            training_candidate <- build_candidate(
              train_data,
              variables,
              exposure_df,
              lag_df
            )
            fitted_model <- fit_candidate(training_candidate)

            # The held-out epidemic receives the templates learned from training.
            test_design <- build_design(
              data = test_data,
              cb_templates = training_candidate$templates,
              lag_max = lag_max,
              include_response = TRUE
            )

            if (!is.null(random_effect)) {
              test_metadata <- unique(
                test_data[, c("epi_id", random_effect), drop = FALSE]
              )
              if (anyDuplicated(test_metadata$epi_id)) {
                stop("`random_effect` must be unique within the held-out epidemic.")
              }
              test_design <- merge(
                test_design,
                test_metadata,
                by = "epi_id",
                all.x = TRUE,
                sort = FALSE
              )
            }

            prepared_test <- prepare_response(
              data = test_design,
              y_var = "y",
              family = family_name
            )
            prediction <- predict_response_engine(
              fitted_model,
              prepared_test,
              family
            )

            list(
              ok = TRUE,
              obs = prepared_test$y_model[1],
              pred = as.numeric(prediction[1])
            )
          },
          warning = function(w) {
            fold_warning_messages <<- c(
              fold_warning_messages,
              conditionMessage(w)
            )
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          list(ok = FALSE, error = conditionMessage(e))
        }
      )

      if (length(fold_warning_messages)) {
        warning_fold[fold_index] <- TRUE
        warning_events[fold_index] <- length(fold_warning_messages)
        candidate_warning_messages[[fold_index]] <- fold_warning_messages

        fold_warning_types <- vapply(
          fold_warning_messages,
          classify_warning,
          character(1)
        )
        fold_type_counts <- table(factor(
          fold_warning_types,
          levels = warning_type_levels
        ))
        warning_type_fold[fold_index, ] <- as.integer(fold_type_counts) > 0L
        warning_type_events[fold_index, ] <- as.integer(fold_type_counts)

        candidate_warnings[[fold_index]] <- data.frame(
          model_id = model_id,
          fold = fold_index,
          group = as.character(test_group),
          stage = "LOOCV",
          df_var = exposure_df,
          df_lag = lag_df,
          vars = variables_label,
          n_warning_events = length(fold_warning_messages),
          warning_types = paste(unique(fold_warning_types), collapse = " | "),
          n_knot_warning_events = as.integer(fold_type_counts[["knot"]]),
          n_convergence_warning_events = as.integer(fold_type_counts[["convergence"]]),
          n_hessian_warning_events = as.integer(fold_type_counts[["hessian"]]),
          n_other_warning_events = as.integer(fold_type_counts[["other"]]),
          warning = paste(unique(fold_warning_messages), collapse = " | "),
          stringsAsFactors = FALSE
        )
      }

      if (!isTRUE(fold_result$ok)) {
        candidate_failures[[fold_index]] <- data.frame(
          model_id = model_id,
          fold = fold_index,
          group = as.character(test_group),
          df_var = exposure_df,
          df_lag = lag_df,
          vars = variables_label,
          error = fold_result$error,
          stringsAsFactors = FALSE
        )
        next
      }

      if (is.finite(fold_result$obs) && is.finite(fold_result$pred)) {
        observed[fold_index] <- fold_result$obs
        predicted[fold_index] <- fold_result$pred
        success[fold_index] <- TRUE

        prediction_row <- data.frame(
          model_id = model_id,
          fold = fold_index,
          group = as.character(test_group),
          df_var = exposure_df,
          df_lag = lag_df,
          vars = variables_label,
          observed = fold_result$obs,
          predicted = fold_result$pred,
          stringsAsFactors = FALSE
        )

        # For binomial outcomes, `predicted` remains the out-of-fold predicted
        # probability. The thresholded class is stored separately so no
        # probability information is discarded.
        if (identical(outcome_type, "binary")) {
          prediction_row$predicted_class <- as.integer(
            fold_result$pred >= threshold
          )
        }

        candidate_predictions[[fold_index]] <- prediction_row
      }
    }

    n_success <- sum(success)
    n_failed <- n_folds - n_success
    n_warning_folds <- sum(warning_fold)
    n_warning_events <- sum(warning_events)
    warning_rate <- n_warning_folds / n_folds

    warning_folds_by_type <- colSums(warning_type_fold)
    warning_events_by_type <- colSums(warning_type_events)

    all_loocv_warning_messages <- unlist(
      candidate_warning_messages,
      recursive = FALSE,
      use.names = FALSE
    )
    all_loocv_warning_messages <- as.character(all_loocv_warning_messages)
    unique_loocv_warning_messages <- unique(all_loocv_warning_messages)

    candidate_result <- NULL
    metric_warning_messages <- character(0)
    full_fit <- NULL
    full_fit_warning_messages <- character(0)
    full_fit_warning_events_by_type <- setNames(
      integer(length(warning_type_levels)),
      warning_type_levels
    )

    if (n_success >= min_success) {
      metrics <- withCallingHandlers(
        .compute_performance_metrics(
          observed = observed[success],
          predicted = predicted[success],
          family = family_name,
          threshold = threshold
        ),
        warning = function(w) {
          metric_warning_messages <<- c(
            metric_warning_messages,
            conditionMessage(w)
          )
          invokeRestart("muffleWarning")
        }
      )

      metric_warning_messages <- unique(metric_warning_messages)

      if (length(metric_warning_messages)) {
        candidate_warnings[[n_folds + 1L]] <- data.frame(
          model_id = model_id,
          fold = NA_integer_,
          group = NA_character_,
          stage = "metrics",
          df_var = exposure_df,
          df_lag = lag_df,
          vars = variables_label,
          n_warning_events = length(metric_warning_messages),
          warning_types = "other",
          n_knot_warning_events = 0L,
          n_convergence_warning_events = 0L,
          n_hessian_warning_events = 0L,
          n_other_warning_events = length(metric_warning_messages),
          warning = paste(metric_warning_messages, collapse = " | "),
          stringsAsFactors = FALSE
        )
      }

      candidate_identity <- data.frame(
        model_id = model_id,
        df_var = exposure_df,
        df_lag = lag_df,
        vars = variables_label,
        n_vars = length(variables),
        stringsAsFactors = FALSE
      )

      candidate_diagnostics <- data.frame(
        n_folds = n_folds,
        n_success = n_success,
        n_failed = n_failed,
        n_warning_folds = n_warning_folds,
        warning_rate = warning_rate,
        n_warning_events = n_warning_events,
        n_knot_warning_folds = as.integer(warning_folds_by_type[["knot"]]),
        n_convergence_warning_folds = as.integer(warning_folds_by_type[["convergence"]]),
        n_hessian_warning_folds = as.integer(warning_folds_by_type[["hessian"]]),
        n_other_warning_folds = as.integer(warning_folds_by_type[["other"]]),
        n_knot_warning_events = as.integer(warning_events_by_type[["knot"]]),
        n_convergence_warning_events = as.integer(warning_events_by_type[["convergence"]]),
        n_hessian_warning_events = as.integer(warning_events_by_type[["hessian"]]),
        n_other_warning_events = as.integer(warning_events_by_type[["other"]]),
        stringsAsFactors = FALSE
      )

      candidate_result <- cbind(
        candidate_identity,
        metrics,
        candidate_diagnostics
      )

      # Preserve the original behavior: every candidate that reaches
      # min_success is refitted on all data when keep_fits = TRUE. Ranking and
      # top_n filtering are still applied later, exactly as before.
      if (keep_fits) {
        full_fit_result <- tryCatch(
          withCallingHandlers(
            {
              list(
                ok = TRUE,
                fit = fit_candidate(
                  build_candidate(
                    data_long,
                    variables,
                    exposure_df,
                    lag_df
                  )
                )
              )
            },
            warning = function(w) {
              full_fit_warning_messages <<- c(
                full_fit_warning_messages,
                conditionMessage(w)
              )
              invokeRestart("muffleWarning")
            }
          ),
          error = function(e) {
            list(ok = FALSE, error = conditionMessage(e))
          }
        )

        if (length(full_fit_warning_messages)) {
          full_fit_warning_types <- vapply(
            full_fit_warning_messages,
            classify_warning,
            character(1)
          )
          full_fit_type_counts <- table(factor(
            full_fit_warning_types,
            levels = warning_type_levels
          ))
          full_fit_warning_events_by_type[] <- as.integer(full_fit_type_counts)

          candidate_warnings[[n_folds + 2L]] <- data.frame(
            model_id = model_id,
            fold = NA_integer_,
            group = NA_character_,
            stage = "full_refit",
            df_var = exposure_df,
            df_lag = lag_df,
            vars = variables_label,
            n_warning_events = length(full_fit_warning_messages),
            warning_types = paste(unique(full_fit_warning_types), collapse = " | "),
            n_knot_warning_events = as.integer(full_fit_type_counts[["knot"]]),
            n_convergence_warning_events = as.integer(full_fit_type_counts[["convergence"]]),
            n_hessian_warning_events = as.integer(full_fit_type_counts[["hessian"]]),
            n_other_warning_events = as.integer(full_fit_type_counts[["other"]]),
            warning = paste(unique(full_fit_warning_messages), collapse = " | "),
            stringsAsFactors = FALSE
          )
        }

        if (isTRUE(full_fit_result$ok)) {
          full_fit <- full_fit_result$fit
        } else {
          candidate_failures[[n_folds + 1L]] <- data.frame(
            model_id = model_id,
            fold = NA_integer_,
            group = NA_character_,
            df_var = exposure_df,
            df_lag = lag_df,
            vars = variables_label,
            error = paste0(
              "Full-data refit failed: ",
              full_fit_result$error
            ),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    warning_summary <- data.frame(
      model_id = model_id,
      df_var = exposure_df,
      df_lag = lag_df,
      vars = variables_label,
      n_folds = n_folds,
      n_warning_folds = n_warning_folds,
      warning_rate = warning_rate,
      n_warning_events = n_warning_events,
      n_knot_warning_folds = as.integer(warning_folds_by_type[["knot"]]),
      n_convergence_warning_folds = as.integer(warning_folds_by_type[["convergence"]]),
      n_hessian_warning_folds = as.integer(warning_folds_by_type[["hessian"]]),
      n_other_warning_folds = as.integer(warning_folds_by_type[["other"]]),
      n_knot_warning_events = as.integer(warning_events_by_type[["knot"]]),
      n_convergence_warning_events = as.integer(warning_events_by_type[["convergence"]]),
      n_hessian_warning_events = as.integer(warning_events_by_type[["hessian"]]),
      n_other_warning_events = as.integer(warning_events_by_type[["other"]]),
      n_warning_types = length(unique_loocv_warning_messages),
      n_warning_classes = length(unique(vapply(
        unique_loocv_warning_messages,
        classify_warning,
        character(1)
      ))),
      warning_classes = if (length(unique_loocv_warning_messages)) {
        paste(unique(vapply(
          unique_loocv_warning_messages,
          classify_warning,
          character(1)
        )), collapse = " | ")
      } else {
        NA_character_
      },
      warning_messages = if (length(unique_loocv_warning_messages)) {
        paste(unique_loocv_warning_messages, collapse = " | ")
      } else {
        NA_character_
      },
      metric_warning_events = length(metric_warning_messages),
      metric_warning_messages = if (length(metric_warning_messages)) {
        paste(metric_warning_messages, collapse = " | ")
      } else {
        NA_character_
      },
      full_refit_warning_events = length(full_fit_warning_messages),
      full_refit_knot_warning_events = as.integer(full_fit_warning_events_by_type[["knot"]]),
      full_refit_convergence_warning_events = as.integer(full_fit_warning_events_by_type[["convergence"]]),
      full_refit_hessian_warning_events = as.integer(full_fit_warning_events_by_type[["hessian"]]),
      full_refit_other_warning_events = as.integer(full_fit_warning_events_by_type[["other"]]),
      full_refit_warning_classes = if (length(full_fit_warning_messages)) {
        paste(unique(vapply(
          full_fit_warning_messages,
          classify_warning,
          character(1)
        )), collapse = " | ")
      } else {
        NA_character_
      },
      full_refit_warning_messages = if (length(full_fit_warning_messages)) {
        paste(unique(full_fit_warning_messages), collapse = " | ")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )

    candidate_elapsed <- proc.time()[["elapsed"]] - candidate_start_time

    list(
      model_id = model_id,
      result = candidate_result,
      predictions = Filter(Negate(is.null), candidate_predictions),
      failures = Filter(Negate(is.null), candidate_failures),
      warnings = Filter(Negate(is.null), candidate_warnings),
      warning_summary = warning_summary,
      fit = full_fit,
      elapsed_seconds = as.numeric(candidate_elapsed)
    )
  }

  # Package functions should not change the user's global future::plan().
  # Therefore parallelism is used whenever the caller has configured >1 worker;
  # otherwise the exact same evaluator runs sequentially.
  use_future <-
    total_candidates > 1L &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE) &&
    future::nbrOfWorkers() > 1L

  if (use_future) {
    if (verbose) {
      message(
        "Evaluating ", total_candidates,
        " candidate models across ", future::nbrOfWorkers(),
        " future workers; LOOCV folds remain sequential within each candidate."
      )
    }

    candidate_outputs <- future.apply::future_lapply(
      X = seq_len(total_candidates),
      FUN = evaluate_candidate,
      future.seed = TRUE,
      # One future per candidate: folds remain sequential inside each future.
      future.scheduling = Inf
    )
  } else {
    if (verbose && total_candidates > 1L) {
      message(
        "Evaluating candidates sequentially. To enable candidate-level ",
        "parallelism, set a future plan with >1 worker before calling ",
        "find_bestfit()."
      )
    }
    candidate_outputs <- lapply(
      seq_len(total_candidates),
      evaluate_candidate
    )
  }

  results <- Filter(
    Negate(is.null),
    lapply(candidate_outputs, function(x) x$result)
  )
  predictions <- unlist(
    lapply(candidate_outputs, function(x) x$predictions),
    recursive = FALSE,
    use.names = FALSE
  )
  failures <- unlist(
    lapply(candidate_outputs, function(x) x$failures),
    recursive = FALSE,
    use.names = FALSE
  )
  warnings <- unlist(
    lapply(candidate_outputs, function(x) x$warnings),
    recursive = FALSE,
    use.names = FALSE
  )

  warning_summary_data <- do.call(
    rbind,
    lapply(candidate_outputs, function(x) x$warning_summary)
  )
  rownames(warning_summary_data) <- NULL

  candidate_times <- data.frame(
    model_id = vapply(
      candidate_outputs,
      function(x) as.integer(x$model_id),
      integer(1)
    ),
    elapsed_seconds = vapply(
      candidate_outputs,
      function(x) as.numeric(x$elapsed_seconds),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  # Warnings have already been muffled inside each fold. Emit at most one
  # concise summary per affected candidate, in deterministic model_id order.
  if (verbose && nrow(warning_summary_data)) {
    affected <- warning_summary_data$n_warning_folds > 0L |
      warning_summary_data$metric_warning_events > 0L |
      warning_summary_data$full_refit_warning_events > 0L

    for (i in which(affected)) {
      ws <- warning_summary_data[i, , drop = FALSE]
      loocv_text <- if (ws$n_warning_folds > 0L) {
        paste0(
          ws$n_warning_folds, "/", ws$n_folds,
          " LOOCV folds (", ws$n_warning_events, " warning event(s))"
        )
      } else {
        "0 LOOCV folds"
      }

      type_fold_counts <- c(
        knot = ws$n_knot_warning_folds,
        convergence = ws$n_convergence_warning_folds,
        hessian = ws$n_hessian_warning_folds,
        other = ws$n_other_warning_folds
      )
      type_fold_counts <- type_fold_counts[type_fold_counts > 0L]
      type_text <- if (length(type_fold_counts)) {
        paste0(
          "; warning classes (affected folds): ",
          paste(
            paste0(names(type_fold_counts), "=", as.integer(type_fold_counts)),
            collapse = ", "
          )
        )
      } else {
        ""
      }

      metric_text <- if (ws$metric_warning_events > 0L) {
        paste0(
          "; performance metrics: ", ws$metric_warning_events,
          " warning event(s)"
        )
      } else {
        ""
      }

      full_fit_text <- if (ws$full_refit_warning_events > 0L) {
        paste0(
          "; full-data refit: ", ws$full_refit_warning_events,
          " warning event(s)"
        )
      } else {
        ""
      }

      message(
        "Warning summary [model ", ws$model_id, "/", total_candidates,
        "] df_var = ", ws$df_var,
        ", df_lag = ", ws$df_lag,
        ", vars = ", ws$vars,
        ": ", loocv_text, type_text, metric_text, full_fit_text,
        ". See attr(result, \"warning_summary\") and attr(result, \"warnings\") for details."
      )
    }
  }

  fits_list <- list()
  if (keep_fits) {
    for (candidate_output in candidate_outputs) {
      if (!is.null(candidate_output$fit)) {
        fits_list[[as.character(candidate_output$model_id)]] <- candidate_output$fit
      }
    }
  }

  if (!length(results)) stop("No candidate model produced enough successful LOOCV predictions.")
  results_data <- do.call(rbind, results)

  # Rank primarily by the user-selected metric. Remaining family-appropriate
  # metrics are deterministic tie-breakers, each using its own known direction.
  ranking_metrics <- c(
    rank_metric,
    setdiff(available_metrics, rank_metric)
  )

  if (all(is.na(results_data[[rank_metric]]))) {
    warning(
      "`rank_metric = ", rank_metric,
      "` is NA for all eligible candidates. Ranking will therefore be ",
      "resolved by the remaining available metrics in canonical order.",
      call. = FALSE
    )
  }

  ranking_vectors <- lapply(
    ranking_metrics,
    function(metric_name) {
      values <- results_data[[metric_name]]
      if (identical(.metric_direction(metric_name), "maximize")) {
        -values
      } else {
        values
      }
    }
  )

  ranking_order <- do.call(
    order,
    c(
      ranking_vectors,
      list(results_data$model_id),
      list(na.last = TRUE)
    )
  )

  results_data <- results_data[ranking_order, , drop = FALSE]
  results_data$rank <- seq_len(nrow(results_data))
  results_data <- results_data[, c("rank", setdiff(names(results_data), "rank")), drop = FALSE]
  if (is.finite(top_n)) results_data <- utils::head(results_data, top_n)
  rownames(results_data) <- NULL

  prediction_data <- if (length(predictions)) do.call(rbind, predictions) else data.frame()
  failure_data <- if (length(failures)) do.call(rbind, failures) else data.frame()
  warning_data <- if (length(warnings)) do.call(rbind, warnings) else data.frame()

  # Preserve the original flat prediction data frame and additionally provide
  # a convenient list split by model_id for direct model-specific extraction.
  predictions_by_model <- if (nrow(prediction_data)) {
    split_predictions <- split(
      prediction_data,
      as.character(prediction_data$model_id),
      drop = TRUE
    )
    lapply(split_predictions, function(x) {
      rownames(x) <- NULL
      x
    })
  } else {
    list()
  }

  total_elapsed_seconds <- as.numeric(
    proc.time()[["elapsed"]] - function_start_time
  )

  format_elapsed <- function(seconds) {
    seconds <- max(0, as.numeric(seconds))
    hours <- floor(seconds / 3600)
    minutes <- floor((seconds %% 3600) / 60)
    secs <- seconds %% 60
    sprintf("%02d:%02d:%05.2f", as.integer(hours), as.integer(minutes), secs)
  }

  timing_data <- data.frame(
    total_elapsed_seconds = total_elapsed_seconds,
    total_elapsed = format_elapsed(total_elapsed_seconds),
    total_candidates = total_candidates,
    n_folds = n_folds,
    parallel = use_future,
    workers = if (use_future) future::nbrOfWorkers() else 1L,
    stringsAsFactors = FALSE
  )

  attr(results_data, "predictions") <- prediction_data
  attr(results_data, "predictions_by_model") <- predictions_by_model
  attr(results_data, "failures") <- failure_data
  attr(results_data, "warnings") <- warning_data
  attr(results_data, "warning_summary") <- warning_summary_data
  attr(results_data, "timing") <- timing_data
  attr(results_data, "candidate_times") <- candidate_times
  attr(results_data, "family") <- family_name
  attr(results_data, "outcome_type") <- outcome_type
  attr(results_data, "rank_metric") <- rank_metric
  attr(results_data, "threshold") <- threshold

  if (keep_fits) {
    retained_ids <- as.character(results_data$model_id)
    attr(results_data, "fits") <- fits_list[names(fits_list) %in% retained_ids]
  }

  if (verbose) {
    message(
      "find_bestfit completed in ",
      timing_data$total_elapsed,
      " (", format(round(total_elapsed_seconds, 2), nsmall = 2),
      " s). Ranked by ", rank_metric, " (", rank_direction, ")."
    )
  }

  results_data
}
