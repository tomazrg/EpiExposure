#' Find the best DLNM model structure using LOOCV and Lin's CCC
#'
#' Tests combinations of exposure variables and basis dimensions using
#' leave-one-group-out cross-validation. Candidate models are ranked by Lin's
#' concordance correlation coefficient (CCC), with Cb, Pearson correlation,
#' RMSE, and MAE also reported.
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
#' @param max_lag Maximum lag. A numeric vector such as `c(0, 85)` is normalized
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
#'   `fit_epidlnm()`.
#' @param random_effect Optional character scalar naming a random-effect column.
#' @param min_success Minimum number of successful LOOCV folds required.
#' @param top_n Positive integer or `Inf`; number of ranked candidates returned.
#' @param keep_fits Logical. Refit successful retained candidates on all data and
#'   store them in `attr(result, "fits")`.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments passed to `fit_epidlnm()`.
#'
#' @return A data frame ranked by decreasing CCC, with attributes
#'   `"predictions"`, `"failures"`, and optionally `"fits"`.
#'
#' @details
#' For each fold, basis templates are estimated from the training groups only
#' and then applied to the held-out group. This prevents the held-out exposure
#' history from determining the training basis specification.
#'
#' Direct engine-specific prediction is attempted first. If unavailable, a
#' fixed-effect fallback uses the official cross-basis column order stored in
#' `epiexposure_cb_cols`, or `epiexposure_data_template` as a fallback. Missing
#' model or design columns produce an explicit error rather than being omitted.
#'
#' The function ranks point predictions and does not propagate coefficient or
#' posterior-draw uncertainty.
#'
#' @export
find_bestfit <- function(
    data,
    response = "y",
    group = "epi_id",
    time = "time",
    vars,
    max_lag,
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
    top_n = Inf,
    keep_fits = FALSE,
    verbose = TRUE,
    ...
) {
  model_engine <- match.arg(
    model_engine,
    choices = c("glm", "glmmTMB", "gam", "gamm", "gls", "spamm", "brms", "inla", "bdlnm")
  )

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
  if (!is.numeric(max_lag) || !length(max_lag) || any(!is.finite(max_lag)) || max(max_lag) < 0) {
    stop("`max_lag` must contain finite non-negative values.")
  }
  max_lag <- as.integer(max(max_lag))

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

  resolve_family_name <- function(family) {

    if (is.character(family)) {
      return(tolower(family[[1]]))
    }

    if (inherits(family, "family")) {
      return(tolower(family$family))
    }

    stop(
      "Unsupported family specification."
    )

  }

  ccc_lins <- function(obs, pred) {
    valid <- is.finite(obs) & is.finite(pred)
    obs <- obs[valid]
    pred <- pred[valid]
    if (length(obs) < 2L) {
      return(list(CCC = NA_real_, Cb = NA_real_, rho = NA_real_,
                  RMSE = NA_real_, MAE = NA_real_))
    }
    mean_obs <- mean(obs)
    mean_pred <- mean(pred)
    var_obs <- stats::var(obs)
    var_pred <- stats::var(pred)
    covariance <- stats::cov(obs, pred)
    rho <- suppressWarnings(stats::cor(obs, pred))
    ccc <- (2 * covariance) / (var_obs + var_pred + (mean_obs - mean_pred)^2)
    cb <- if (is.finite(rho) && abs(rho) > .Machine$double.eps) ccc / rho else NA_real_
    list(
      CCC = as.numeric(ccc), Cb = as.numeric(cb), rho = as.numeric(rho),
      RMSE = sqrt(mean((obs - pred)^2)), MAE = mean(abs(obs - pred))
    )
  }

  get_linkinv_from_family <- function(family_object) {
    if (inherits(family_object, "family") && !is.null(family_object$linkinv)) {
      return(family_object$linkinv)
    }
    if (is.character(family_object)) {
      family_name <- tolower(family_object[1])
      if (family_name %in% c("beta", "binomial")) return(stats::plogis)
      if (family_name %in% c("poisson", "gamma", "negbin", "negative_binomial")) return(exp)
      if (family_name == "gaussian") return(identity)
    }
    stop("Could not determine the inverse-link function from `family`.")
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

  build_candidate <- function(input_data, variables, exposure_df, lag_df) {
    templates <- define_exposure(
      data = input_data,
      vars = variables,
      max_lag = max_lag,
      df_var = exposure_df,
      df_lag = lag_df,
      fun_var = fun_var,
      fun_lag = fun_lag
    )
    design <- build_design(
      data = input_data,
      cb_templates = templates,
      max_lag = max_lag,
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
      family = resolve_family_name(family)
    )
    list(templates = templates, design = design, prepared = prepared)
  }

  fit_candidate <- function(candidate) {
    fit_epidlnm(
      data = candidate$prepared,
      model_engine = model_engine,
      family = family,
      random_effect = random_effect,
      epiexposure_spec = attr(candidate$templates, "spec"),
      basis_objects = candidate$templates,
      ...
    )
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

  check_temporal_coverage(data_long, vars, max_lag)
  fold_ids <- unique(data_long$epi_id)
  if (length(fold_ids) < 2L) stop("At least two groups are required for LOOCV.")
  if (min_success > length(fold_ids)) {
    stop("`min_success` cannot exceed the number of LOOCV groups.")
  }

  results <- list()
  predictions <- list()
  failures <- list()
  fits_list <- list()
  model_id <- 1L
  prediction_id <- 1L
  failure_id <- 1L
  total_candidates <- length(df_var_grid) * length(df_lag_grid) * length(var_sets)
  candidate_id <- 1L

  for (exposure_df in df_var_grid) {
    for (lag_df in df_lag_grid) {
      for (variables in var_sets) {
        if (verbose) {
          message("[", candidate_id, "/", total_candidates, "] df_var = ",
                  exposure_df, ", df_lag = ", lag_df,
                  ", vars = ", paste(variables, collapse = " + "))
        }
        candidate_id <- candidate_id + 1L
        observed <- rep(NA_real_, length(fold_ids))
        predicted <- rep(NA_real_, length(fold_ids))
        success <- rep(FALSE, length(fold_ids))

        for (fold_index in seq_along(fold_ids)) {
          test_group <- fold_ids[fold_index]
          train_data <- data_long[data_long$epi_id != test_group, , drop = FALSE]
          test_data <- data_long[data_long$epi_id == test_group, , drop = FALSE]

          fold_result <- tryCatch({
            training_candidate <- build_candidate(train_data, variables, exposure_df, lag_df)
            fitted_model <- fit_candidate(training_candidate)
            test_templates <- training_candidate$templates
            test_design <- build_design(
              data = test_data,
              cb_templates = test_templates,
              max_lag = max_lag,
              include_response = TRUE
            )
            if (!is.null(random_effect)) {
              test_metadata <- unique(test_data[, c("epi_id", random_effect), drop = FALSE])
              if (anyDuplicated(test_metadata$epi_id)) {
                stop("`random_effect` must be unique within the held-out epidemic.")
              }
              test_design <- merge(test_design, test_metadata, by = "epi_id",
                                   all.x = TRUE, sort = FALSE)
            }
            prepared_test <- prepare_response(
              data = test_design,
              y_var = "y",
              family = resolve_family_name(family)
            )
            prediction <- predict_response_engine(fitted_model, prepared_test, family)
            list(obs = prepared_test$y_model[1], pred = as.numeric(prediction[1]))
          }, error = function(e) {
            failures[[failure_id]] <<- data.frame(
              model_id = model_id, fold = fold_index, group = as.character(test_group),
              df_var = exposure_df, df_lag = lag_df,
              vars = paste(variables, collapse = " + "),
              error = conditionMessage(e), stringsAsFactors = FALSE
            )
            failure_id <<- failure_id + 1L
            NULL
          })

          if (!is.null(fold_result) && is.finite(fold_result$obs) && is.finite(fold_result$pred)) {
            observed[fold_index] <- fold_result$obs
            predicted[fold_index] <- fold_result$pred
            success[fold_index] <- TRUE
            predictions[[prediction_id]] <- data.frame(
              model_id = model_id, fold = fold_index, group = as.character(test_group),
              df_var = exposure_df, df_lag = lag_df,
              vars = paste(variables, collapse = " + "),
              observed = fold_result$obs, predicted = fold_result$pred,
              stringsAsFactors = FALSE
            )
            prediction_id <- prediction_id + 1L
          }
        }

        n_success <- sum(success)
        n_failed <- length(fold_ids) - n_success
        if (n_success >= min_success) {
          metrics <- ccc_lins(observed[success], predicted[success])
          results[[length(results) + 1L]] <- data.frame(
            model_id = model_id, df_var = exposure_df, df_lag = lag_df,
            vars = paste(variables, collapse = " + "), n_vars = length(variables),
            CCC = metrics$CCC, Cb = metrics$Cb, rho = metrics$rho,
            RMSE = metrics$RMSE, MAE = metrics$MAE,
            n_folds = length(fold_ids), n_success = n_success, n_failed = n_failed,
            stringsAsFactors = FALSE
          )
          if (keep_fits) {
            full_fit <- tryCatch({
              fit_candidate(build_candidate(data_long, variables, exposure_df, lag_df))
            }, error = function(e) {
              failures[[failure_id]] <<- data.frame(
                model_id = model_id, fold = NA_integer_, group = NA_character_,
                df_var = exposure_df, df_lag = lag_df,
                vars = paste(variables, collapse = " + "),
                error = paste0("Full-data refit failed: ", conditionMessage(e)),
                stringsAsFactors = FALSE
              )
              failure_id <<- failure_id + 1L
              NULL
            })
            if (!is.null(full_fit)) fits_list[[as.character(model_id)]] <- full_fit
          }
        }
        model_id <- model_id + 1L
      }
    }
  }

  if (!length(results)) stop("No candidate model produced enough successful LOOCV predictions.")
  results_data <- do.call(rbind, results)
  results_data <- results_data[order(
    -results_data$CCC, -results_data$Cb, -results_data$rho,
    results_data$RMSE, results_data$MAE, na.last = TRUE
  ), , drop = FALSE]
  results_data$rank <- seq_len(nrow(results_data))
  results_data <- results_data[, c("rank", setdiff(names(results_data), "rank")), drop = FALSE]
  if (is.finite(top_n)) results_data <- utils::head(results_data, top_n)
  rownames(results_data) <- NULL

  prediction_data <- if (length(predictions)) do.call(rbind, predictions) else data.frame()
  failure_data <- if (length(failures)) do.call(rbind, failures) else data.frame()
  attr(results_data, "predictions") <- prediction_data
  attr(results_data, "failures") <- failure_data
  if (keep_fits) {
    retained_ids <- as.character(results_data$model_id)
    attr(results_data, "fits") <- fits_list[names(fits_list) %in% retained_ids]
  }

  results_data
}
