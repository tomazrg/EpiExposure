#' Summarise DLNM effects
#'
#' Summarises DLNM effects on either the daily or accumulated scale using
#' exposure grids derived from observed data and the DLNM specification stored
#' in the fitted model.
#'
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param data Long-format exposure data. Must contain `epi_id`, `time`, and
#'   the exposure variables used in the fitted model. Within each epidemic,
#'   observations must be supplied chronologically, from earliest to most recent.
#' @param var Exposure-variable name, character vector of names, or `NULL`.
#'   If `NULL`, all variables stored in the fitted-model metadata are used.
#' @param lag_max Optional maximum lag. Ignored when available in model metadata.
#' @param df_var Optional exposure-basis degrees of freedom used as fallback.
#' @param df_lag Optional lag-basis degrees of freedom used as fallback.
#' @param fun_var Optional exposure-basis function used as fallback.
#' @param fun_lag Optional lag-basis function used as fallback.
#' @param scale Character. `"daily"` or `"accumulated"`.
#' @param lag_periods Optional table containing `period`, `lag_start`, and
#'   `lag_end`; required for non-incremental accumulated summaries.
#' @param probs Quantiles used to define the exposure grid.
#' @param ref Reference definition list with `method` and optional `value`.
#' @param effect_measure Character. `"linear"`, `"exponentiated"`, or
#'   `"percent"`.
#' @param incremental Logical. If `TRUE` with accumulated effects, returns
#'   cumulative effects from lag 0 through each lag.
#' @param uncertainty Logical. If `TRUE`, propagates coefficient uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Positive integer number of coefficient draws.
#'
#' @details
#' `time` is chronological, not lag. Observations are ordered internally by
#' `time` before reconstructing the cross-basis. The DLNM representation maps
#' the most recent observation to lag 0. Returned effects remain indexed by lag.
#'
#' Coefficient names are aligned with cross-basis columns using the ordered
#' `epiexposure_cb_cols` metadata stored by `fit_epidlnm()`. The function checks
#' both the coefficient vector and covariance-matrix names before prediction.
#'
#' @return A data.frame.
#' @export
summarise_effects <- function(
    fit,
    data,
    var = NULL,
    lag_max = NULL,
    df_var = NULL,
    df_lag = NULL,
    fun_var = NULL,
    fun_lag = NULL,
    scale = c("daily", "accumulated"),
    lag_periods = NULL,
    probs = seq(0.05, 0.95, by = 0.01),
    ref = list(method = "median", value = NULL),
    effect_measure = c("linear", "exponentiated", "percent"),
    incremental = FALSE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {
  scale <- match.arg(scale)
  effect_measure <- match.arg(effect_measure)
  output <- match.arg(output)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!all(c("epi_id", "time") %in% names(data))) {
    stop("`data` must contain at least 'epi_id' and 'time'.")
  }
  if (!is.logical(incremental) || length(incremental) != 1L || is.na(incremental)) {
    stop("`incremental` must be TRUE or FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L || is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1)) {
    stop("`probs` must contain finite probabilities between 0 and 1.")
  }
  if (!is.list(ref) || is.null(ref$method) || length(ref$method) != 1L) {
    stop("`ref` must be a list containing a single `method` value.")
  }
  if (scale == "accumulated" && !incremental && is.null(lag_periods)) {
    stop("`lag_periods` is required when scale = 'accumulated' and incremental = FALSE.")
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b
  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    x <- x[is.finite(x)]
    if (!length(x)) return(rep(NA_real_, length(probs)))
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }
  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(0)
    stats::sd(x)
  }
  sort_cb_names <- function(x) {
    if (!length(x)) return(x)
    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    missing_idx <- is.na(idx)
    idx[missing_idx] <- seq_along(x)[missing_idx]
    x[order(idx)]
  }
  match_brms_draw_names <- function(cb_names_ref, draw_names) {
    prefixed <- paste0("b_", cb_names_ref)
    if (all(prefixed %in% draw_names)) return(prefixed)
    if (all(cb_names_ref %in% draw_names)) return(cb_names_ref)
    stop("Could not match brms posterior draws to cross-basis coefficients.")
  }
  get_cb_names <- function(coefficient_names, variable, expected, stored = NULL) {
    matched <- character(0)
    if (!is.null(stored)) {
      ordered <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"), stored, value = TRUE
      ))
      matched <- ordered[ordered %in% coefficient_names]
    }
    if (!length(matched)) {
      matched <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"), coefficient_names, value = TRUE
      ))
    }
    if (length(matched) != expected) {
      stop(
        "Could not align coefficients with the cross-basis for variable '",
        variable, "'. Expected ", expected, " coefficients but found ",
        length(matched), "."
      )
    }
    matched
  }
  check_vcov_names <- function(V, coefficient_names) {
    if (is.null(rownames(V)) || is.null(colnames(V))) {
      stop("The coefficient covariance matrix has no row or column names.")
    }
    missing_names <- setdiff(
      coefficient_names,
      intersect(rownames(V), colnames(V))
    )
    if (length(missing_names)) {
      stop(
        "Cross-basis coefficients missing from the covariance matrix: ",
        paste(missing_names, collapse = ", "), "."
      )
    }
    invisible(TRUE)
  }
  check_lag_coverage <- function(dat, maximum_lag) {
    required_n <- maximum_lag + 1L
    temporal_summary <- dat |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        n_time = dplyr::n_distinct(time),
        has_missing_time = any(is.na(time)),
        .groups = "drop"
      )
    missing_time <- temporal_summary |> dplyr::filter(has_missing_time)
    if (nrow(missing_time)) {
      stop("Missing `time` values detected. Example epi_id: ",
           paste(utils::head(missing_time$epi_id, 5), collapse = ", "), ".")
    }
    duplicate_time <- temporal_summary |> dplyr::filter(n_rows != n_time)
    if (nrow(duplicate_time)) {
      stop("Duplicated `time` values detected within epidemics. Example epi_id: ",
           paste(utils::head(duplicate_time$epi_id, 5), collapse = ", "), ".")
    }
    insufficient <- temporal_summary |> dplyr::filter(n_time < required_n)
    if (nrow(insufficient)) {
      stop(
        "Insufficient temporal coverage. Required observations per epi_id: ",
        required_n, ". Example epi_id: ",
        paste(utils::head(insufficient$epi_id, 5), collapse = ", "), "."
      )
    }
    invisible(TRUE)
  }
  build_pooled_series <- function(dat, variable, separator_n) {
    ids <- unique(dat$epi_id)
    pooled <- vector("list", length(ids))
    for (i in seq_along(ids)) {
      values <- dat |>
        dplyr::filter(epi_id == ids[i]) |>
        dplyr::arrange(time) |>
        dplyr::pull(.data[[variable]])
      pooled[[i]] <- c(values, rep(NA_real_, separator_n))
    }
    unlist(pooled)
  }

  extract_coef_vcov <- function(model) {
    if (inherits(model, "glmmTMB")) return(list(
      beta = glmmTMB::fixef(model)$cond,
      vcov = as.matrix(stats::vcov(model)$cond)
    ))
    if (inherits(model, "merMod")) return(list(
      beta = lme4::fixef(model), vcov = as.matrix(stats::vcov(model))
    ))
    if (inherits(model, "lme")) return(list(
      beta = nlme::fixef(model), vcov = as.matrix(stats::vcov(model))
    ))
    if (inherits(model, "gls") || inherits(model, "gam")) return(list(
      beta = stats::coef(model), vcov = as.matrix(stats::vcov(model))
    ))
    if (inherits(model, "HLfit")) {
      V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
      if (is.null(V)) stop("Could not extract vcov from the spaMM model.")
      return(list(beta = spaMM::fixef(model), vcov = V))
    }
    if (inherits(model, "brmsfit")) {
      fixed <- brms::fixef(model)
      beta <- fixed[, "Estimate"]
      names(beta) <- rownames(fixed)
      return(list(beta = beta, vcov = as.matrix(stats::vcov(model))))
    }
    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      if (is.null(names(beta))) names(beta) <- rownames(model$summary.fixed)
      V <- diag(model$summary.fixed$sd^2)
      dimnames(V) <- list(names(beta), names(beta))
      return(list(beta = beta, vcov = V))
    }
    if (inherits(model, "bdlnm")) {
      beta <- model$coefficients.summary[, "mean"]
      return(list(beta = beta, vcov = stats::cov(t(model$coefficients))))
    }
    beta <- stats::coef(model)
    if (!is.numeric(beta)) stop("Could not extract numeric model coefficients.")
    list(beta = beta, vcov = as.matrix(stats::vcov(model)))
  }

  extract_beta_draws <- function(model, names_ref, n, coefficient_info) {
    if (inherits(model, "brmsfit")) {
      if (!requireNamespace("posterior", quietly = TRUE)) {
        stop("Package 'posterior' is required for brms uncertainty.")
      }
      draws <- posterior::as_draws_matrix(model)
      draw_names <- match_brms_draw_names(names_ref, colnames(draws))
      draws <- as.matrix(draws[, draw_names, drop = FALSE])
      if (nrow(draws) > n) {
        draws <- draws[sample(seq_len(nrow(draws)), n), , drop = FALSE]
      }
      colnames(draws) <- names_ref
      return(draws)
    }
    if (inherits(model, "inla")) {
      if (!requireNamespace("INLA", quietly = TRUE)) {
        stop("Package 'INLA' is required for INLA uncertainty.")
      }
      posterior_samples <- tryCatch(
        INLA::inla.posterior.sample(n = n, result = model),
        error = function(e) NULL
      )
      if (is.null(posterior_samples)) {
        stop("Could not draw INLA posterior samples; fit with config = TRUE.")
      }
      draws <- do.call(rbind, lapply(posterior_samples, function(sample_object) {
        latent <- sample_object$latent
        names(latent) <- gsub(":1$", "", names(latent))
        values <- latent[names_ref]
        if (anyNA(values)) {
          stop("Could not match INLA draws to cross-basis coefficients.")
        }
        as.numeric(values)
      }))
      colnames(draws) <- names_ref
      return(draws)
    }
    if (inherits(model, "bdlnm")) {
      draws <- model$coefficients
      if (is.null(dim(draws))) draws <- matrix(draws, ncol = 1L)
      if (is.null(rownames(draws)) || !all(names_ref %in% rownames(draws))) {
        stop("Could not match bdlnm draws to cross-basis coefficients.")
      }
      draws <- draws[names_ref, , drop = FALSE]
      if (ncol(draws) > n) {
        draws <- draws[, sample(seq_len(ncol(draws)), n), drop = FALSE]
      }
      draws <- t(draws)
      colnames(draws) <- names_ref
      return(draws)
    }
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty.")
    }
    check_vcov_names(coefficient_info$vcov, names_ref)
    draws <- MASS::mvrnorm(
      n = n,
      mu = coefficient_info$beta[names_ref],
      Sigma = coefficient_info$vcov[names_ref, names_ref, drop = FALSE]
    )
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1L)
    colnames(draws) <- names_ref
    draws
  }

  get_linkfun <- function(family_object) {
    if (inherits(family_object, "family") && !is.null(family_object$linkfun)) {
      return(family_object$linkfun)
    }
    if (is.character(family_object)) {
      family_name <- tolower(family_object[1])
      if (family_name %in% c("beta", "binomial")) return(stats::qlogis)
      if (family_name %in% c("poisson", "gamma", "negbin", "negative_binomial")) return(log)
      if (family_name == "gaussian") return(identity)
    }
    stop("Could not determine the model link function from metadata.")
  }
  get_linkinv <- function(family_object) {
    if (inherits(family_object, "family") && !is.null(family_object$linkinv)) {
      return(family_object$linkinv)
    }
    if (is.character(family_object)) {
      family_name <- tolower(family_object[1])
      if (family_name %in% c("beta", "binomial")) return(stats::plogis)
      if (family_name %in% c("poisson", "gamma", "negbin", "negative_binomial")) return(exp)
      if (family_name == "gaussian") return(identity)
    }
    stop("Could not determine the model inverse-link function from metadata.")
  }

  fit_spec <- attr(fit, "epiexposure_spec")
  fit_vars <- attr(fit, "epiexposure_vars")
  family_fit <- attr(fit, "epiexposure_family")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")
  data_template <- attr(fit, "epiexposure_data_template")

  if (is.null(var)) {
    variables <- fit_vars
    if (is.null(variables) || !length(variables)) {
      stop("`var` is NULL and no exposure variables were found in model metadata.")
    }
  } else {
    if (!is.character(var) || !length(var) || anyNA(var) || any(var == "")) {
      stop("`var` must be NULL or a non-empty character vector.")
    }
    variables <- var
  }
  if (anyDuplicated(variables)) stop("`var` must contain unique variable names.")
  unknown_model_vars <- if (is.null(fit_vars)) character(0) else setdiff(variables, fit_vars)
  if (length(unknown_model_vars)) {
    stop("Variables not found in fitted-model metadata: ",
         paste(unknown_model_vars, collapse = ", "), ".")
  }
  missing_data_vars <- setdiff(variables, names(data))
  if (length(missing_data_vars)) {
    stop("Variables not present in `data`: ",
         paste(missing_data_vars, collapse = ", "), ".")
  }

  transform_effect <- function(eta) {
    if (effect_measure == "linear") return(eta)
    if (effect_measure == "exponentiated") return(exp(eta))
    (exp(eta) - 1) * 100
  }

  summarise_one_variable <- function(variable) {
    if (!is.null(fit_spec) && !is.null(fit_spec[[variable]])) {
      variable_spec <- fit_spec[[variable]]
      lag_max_use <- as.integer(max(variable_spec$lag_max))
      argvar <- variable_spec$argvar
      arglag <- variable_spec$arglag
    } else {
      if (is.null(lag_max)) {
        stop("No stored specification for '", variable, "'; provide `lag_max`.")
      }
      lag_max_use <- as.integer(max(lag_max))
      exposure_function <- fun_var %||% "ns"
      lag_function <- fun_lag %||% "ns"
      exposure_df <- df_var %||% 4
      lag_df <- df_lag %||% 4
      argvar <- switch(
        exposure_function,
        ns = list(fun = "ns", df = exposure_df),
        bs = list(fun = "bs", df = exposure_df),
        poly = list(fun = "poly", degree = exposure_df),
        lin = list(fun = "lin"),
        stop("Unsupported fun_var: ", exposure_function)
      )
      if (argvar$fun != "lin") argvar$intercept <- FALSE
      arglag <- switch(
        lag_function,
        ns = list(fun = "ns", df = lag_df),
        ps = list(fun = "ps", df = lag_df),
        lin = list(fun = "lin"),
        stop("Unsupported fun_lag: ", lag_function)
      )
    }

    check_lag_coverage(data, lag_max_use)
    pooled <- build_pooled_series(data, variable, lag_max_use)
    cb <- dlnm::crossbasis(
      pooled,
      lag = lag_max_use,
      argvar = argvar,
      arglag = arglag
    )

    coefficient_info <- extract_coef_vcov(fit)
    coefficients <- coefficient_info$beta
    covariance <- coefficient_info$vcov
    if (is.null(names(coefficients))) {
      stop("The fitted coefficient vector has no names.")
    }

    stored_columns <- cb_cols_fit
    if (is.null(stored_columns) && !is.null(data_template)) {
      stored_columns <- grep("^cb_", names(data_template), value = TRUE)
    }
    cb_names_ref <- get_cb_names(
      coefficient_names = names(coefficients),
      variable = variable,
      expected = ncol(cb),
      stored = stored_columns
    )
    check_vcov_names(covariance, cb_names_ref)
    beta <- coefficients[cb_names_ref]
    vcov_cb <- covariance[cb_names_ref, cb_names_ref, drop = FALSE]

    x_all <- data[[variable]]
    if (!is.numeric(x_all) || all(!is.finite(x_all))) {
      stop("Exposure variable '", variable, "' must contain finite numeric values.")
    }
    at_values <- sort(unique(as.numeric(stats::quantile(
      x_all, probs = probs, na.rm = TRUE
    ))))
    center_value <- switch(
      ref$method,
      median = stats::median(x_all, na.rm = TRUE),
      percentile = {
        if (is.null(ref$value) || length(ref$value) != 1L ||
            !is.numeric(ref$value) || ref$value < 0 || ref$value > 1) {
          stop("For ref$method = 'percentile', `ref$value` must be one probability.")
        }
        stats::quantile(x_all, ref$value, na.rm = TRUE)
      },
      fixed = {
        if (is.null(ref$value) || length(ref$value) != 1L ||
            !is.numeric(ref$value) || !is.finite(ref$value)) {
          stop("For ref$method = 'fixed', `ref$value` must be one finite number.")
        }
        ref$value
      },
      stop("Invalid ref$method. Use 'median', 'percentile', or 'fixed'.")
    )
    center_value <- as.numeric(center_value)

    point_prediction <- dlnm::crosspred(
      cb,
      coef = beta,
      vcov = vcov_cb,
      at = at_values,
      cen = center_value,
      bylag = 1
    )
    point_eta <- point_prediction$matfit
    lag_index <- suppressWarnings(as.integer(gsub("lag", "", colnames(point_eta))))
    if (anyNA(lag_index)) lag_index <- 0:(ncol(point_eta) - 1L)
    lag_order <- order(lag_index)
    lag_index <- lag_index[lag_order]
    point_eta <- point_eta[, lag_order, drop = FALSE]

    linkfun <- get_linkfun(family_fit)
    linkinv <- get_linkinv(family_fit)
    baseline_response <- NA_real_
    if (!is.null(fit_vars) && !is.null(fit_spec) &&
        all(fit_vars %in% names(data))) {
      reference_profiles <- lapply(fit_vars, function(current_variable) {
        current_spec <- fit_spec[[current_variable]]
        current_lag <- as.integer(max(current_spec$lag_max))
        reference_value <- stats::median(data[[current_variable]], na.rm = TRUE)
        rep(reference_value, current_lag + 1L)
      })
      names(reference_profiles) <- fit_vars
      baseline_response <- tryCatch(
        as.numeric(predict_outcome(
          fit = fit,
          profiles = reference_profiles,
          re = "population",
          type = "response",
          uncertainty = FALSE
        )$prediction[1]),
        error = function(e) {
          warning("Baseline prediction failed: ", conditionMessage(e), call. = FALSE)
          NA_real_
        }
      )
    }

    build_output <- function(eta_matrix) {
      if (scale == "daily") {
        transformed_eta <- eta_matrix
        grid <- expand.grid(
          value = at_values,
          lag = lag_index,
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )
        result <- data.frame(
          var = variable,
          lag = grid$lag,
          scale = "daily",
          value = grid$value,
          eta = as.vector(transformed_eta),
          stringsAsFactors = FALSE
        )
      } else if (incremental) {
        transformed_eta <- t(apply(eta_matrix, 1, cumsum))
        grid <- expand.grid(
          value = at_values,
          lag = lag_index,
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )
        result <- data.frame(
          var = variable,
          lag = grid$lag,
          period = paste0("0-", grid$lag),
          scale = "accumulated",
          value = grid$value,
          eta = as.vector(transformed_eta),
          stringsAsFactors = FALSE
        )
      } else {
        required_period_columns <- c("period", "lag_start", "lag_end")
        if (!all(required_period_columns %in% names(lag_periods))) {
          stop("`lag_periods` must contain period, lag_start, and lag_end.")
        }
        return(purrr::map_dfr(seq_len(nrow(lag_periods)), function(i) {
          selected_lags <- which(
            lag_index >= lag_periods$lag_start[i] &
              lag_index <= lag_periods$lag_end[i]
          )
          if (!length(selected_lags)) {
            stop(
              "No lags found for period '", lag_periods$period[i],
              "' (", lag_periods$lag_start[i], "-", lag_periods$lag_end[i], ")."
            )
          }
          period_eta <- rowSums(eta_matrix[, selected_lags, drop = FALSE])
          period_result <- data.frame(
            var = variable,
            period = lag_periods$period[i],
            scale = "accumulated",
            value = at_values,
            eta = period_eta,
            stringsAsFactors = FALSE
          )
          period_result$effect <- transform_effect(period_result$eta)
          period_result$baseline <- baseline_response
          if (is.na(baseline_response)) {
            period_result$predicted <- NA_real_
            period_result$delta <- NA_real_
            period_result$delta_pp <- NA_real_
          } else {
            period_result$predicted <- linkinv(
              linkfun(baseline_response) + period_result$eta
            )
            period_result$delta <- period_result$predicted - baseline_response
            period_result$delta_pp <- 100 * period_result$delta
          }
          period_result
        }))
      }

      result$effect <- transform_effect(result$eta)
      result$baseline <- baseline_response
      if (is.na(baseline_response)) {
        result$predicted <- NA_real_
        result$delta <- NA_real_
        result$delta_pp <- NA_real_
      } else {
        result$predicted <- linkinv(linkfun(baseline_response) + result$eta)
        result$delta <- result$predicted - baseline_response
        result$delta_pp <- 100 * result$delta
      }
      result
    }

    if (!uncertainty) return(build_output(point_eta))

    beta_draws <- extract_beta_draws(
      model = fit,
      names_ref = cb_names_ref,
      n = n_samples,
      coefficient_info = coefficient_info
    )
    draw_sd <- apply(beta_draws, 2, stats::sd)
    if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
      warning("Near-zero coefficient-draw variability for variable '",
              variable, "'. Intervals may collapse.", call. = FALSE)
    }

    eta_list <- lapply(seq_len(nrow(beta_draws)), function(i) {
      draw_prediction <- dlnm::crosspred(
        cb,
        coef = beta_draws[i, ],
        vcov = diag(0, ncol(beta_draws)),
        at = at_values,
        cen = center_value,
        bylag = 1
      )
      draw_prediction$matfit[, lag_order, drop = FALSE]
    })
    sample_output <- dplyr::bind_rows(lapply(seq_along(eta_list), function(i) {
      current <- build_output(eta_list[[i]])
      current$sample <- i
      current[, c("sample", setdiff(names(current), "sample")), drop = FALSE]
    }))
    if (output == "samples") return(sample_output)

    grouping_columns <- intersect(
      c("var", "lag", "period", "scale", "value"),
      names(sample_output)
    )
    metric_columns <- intersect(
      c("eta", "effect", "delta", "delta_pp", "predicted"),
      names(sample_output)
    )
    summary_output <- sample_output |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(metric_columns),
          list(
            estimate = ~stats::median(.x, na.rm = TRUE),
            sd = ~safe_sd(.x),
            lower = ~safe_quantile(.x)[1],
            upper = ~safe_quantile(.x)[2]
          ),
          .names = "{.col}_{.fn}"
        ),
        baseline = dplyr::first(.data$baseline),
        .groups = "drop"
      )
    for (metric in metric_columns) {
      estimate_column <- paste0(metric, "_estimate")
      names(summary_output)[names(summary_output) == estimate_column] <- metric
    }
    as.data.frame(summary_output)
  }

  output_list <- lapply(variables, summarise_one_variable)
  output_data <- dplyr::bind_rows(output_list)
  rownames(output_data) <- NULL
  output_data
}
