#' Compute lag-specific decomposition of Exposure Cumulative Impact (ECI)
#'
#' Decomposes the model-weighted Exposure Cumulative Impact into lag-specific
#' weights and contributions using numerical derivatives of the DLNM linear
#' predictor.
#'
#' Profiles must be supplied in chronological order, from the earliest to the
#' most recent observation. When `data` is used, observations are ordered
#' internally by `time`. The returned decomposition is indexed by lag, where
#' lag 0 is the most recent observation and lag max is the oldest observation.
#'
#' @param profile Numeric vector or list of numeric vectors containing exposure
#'   profiles in chronological order. For multiple variables, use a named list
#'   or an unnamed list whose order matches `var`.
#' @param data Optional long-format data frame containing the grouping column,
#'   `time`, and the exposure variables. Within each group, `time` must uniquely
#'   identify observations from earliest to most recent.
#' @param group Optional character scalar naming the grouping column in `data`.
#' @param group_level Optional grouping value or vector of values to evaluate.
#'   If `NULL`, all levels are evaluated.
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param var Optional character scalar or vector naming fitted exposure
#'   variables. Required when the fitted model contains multiple exposures.
#' @param eps Positive finite-difference perturbation. Default is `1e-6`.
#' @param center Logical. If `TRUE`, use central differences; otherwise use
#'   forward differences.
#' @param absolute Logical. If `TRUE`, include absolute contributions and
#'   absolute weights.
#' @param uncertainty Logical. If `TRUE`, propagate coefficient uncertainty.
#' @param output Character. `"summary"` or `"samples"`. This argument only
#'   changes the returned object when `uncertainty = TRUE`.
#' @param n_samples Positive integer number of coefficient draws.
#'
#' @return For one profile-variable combination, a list containing `ECI_raw`,
#'   `ECI_weighted`, and `by_lag`. With uncertainty, ECI uncertainty summaries
#'   are also returned. For multiple variables or groups, a list containing
#'   `eci_summary` and `by_lag` is returned; sample-level outputs are added when
#'   requested.
#'
#' @details
#' For a chronological profile of length `lag_max + 1`, profile position 1 is
#' associated with lag max and the final profile position is associated with
#' lag 0. The function does not reverse the profile before constructing the
#' cross-basis. The final cross-basis row performs the retrospective lag mapping.
#'
#' All weighted ECI estimates, numerical-derivative weights, and lag-specific
#' contributions are returned on the linear-predictor scale.
#'
#' Lag-specific weights are numerical derivatives of the DLNM contribution:
#' \deqn{w_l(x) \approx [\eta(x + \epsilon e_l)-\eta(x)]/\epsilon}
#' or the corresponding central difference. Contributions are calculated as
#' \deqn{contribution_l = x_l w_l(x)}.
#'
#' @export
compute_ecilag <- function(
    profile = NULL,
    data = NULL,
    group = NULL,
    group_level = NULL,
    fit,
    var = NULL,
    eps = 1e-6,
    center = TRUE,
    absolute = TRUE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {
  output <- match.arg(output)

  if (missing(fit) || is.null(fit)) stop("`fit` must be provided.")
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a positive finite numeric scalar.")
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("`center` must be TRUE or FALSE.")
  }
  if (!is.logical(absolute) || length(absolute) != 1L || is.na(absolute)) {
    stop("`absolute` must be TRUE or FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L || is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

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
    index <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    missing_index <- is.na(index)
    index[missing_index] <- seq_along(x)[missing_index]
    x[order(index)]
  }

  match_brms_draw_names <- function(cb_names_ref, draw_names) {
    prefixed <- paste0("b_", cb_names_ref)
    if (all(prefixed %in% draw_names)) return(prefixed)
    if (all(cb_names_ref %in% draw_names)) return(cb_names_ref)
    stop("Could not match brms posterior draws to cross-basis coefficients.")
  }

  get_cb_names <- function(
    coefficient_names,
    variable,
    expected_columns,
    stored_cb_columns = NULL,
    data_template = NULL
  ) {
    matched <- character(0)

    if (!is.null(stored_cb_columns)) {
      stored <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"),
        stored_cb_columns,
        value = TRUE
      ))
      matched <- stored[stored %in% coefficient_names]
    }

    if (!length(matched) && !is.null(data_template)) {
      template <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"),
        names(data_template),
        value = TRUE
      ))
      matched <- template[template %in% coefficient_names]
    }

    if (!length(matched)) {
      matched <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"),
        coefficient_names,
        value = TRUE
      ))
    }

    if (length(matched) != expected_columns) {
      stop(
        "Could not align fitted coefficients with the cross-basis for variable '",
        variable, "'. Expected ", expected_columns, " coefficients but found ",
        length(matched), "."
      )
    }

    matched
  }

  check_vcov_names <- function(covariance_matrix, coefficient_names) {
    if (is.null(rownames(covariance_matrix)) || is.null(colnames(covariance_matrix))) {
      stop("The coefficient covariance matrix must contain row and column names.")
    }
    missing_names <- setdiff(
      coefficient_names,
      intersect(rownames(covariance_matrix), colnames(covariance_matrix))
    )
    if (length(missing_names)) {
      stop(
        "The following cross-basis coefficients are absent from the covariance matrix: ",
        paste(missing_names, collapse = ", "), "."
      )
    }
    invisible(TRUE)
  }

  extract_coef_vcov <- function(model) {
    if (inherits(model, "glmmTMB")) return(list(
      beta = glmmTMB::fixef(model)$cond,
      vcov = as.matrix(stats::vcov(model)$cond)
    ))

    if (inherits(model, "merMod")) return(list(
      beta = lme4::fixef(model),
      vcov = as.matrix(stats::vcov(model))
    ))

    if (inherits(model, "lme")) return(list(
      beta = nlme::fixef(model),
      vcov = as.matrix(stats::vcov(model))
    ))

    if (inherits(model, "gls")) return(list(
      beta = stats::coef(model),
      vcov = as.matrix(stats::vcov(model))
    ))

    if (is.list(model) && !is.null(model$gam) && inherits(model$gam, "gam")) {
      return(list(
        beta = stats::coef(model$gam),
        vcov = as.matrix(stats::vcov(model$gam))
      ))
    }

    if (inherits(model, "glm") || inherits(model, "gam")) return(list(
      beta = stats::coef(model),
      vcov = as.matrix(stats::vcov(model))
    ))

    if (inherits(model, "HLfit")) {
      covariance_matrix <- tryCatch(
        as.matrix(stats::vcov(model)),
        error = function(e) NULL
      )
      if (is.null(covariance_matrix)) {
        stop("Could not extract the covariance matrix from the spaMM model.")
      }
      return(list(beta = spaMM::fixef(model), vcov = covariance_matrix))
    }

    if (inherits(model, "brmsfit")) {
      fixed_effects <- brms::fixef(model)
      beta <- fixed_effects[, "Estimate"]
      names(beta) <- rownames(fixed_effects)
      return(list(beta = beta, vcov = as.matrix(stats::vcov(model))))
    }

    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      if (is.null(names(beta))) names(beta) <- rownames(model$summary.fixed)
      covariance_matrix <- diag(model$summary.fixed$sd^2)
      dimnames(covariance_matrix) <- list(names(beta), names(beta))
      return(list(beta = beta, vcov = covariance_matrix))
    }

    if (inherits(model, "bdlnm")) {
      if (is.null(model$coefficients.summary) || is.null(model$coefficients)) {
        stop("The bdlnm model lacks coefficient summaries or posterior draws.")
      }
      beta <- model$coefficients.summary[, "mean"]
      covariance_matrix <- stats::cov(t(model$coefficients))
      return(list(beta = beta, vcov = covariance_matrix))
    }

    stop("Unsupported model class for lag-specific ECI decomposition.")
  }

  extract_beta_draws <- function(
    model,
    names_ref,
    requested_samples,
    coefficient_info
  ) {
    if (inherits(model, "brmsfit")) {
      if (!requireNamespace("posterior", quietly = TRUE)) {
        stop("Package 'posterior' is required for brms uncertainty.")
      }
      draws <- posterior::as_draws_matrix(model)
      draw_names <- match_brms_draw_names(names_ref, colnames(draws))
      draws <- as.matrix(draws[, draw_names, drop = FALSE])
      if (nrow(draws) > requested_samples) {
        draws <- draws[
          sample(seq_len(nrow(draws)), requested_samples),
          ,
          drop = FALSE
        ]
      }
      colnames(draws) <- names_ref
      return(draws)
    }

    if (inherits(model, "inla")) {
      if (!requireNamespace("INLA", quietly = TRUE)) {
        stop("Package 'INLA' is required for INLA uncertainty.")
      }
      posterior_samples <- tryCatch(
        INLA::inla.posterior.sample(n = requested_samples, result = model),
        error = function(e) NULL
      )
      if (is.null(posterior_samples)) {
        stop("Could not sample the INLA posterior; fit with config = TRUE.")
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
      if (ncol(draws) > requested_samples) {
        draws <- draws[
          ,
          sample(seq_len(ncol(draws)), requested_samples),
          drop = FALSE
        ]
      }
      output_draws <- t(draws)
      colnames(output_draws) <- names_ref
      return(output_draws)
    }

    check_vcov_names(coefficient_info$vcov, names_ref)

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty.")
    }

    draws <- MASS::mvrnorm(
      n = requested_samples,
      mu = coefficient_info$beta[names_ref],
      Sigma = coefficient_info$vcov[names_ref, names_ref, drop = FALSE]
    )
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1L)
    colnames(draws) <- names_ref
    draws
  }

  specification <- attr(fit, "epiexposure_spec")
  fitted_variables <- attr(fit, "epiexposure_vars")
  fitted_cb_columns <- attr(fit, "epiexposure_cb_cols")
  fitted_data_template <- attr(fit, "epiexposure_data_template")

  if (is.null(specification) || !is.list(specification)) {
    stop("`fit` is missing a valid `epiexposure_spec` attribute.")
  }
  if (is.null(fitted_variables) || !length(fitted_variables)) {
    stop("`fit` is missing `epiexposure_vars` metadata.")
  }

  if (is.null(var)) {
    if (length(fitted_variables) == 1L) {
      var <- fitted_variables[1]
    } else {
      stop(
        "The fitted model contains multiple exposures (",
        paste(fitted_variables, collapse = ", "),
        "). Provide `var` explicitly."
      )
    }
  }

  if (!is.character(var) || !length(var) || anyNA(var) || any(var == "")) {
    stop("`var` must be a non-empty character scalar or vector.")
  }
  if (anyDuplicated(var)) stop("`var` must contain unique names.")

  invalid_variables <- setdiff(var, fitted_variables)
  if (length(invalid_variables)) {
    stop(
      "Variables not found in the fitted model: ",
      paste(invalid_variables, collapse = ", "),
      "."
    )
  }

  using_profile <- !is.null(profile)
  using_data <- !is.null(data)
  if (using_profile && using_data) stop("Provide either `profile` or `data`, not both.")
  if (!using_profile && !using_data) stop("Provide either `profile` or `data` with `group`.")

  profiles <- NULL
  levels_to_use <- NULL

  if (using_profile) {
    if (is.numeric(profile)) {
      if (length(var) != 1L) {
        stop("A numeric `profile` can only be used with one variable.")
      }
      if (any(!is.finite(profile))) stop("`profile` must contain finite values.")
      profiles <- stats::setNames(list(as.numeric(profile)), var)
    } else if (is.list(profile)) {
      if (!length(profile)) stop("`profile` cannot be empty.")

      for (i in seq_along(profile)) {
        if (!is.numeric(profile[[i]]) || any(!is.finite(profile[[i]]))) {
          stop("All profiles must be finite numeric vectors.")
        }
        profile[[i]] <- as.numeric(profile[[i]])
      }

      if (is.null(names(profile))) {
        if (length(profile) != length(var)) {
          stop("The length of an unnamed `profile` list must match `var`.")
        }
        profiles <- profile
        names(profiles) <- var
      } else {
        if (anyNA(names(profile)) || any(names(profile) == "")) {
          stop("Every element of a named `profile` list must have a name.")
        }
        if (anyDuplicated(names(profile))) {
          stop("Names in `profile` must be unique.")
        }
        missing_profiles <- setdiff(var, names(profile))
        if (length(missing_profiles)) {
          stop("Profiles are missing for: ", paste(missing_profiles, collapse = ", "), ".")
        }
        profiles <- profile[var]
      }
    } else {
      stop("`profile` must be a numeric vector or list of numeric vectors.")
    }
  }

  if (using_data) {
    if (!is.data.frame(data)) stop("`data` must be a data.frame.")
    if (is.null(group) || !is.character(group) || length(group) != 1L ||
        is.na(group) || group == "") {
      stop("When using `data`, `group` must be one non-empty column name.")
    }

    required_columns <- c(group, "time", var)
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns)) {
      stop("Columns missing from `data`: ", paste(missing_columns, collapse = ", "), ".")
    }

    if (anyNA(data[[group]])) stop("The grouping column cannot contain missing values.")
    if (!is.numeric(data$time) || any(!is.finite(data$time))) {
      stop("`time` must contain finite numeric values.")
    }

    for (variable in var) {
      if (!is.numeric(data[[variable]]) || any(!is.finite(data[[variable]]))) {
        stop("Exposure variable '", variable, "' must contain finite numeric values.")
      }
    }

    available_levels <- unique(data[[group]])
    if (is.null(group_level)) {
      levels_to_use <- available_levels
    } else {
      missing_levels <- setdiff(group_level, available_levels)
      if (length(missing_levels)) {
        stop("Group levels not found: ", paste(missing_levels, collapse = ", "), ".")
      }
      levels_to_use <- group_level
    }

    if (!length(levels_to_use)) {
      stop("No group levels are available for lag-specific ECI calculation.")
    }
  }

  coefficient_info <- extract_coef_vcov(fit)
  beta_full <- coefficient_info$beta
  if (is.null(names(beta_full)) || anyNA(names(beta_full)) || any(names(beta_full) == "")) {
    stop("Could not determine valid coefficient names.")
  }

  build_cb_row <- function(profile_values, variable_specification, maximum_lag) {
    cross_basis <- dlnm::crossbasis(
      profile_values,
      lag = maximum_lag,
      argvar = variable_specification$argvar,
      arglag = variable_specification$arglag
    )
    row <- as.numeric(cross_basis[maximum_lag + 1L, , drop = TRUE])
    if (any(!is.finite(row))) {
      stop("The reconstructed cross-basis row contains non-finite values.")
    }
    row
  }

  calculate_eta <- function(
    profile_values,
    beta_cb,
    variable_specification,
    maximum_lag
  ) {
    row <- build_cb_row(profile_values, variable_specification, maximum_lag)
    if (length(row) != length(beta_cb)) {
      stop("Cross-basis and coefficient dimensions do not match.")
    }
    sum(row * beta_cb)
  }

  compute_one <- function(variable, chronology) {
    variable_specification <- specification[[variable]]
    if (is.null(variable_specification)) {
      stop("Missing exposure specification for variable '", variable, "'.")
    }

    maximum_lag <- as.integer(max(variable_specification$lag_max))
    expected_length <- maximum_lag + 1L

    if (length(chronology) != expected_length) {
      stop(
        "Temporal coverage for variable '", variable,
        "' must equal lag_max + 1. Expected ", expected_length,
        " observations but received ", length(chronology), "."
      )
    }
    if (!is.numeric(chronology) || any(!is.finite(chronology))) {
      stop("Profile for '", variable, "' must contain finite numeric values.")
    }

    chronology <- as.numeric(chronology)
    eci_raw <- sum(chronology)
    reference_row <- build_cb_row(chronology, variable_specification, maximum_lag)

    coefficient_names <- get_cb_names(
      coefficient_names = names(beta_full),
      variable = variable,
      expected_columns = length(reference_row),
      stored_cb_columns = fitted_cb_columns,
      data_template = fitted_data_template
    )

    lag_for_position <- maximum_lag:0

    calculate_decomposition <- function(beta_cb) {
      eta0 <- calculate_eta(
        chronology,
        beta_cb,
        variable_specification,
        maximum_lag
      )
      weights <- numeric(expected_length)

      for (i in seq_len(expected_length)) {
        plus <- chronology
        plus[i] <- plus[i] + eps

        if (center) {
          minus <- chronology
          minus[i] <- minus[i] - eps
          weights[i] <- (
            calculate_eta(plus, beta_cb, variable_specification, maximum_lag) -
              calculate_eta(minus, beta_cb, variable_specification, maximum_lag)
          ) / (2 * eps)
        } else {
          weights[i] <- (
            calculate_eta(plus, beta_cb, variable_specification, maximum_lag) - eta0
          ) / eps
        }
      }

      contribution <- chronology * weights
      denominator <- sum(abs(contribution))
      percent_contribution <- if (denominator > 0) {
        100 * abs(contribution) / denominator
      } else {
        rep(NA_real_, expected_length)
      }

      list(
        eta = eta0,
        weight = weights,
        contribution = contribution,
        percent = percent_contribution
      )
    }

    if (!uncertainty) {
      result <- calculate_decomposition(as.numeric(beta_full[coefficient_names]))
      by_lag <- data.frame(
        var = variable,
        lag = lag_for_position,
        exposure = chronology,
        weight = result$weight,
        contribution = result$contribution,
        percent_contribution = result$percent,
        stringsAsFactors = FALSE
      )

      if (absolute) {
        by_lag$abs_contribution <- abs(by_lag$contribution)
        by_lag$abs_weight <- abs(by_lag$weight)
      }

      by_lag <- by_lag[order(by_lag$lag), , drop = FALSE]
      rownames(by_lag) <- NULL

      return(list(
        ECI_raw = eci_raw,
        ECI_weighted = result$eta,
        by_lag = by_lag
      ))
    }

    beta_draws <- extract_beta_draws(
      model = fit,
      names_ref = coefficient_names,
      requested_samples = n_samples,
      coefficient_info = coefficient_info
    )

    draw_sd <- apply(beta_draws, 2, stats::sd)
    if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
      warning(
        "Near-zero coefficient-draw variability for variable '", variable,
        "'. ECI-lag intervals may collapse.",
        call. = FALSE
      )
    }

    draw_results <- lapply(seq_len(nrow(beta_draws)), function(draw_index) {
      calculate_decomposition(beta_draws[draw_index, ])
    })

    eta_draws <- vapply(draw_results, `[[`, numeric(1), "eta")
    weight_matrix <- do.call(cbind, lapply(draw_results, `[[`, "weight"))
    contribution_matrix <- do.call(cbind, lapply(draw_results, `[[`, "contribution"))
    percent_matrix <- do.call(cbind, lapply(draw_results, `[[`, "percent"))

    summarise_rows <- function(matrix_object) {
      list(
        estimate = apply(matrix_object, 1, stats::median, na.rm = TRUE),
        sd = apply(matrix_object, 1, safe_sd),
        lower = apply(matrix_object, 1, function(x) safe_quantile(x)[1]),
        upper = apply(matrix_object, 1, function(x) safe_quantile(x)[2])
      )
    }

    weight_summary <- summarise_rows(weight_matrix)
    contribution_summary <- summarise_rows(contribution_matrix)
    percent_summary <- summarise_rows(percent_matrix)

    by_lag <- data.frame(
      var = variable,
      lag = lag_for_position,
      exposure = chronology,
      weight = weight_summary$estimate,
      weight_sd = weight_summary$sd,
      weight_lower = weight_summary$lower,
      weight_upper = weight_summary$upper,
      contribution = contribution_summary$estimate,
      contribution_sd = contribution_summary$sd,
      contribution_lower = contribution_summary$lower,
      contribution_upper = contribution_summary$upper,
      percent_contribution = percent_summary$estimate,
      percent_contribution_sd = percent_summary$sd,
      percent_contribution_lower = percent_summary$lower,
      percent_contribution_upper = percent_summary$upper,
      stringsAsFactors = FALSE
    )

    if (absolute) {
      absolute_contribution_summary <- summarise_rows(abs(contribution_matrix))
      absolute_weight_summary <- summarise_rows(abs(weight_matrix))

      by_lag$abs_contribution <- absolute_contribution_summary$estimate
      by_lag$abs_contribution_sd <- absolute_contribution_summary$sd
      by_lag$abs_contribution_lower <- absolute_contribution_summary$lower
      by_lag$abs_contribution_upper <- absolute_contribution_summary$upper
      by_lag$abs_weight <- absolute_weight_summary$estimate
      by_lag$abs_weight_sd <- absolute_weight_summary$sd
      by_lag$abs_weight_lower <- absolute_weight_summary$lower
      by_lag$abs_weight_upper <- absolute_weight_summary$upper
    }

    by_lag <- by_lag[order(by_lag$lag), , drop = FALSE]
    rownames(by_lag) <- NULL
    interval <- safe_quantile(eta_draws)

    result <- list(
      ECI_raw = eci_raw,
      ECI_weighted = stats::median(eta_draws, na.rm = TRUE),
      ECI_weighted_sd = safe_sd(eta_draws),
      ECI_weighted_lower = interval[1],
      ECI_weighted_upper = interval[2],
      by_lag = by_lag
    )

    if (output == "samples") {
      sample_rows <- lapply(seq_along(draw_results), function(draw_index) {
        sample_data <- data.frame(
          sample = draw_index,
          var = variable,
          lag = lag_for_position,
          exposure = chronology,
          weight = draw_results[[draw_index]]$weight,
          contribution = draw_results[[draw_index]]$contribution,
          percent_contribution = draw_results[[draw_index]]$percent,
          stringsAsFactors = FALSE
        )
        if (absolute) {
          sample_data$abs_contribution <- abs(sample_data$contribution)
          sample_data$abs_weight <- abs(sample_data$weight)
        }
        sample_data[order(sample_data$lag), , drop = FALSE]
      })

      result$by_lag_samples <- do.call(rbind, sample_rows)
      rownames(result$by_lag_samples) <- NULL
      result$ECI_weighted_samples <- data.frame(
        sample = seq_along(eta_draws),
        ECI_weighted = eta_draws,
        stringsAsFactors = FALSE
      )
    }

    result
  }

  package_results <- function(results, variables, groups = NULL) {
    if (length(results) == 1L && is.null(groups)) return(results[[1]])

    summary_data <- do.call(rbind, lapply(seq_along(results), function(i) {
      current <- results[[i]]
      data.frame(
        var = variables[i],
        ECI_raw = current$ECI_raw,
        ECI_weighted = current$ECI_weighted,
        ECI_weighted_sd = current$ECI_weighted_sd %||% NA_real_,
        ECI_weighted_lower = current$ECI_weighted_lower %||% NA_real_,
        ECI_weighted_upper = current$ECI_weighted_upper %||% NA_real_,
        stringsAsFactors = FALSE
      )
    }))

    by_lag <- do.call(rbind, lapply(results, `[[`, "by_lag"))

    if (!is.null(groups)) {
      summary_data[[group]] <- groups
      summary_data <- summary_data[
        ,
        c(group, setdiff(names(summary_data), group)),
        drop = FALSE
      ]

      by_lag[[group]] <- rep(
        groups,
        vapply(results, function(x) nrow(x$by_lag), integer(1))
      )
      by_lag <- by_lag[, c(group, setdiff(names(by_lag), group)), drop = FALSE]
    }

    output_object <- list(
      eci_summary = summary_data,
      by_lag = by_lag
    )

    if (uncertainty && output == "samples") {
      by_lag_samples <- do.call(rbind, lapply(seq_along(results), function(i) {
        sample_data <- results[[i]]$by_lag_samples
        if (!is.null(groups)) sample_data[[group]] <- groups[i]
        sample_data
      }))

      weighted_samples <- do.call(rbind, lapply(seq_along(results), function(i) {
        sample_data <- results[[i]]$ECI_weighted_samples
        sample_data$var <- variables[i]
        if (!is.null(groups)) sample_data[[group]] <- groups[i]
        sample_data
      }))

      if (!is.null(groups)) {
        by_lag_samples <- by_lag_samples[
          ,
          c(group, setdiff(names(by_lag_samples), group)),
          drop = FALSE
        ]
        weighted_samples <- weighted_samples[
          ,
          c(group, setdiff(names(weighted_samples), group)),
          drop = FALSE
        ]
      }

      output_object$by_lag_samples <- by_lag_samples
      output_object$ECI_weighted_samples <- weighted_samples
    }

    output_object
  }

  if (using_profile) {
    results <- lapply(var, function(variable) {
      compute_one(variable, profiles[[variable]])
    })
    return(package_results(results, var))
  }

  results <- list()
  result_variables <- character(0)
  result_groups <- rep(levels_to_use, each = length(var))
  output_index <- 1L

  for (current_level in levels_to_use) {
    group_data <- data[data[[group]] == current_level, , drop = FALSE]
    if (!nrow(group_data)) {
      stop("No observations found for group level '", current_level, "'.")
    }

    group_data <- group_data[order(group_data$time), , drop = FALSE]

    if (anyDuplicated(group_data$time)) {
      stop("Duplicated `time` values for group level '", current_level, "'.")
    }
    if (is.unsorted(group_data$time, strictly = TRUE)) {
      stop("`time` is not strictly increasing for group level '", current_level, "'.")
    }

    for (variable in var) {
      results[[output_index]] <- compute_one(
        variable,
        as.numeric(group_data[[variable]])
      )
      result_variables[output_index] <- variable
      output_index <- output_index + 1L
    }
  }

  package_results(results, result_variables, result_groups)
}
