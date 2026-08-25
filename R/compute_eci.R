#' Compute Exposure Cumulative Impact
#'
#' Computes the raw Exposure Cumulative Impact (`ECI_raw`) and, when a fitted
#' model is provided, the model-weighted Exposure Cumulative Impact
#' (`ECI_weighted`) for one or more exposure variables.
#'
#' Exposure profiles must always be supplied in chronological order, from the
#' earliest observation to the most recent observation before disease
#' assessment. The function preserves this ordering when reconstructing the
#' cross-basis. In the final cross-basis row, the most recent observation is
#' internally associated with lag 0.
#'
#' @param profile Numeric vector or list of numeric vectors representing one or
#'   more chronological exposure profiles. If `fit = NULL`, `profile` must be a
#'   single numeric vector. For multiple fitted variables, use a named list or
#'   an unnamed list whose order matches `var`.
#' @param data Optional long-format data frame containing observed exposure
#'   histories. Must contain the grouping column specified in `group`, `time`,
#'   and all exposure variables requested in `var`. Observations are ordered
#'   internally by `time` within each group.
#' @param group Optional character scalar naming the grouping column in `data`.
#' @param group_level Optional grouping value or vector of values to evaluate.
#'   If `NULL`, all observed levels are evaluated.
#' @param fit Fitted model returned by `fit_epidlnm()`. If `NULL`, only
#'   `ECI_raw` is calculated.
#' @param var Optional character scalar or vector naming fitted exposure
#'   variables. Required when the fitted model contains multiple exposures.
#' @param scale Character. Scale used for `ECI_weighted`: `"link"`,
#'   `"response"`, or `"percent"`.
#' @param uncertainty Logical. If `TRUE`, propagates coefficient uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Positive integer number of coefficient draws.
#'
#' @return A data.frame containing `ECI_raw`, `ECI_weighted`, and `scale`.
#'   Depending on the requested input and uncertainty settings, the output may
#'   also contain `var`, the grouping column, `sample`, `sd`, `lower`, and
#'   `upper`.
#'
#' @details
#' `ECI_raw` is the unweighted sum of the chronological exposure profile:
#' \deqn{ECI_{raw}=\sum_{t=0}^{T}x_t}
#'
#' `ECI_weighted` is the contribution of the selected exposure on the requested
#' scale, calculated from the final cross-basis row and the aligned fitted
#' cross-basis coefficients:
#' \deqn{ECI_{weighted}=\mathbf{cb}(x)^\top\boldsymbol{\beta}}
#'
#' For `scale = "response"`, the inverse link is applied to the exposure
#' contribution only. The model intercept, other predictors, and random effects
#' are not included.
#'
#' For `scale = "percent"`, the transformation is:
#' `(exp(eta) - 1) * 100`.
#'
#' For models using a log link (e.g. Poisson, Gamma and Negative Binomial),
#' this quantity can be interpreted as a relative percentage change with
#' respect to the reference exposure condition.
#'
#' For models using a logit link (e.g. Beta and Binomial), this quantity
#' represents a relative exposure effect on the linear-predictor scale and
#' should not be interpreted as a direct percentage change in the response
#' variable.
#'
#' When `uncertainty = FALSE`, ECI is calculated using the central coefficient
#' estimates provided by the fitted model.
#'
#' When `uncertainty = TRUE`, ECI is calculated for each coefficient draw and
#' summarized using the median together with empirical uncertainty intervals.
#'
#' Frequentist uncertainty is approximated with draws from the asymptotic
#' multivariate normal coefficient distribution. Bayesian engines use posterior
#' coefficient draws when available.
#'
#' @export
compute_eci <- function(
    profile = NULL,
    data = NULL,
    group = NULL,
    group_level = NULL,
    fit = NULL,
    var = NULL,
    scale = c("link", "response", "percent"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {
  scale <- match.arg(scale)
  output <- match.arg(output)

  if (!is.logical(uncertainty) || length(uncertainty) != 1L || is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }

  n_samples <- as.integer(n_samples)

  if (uncertainty && n_samples < 2L) {
    stop(
      "`n_samples` must be at least 2 when uncertainty = TRUE."
    )
  }

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

  match_brms_draw_names <- function(cb_names_ref, draw_colnames) {
    prefixed_names <- paste0("b_", cb_names_ref)
    if (all(prefixed_names %in% draw_colnames)) return(prefixed_names)
    if (all(cb_names_ref %in% draw_colnames)) return(cb_names_ref)
    stop("Could not match brms posterior draws to cross-basis coefficients.")
  }

  get_cb_names <- function(
    coefficient_names,
    var_sel,
    expected_columns,
    cb_cols_fit = NULL,
    data_template = NULL
  ) {
    cb_names <- character(0)

    if (!is.null(cb_cols_fit)) {
      stored_names <- sort_cb_names(grep(
        paste0("^cb_", var_sel, "_"),
        cb_cols_fit,
        value = TRUE
      ))
      cb_names <- stored_names[stored_names %in% coefficient_names]
    }

    if (!length(cb_names) && !is.null(data_template)) {
      template_names <- sort_cb_names(grep(
        paste0("^cb_", var_sel, "_"),
        names(data_template),
        value = TRUE
      ))
      cb_names <- template_names[template_names %in% coefficient_names]
    }

    if (!length(cb_names)) {
      cb_names <- sort_cb_names(grep(
        paste0("^cb_", var_sel, "_"),
        coefficient_names,
        value = TRUE
      ))
    }

    if (length(cb_names) != expected_columns) {
      stop(
        "Could not align fitted coefficients with the cross-basis for variable '",
        var_sel, "'. Expected ", expected_columns, " coefficients but found ",
        length(cb_names), "."
      )
    }

    cb_names
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

  get_linkinv <- function(family_fit) {
    if (inherits(family_fit, "family") && !is.null(family_fit$linkinv)) {
      return(family_fit$linkinv)
    }

    if (is.character(family_fit)) {
      family_name <- tolower(family_fit[1])
      if (family_name %in% c("beta", "binomial")) return(stats::plogis)
      if (family_name %in% c(
        "poisson", "gamma", "negbin", "negative_binomial"
      )) return(exp)
      if (family_name == "gaussian") return(identity)
    }

    stop("Could not determine the inverse-link function from model metadata.")
  }

  transform_scale <- function(x, requested_scale, linkinv) {
    if (requested_scale == "link") return(x)
    if (requested_scale == "response") return(linkinv(x))
    if (requested_scale == "percent") return((exp(x) - 1) * 100)
    stop("Unsupported ECI scale.")
  }

  extract_coef_vcov <- function(model) {
    if (inherits(model, "glmmTMB")) {
      return(list(
        beta = glmmTMB::fixef(model)$cond,
        vcov = as.matrix(stats::vcov(model)$cond)
      ))
    }

    if (inherits(model, "merMod")) {
      return(list(
        beta = lme4::fixef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "lme")) {
      return(list(
        beta = nlme::fixef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "gls")) {
      return(list(
        beta = stats::coef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (is.list(model) && !is.null(model$gam) && inherits(model$gam, "gam")) {
      return(list(
        beta = stats::coef(model$gam),
        vcov = as.matrix(stats::vcov(model$gam))
      ))
    }

    if (inherits(model, "gam")) {
      return(list(
        beta = stats::coef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "HLfit")) {
      covariance_matrix <- tryCatch(
        as.matrix(stats::vcov(model)),
        error = function(e) NULL
      )
      if (is.null(covariance_matrix)) {
        stop("Could not extract the covariance matrix from the spaMM model.")
      }
      return(list(
        beta = spaMM::fixef(model),
        vcov = covariance_matrix
      ))
    }

    if (inherits(model, "brmsfit")) {
      fixed_effects <- brms::fixef(model)
      beta <- fixed_effects[, "Estimate"]
      names(beta) <- rownames(fixed_effects)
      return(list(
        beta = beta,
        vcov = as.matrix(stats::vcov(model))
      ))
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

    beta <- stats::coef(model)
    if (!is.numeric(beta)) {
      stop("Could not extract numeric coefficients from the fitted model.")
    }
    covariance_matrix <- tryCatch(
      as.matrix(stats::vcov(model)),
      error = function(e) NULL
    )
    if (is.null(covariance_matrix)) {
      stop("Could not extract the coefficient covariance matrix.")
    }
    list(beta = beta, vcov = covariance_matrix)
  }

  extract_beta_draws <- function(
    model,
    cb_names_ref,
    requested_samples,
    coefficient_info
  ) {
    if (inherits(model, "brmsfit")) {
      if (!requireNamespace("posterior", quietly = TRUE)) {
        stop("Package 'posterior' is required for brms uncertainty.")
      }
      draws <- posterior::as_draws_matrix(model)
      draw_names <- match_brms_draw_names(cb_names_ref, colnames(draws))
      draws <- as.matrix(draws[, draw_names, drop = FALSE])
      if (nrow(draws) > requested_samples) {
        keep <- sample(seq_len(nrow(draws)), requested_samples, replace = FALSE)
        draws <- draws[keep, , drop = FALSE]
      }
      colnames(draws) <- cb_names_ref
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
        stop(
          "INLA posterior samples could not be generated. Fit the model with ",
          "control.compute = list(config = TRUE)."
        )
      }
      draws <- do.call(rbind, lapply(posterior_samples, function(sample_object) {
        latent <- sample_object$latent
        names(latent) <- gsub(":1$", "", names(latent))
        values <- latent[cb_names_ref]
        if (anyNA(values)) {
          stop("Could not match INLA draws to cross-basis coefficients.")
        }
        as.numeric(values)
      }))
      colnames(draws) <- cb_names_ref
      return(draws)
    }

    if (inherits(model, "bdlnm")) {
      draws <- model$coefficients
      if (is.null(dim(draws))) draws <- matrix(draws, ncol = 1L)
      if (is.null(rownames(draws)) || !all(cb_names_ref %in% rownames(draws))) {
        stop("Could not match bdlnm draws to cross-basis coefficients.")
      }
      draws <- draws[cb_names_ref, , drop = FALSE]
      if (ncol(draws) > requested_samples) {
        keep <- sample(seq_len(ncol(draws)), requested_samples, replace = FALSE)
        draws <- draws[, keep, drop = FALSE]
      }
      output_draws <- t(draws)
      colnames(output_draws) <- cb_names_ref
      return(output_draws)
    }

    check_vcov_names(
      coefficient_info$vcov,
      cb_names_ref
    )

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty.")
    }

    beta_hat <- coefficient_info$beta[cb_names_ref]
    covariance_matrix <- coefficient_info$vcov[
      cb_names_ref,
      cb_names_ref,
      drop = FALSE
    ]

    draws <- MASS::mvrnorm(
      n = requested_samples,
      mu = beta_hat,
      Sigma = covariance_matrix
    )
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1L)
    colnames(draws) <- cb_names_ref
    draws
  }

  # ------------------------------------------------------------
  # Raw ECI without a fitted model
  # ------------------------------------------------------------

  if (is.null(fit)) {
    if (is.null(profile)) {
      stop("When `fit = NULL`, `profile` must be provided.")
    }
    if (!is.numeric(profile) || any(!is.finite(profile))) {
      stop("When `fit = NULL`, `profile` must be a finite numeric vector.")
    }
    if (!is.null(data)) {
      warning(
        "`data`, `group`, and `group_level` are ignored when `fit = NULL`.",
        call. = FALSE
      )
    }
    if (!is.null(var)) {
      warning("`var` is ignored when `fit = NULL`.", call. = FALSE)
    }
    return(data.frame(
      ECI_raw = sum(profile),
      ECI_weighted = NA_real_,
      scale = scale,
      stringsAsFactors = FALSE
    ))
  }

  # ------------------------------------------------------------
  # Model metadata
  # ------------------------------------------------------------

  specification <- attr(fit, "epiexposure_spec")
  fitted_variables <- attr(fit, "epiexposure_vars")
  fitted_cb_columns <- attr(fit, "epiexposure_cb_cols")
  fitted_data_template <- attr(fit, "epiexposure_data_template")
  fitted_family <- attr(fit, "epiexposure_family")

  if (is.null(specification) || !is.list(specification)) {
    stop("`fit` does not contain a valid `epiexposure_spec` attribute.")
  }
  if (is.null(fitted_variables) || !length(fitted_variables)) {
    stop("`fit` does not contain exposure-variable metadata.")
  }
  if (is.null(fitted_family)) {
    stop("`fit` does not contain `epiexposure_family` metadata.")
  }

  # ------------------------------------------------------------
  # Resolve variables
  # ------------------------------------------------------------

  if (is.null(var)) {
    if (length(fitted_variables) == 1L) {
      var <- fitted_variables[1]
    } else {
      stop(
        "The fitted model contains multiple exposure variables: ",
        paste(fitted_variables, collapse = ", "), ". Provide `var` explicitly."
      )
    }
  }

  if (!is.character(var) || !length(var) || anyNA(var) || any(var == "")) {
    stop("`var` must be a non-empty character scalar or vector.")
  }
  if (anyDuplicated(var)) {
    stop("`var` must contain unique exposure-variable names.")
  }
  invalid_variables <- setdiff(var, fitted_variables)
  if (length(invalid_variables)) {
    stop(
      "Variables not found in the fitted model: ",
      paste(invalid_variables, collapse = ", "), "."
    )
  }

  # ------------------------------------------------------------
  # Input strategy
  # ------------------------------------------------------------

  using_profile <- !is.null(profile)
  using_data <- !is.null(data)

  if (using_profile && using_data) {
    stop("Provide either `profile` or `data`, not both.")
  }
  if (!using_profile && !using_data) {
    stop("Supply either a chronological `profile` or long-format `data`.")
  }

  profiles <- NULL

  if (using_profile) {
    if (is.numeric(profile)) {
      if (length(var) != 1L) {
        stop("A numeric `profile` can only be used with one exposure variable.")
      }
      if (any(!is.finite(profile))) {
        stop("`profile` must contain only finite numeric values.")
      }
      profiles <- stats::setNames(list(as.numeric(profile)), var)
    } else if (is.list(profile)) {
      if (!length(profile)) stop("`profile` cannot be an empty list.")
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

  # ------------------------------------------------------------
  # Long-format data validation
  # ------------------------------------------------------------

  levels_to_use <- NULL

  if (using_data) {
    if (!is.data.frame(data)) stop("`data` must be a data.frame.")

    if (is.null(group) || !is.character(group) || length(group) != 1L ||
        is.na(group) || group == "") {
      stop("When using `data`, `group` must be one column name.")
    }

    required_columns <- c(group, "time", var)
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns)) {
      stop("Required columns missing from `data`: ",
           paste(missing_columns, collapse = ", "), ".")
    }

    if (anyNA(data[[group]])) {
      stop("The grouping column cannot contain missing values.")
    }

    if (!is.numeric(data$time) || any(!is.finite(data$time))) {
      stop("`time` must contain only finite numeric values.")
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
        stop("Group levels not found in `data`: ",
             paste(missing_levels, collapse = ", "), ".")
      }
      levels_to_use <- group_level
    }

    if (!length(levels_to_use)) {
      stop("No group levels are available for ECI calculation.")
    }
  }

  # ------------------------------------------------------------
  # Coefficients and link
  # ------------------------------------------------------------

  coefficient_info <- extract_coef_vcov(fit)
  beta_full <- coefficient_info$beta

  if (is.null(names(beta_full)) || anyNA(names(beta_full)) || any(names(beta_full) == "")) {
    stop("Could not determine valid names for fitted coefficients.")
  }

  linkinv <- get_linkinv(fitted_family)

  # ------------------------------------------------------------
  # Worker for one variable and one chronological profile
  # ------------------------------------------------------------

  compute_one_variable <- function(variable, profile_values) {
    variable_specification <- specification[[variable]]
    if (is.null(variable_specification)) {
      stop("Missing exposure specification for variable '", variable, "'.")
    }

    lag_max_use <- as.integer(max(variable_specification$lag_max))
    expected_length <- lag_max_use + 1L

    if (length(profile_values) != expected_length) {
      stop(
        "Profile length for variable '", variable,
        "' must equal lag_max + 1. Expected ", expected_length,
        " values but received ", length(profile_values), "."
      )
    }
    if (!is.numeric(profile_values) || any(!is.finite(profile_values))) {
      stop("Profile for variable '", variable, "' must contain finite numeric values.")
    }

    profile_values <- as.numeric(profile_values)
    eci_raw <- sum(profile_values)

    cross_basis <- dlnm::crossbasis(
      profile_values,
      lag = lag_max_use,
      argvar = variable_specification$argvar,
      arglag = variable_specification$arglag
    )

    cross_basis_row <- as.numeric(cross_basis[expected_length, , drop = TRUE])
    if (any(!is.finite(cross_basis_row))) {
      stop("The reconstructed cross-basis row contains non-finite values for '",
           variable, "'.")
    }

    coefficient_names <- get_cb_names(
      coefficient_names = names(beta_full),
      var_sel = variable,
      expected_columns = length(cross_basis_row),
      cb_cols_fit = fitted_cb_columns,
      data_template = fitted_data_template
    )

    if (!uncertainty) {
      beta_cb <- as.numeric(beta_full[coefficient_names])
      eta_weighted <- sum(cross_basis_row * beta_cb)
      eci_weighted <- transform_scale(eta_weighted, scale, linkinv)

      return(data.frame(
        var = variable,
        ECI_raw = eci_raw,
        ECI_weighted = eci_weighted,
        scale = scale,
        stringsAsFactors = FALSE
      ))
    }

    beta_draws <- extract_beta_draws(
      model = fit,
      cb_names_ref = coefficient_names,
      requested_samples = n_samples,
      coefficient_info = coefficient_info
    )

    coefficient_draw_sd <- apply(beta_draws, 2, stats::sd)
    if (all(!is.finite(coefficient_draw_sd)) ||
        all(coefficient_draw_sd < 1e-12, na.rm = TRUE)) {
      warning(
        "Near-zero coefficient-draw variability for variable '", variable,
        "'. ECI intervals may collapse.",
        call. = FALSE
      )
    }

    eta_draws <- as.numeric(beta_draws %*% cross_basis_row)
    eci_draws <- transform_scale(eta_draws, scale, linkinv)

    if (output == "samples") {
      return(data.frame(
        var = variable,
        sample = seq_along(eci_draws),
        ECI_raw = eci_raw,
        ECI_weighted = eci_draws,
        scale = scale,
        stringsAsFactors = FALSE
      ))
    }

    interval <- safe_quantile(eci_draws)
    data.frame(
      var = variable,
      ECI_raw = eci_raw,
      ECI_weighted = stats::median(eci_draws, na.rm = TRUE),
      sd = safe_sd(eci_draws),
      lower = interval[1],
      upper = interval[2],
      scale = scale,
      stringsAsFactors = FALSE
    )
  }

  # ------------------------------------------------------------
  # Direct profiles
  # ------------------------------------------------------------

  if (using_profile) {
    output_list <- lapply(var, function(variable) {
      compute_one_variable(variable, profiles[[variable]])
    })
    output_data <- do.call(rbind, output_list)
    rownames(output_data) <- NULL

    if (length(var) == 1L && is.numeric(profile)) {
      output_data$var <- NULL
    }

    return(output_data)
  }

  # ------------------------------------------------------------
  # Profiles extracted from data
  # ------------------------------------------------------------

  output_by_group <- lapply(levels_to_use, function(current_level) {
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

    group_output_list <- lapply(var, function(variable) {
      compute_one_variable(variable, as.numeric(group_data[[variable]]))
    })

    group_output <- do.call(rbind, group_output_list)
    rownames(group_output) <- NULL
    group_output[[group]] <- current_level
    group_output[, c(group, setdiff(names(group_output), group)), drop = FALSE]
  })

  output_data <- do.call(rbind, output_by_group)
  rownames(output_data) <- NULL
  output_data
}
