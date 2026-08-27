#' Predict outcomes under one or more chronological exposure profiles
#'
#' Predicts the expected outcome associated with one or more user-defined
#' exposure histories using the DLNM specification and model metadata stored by
#' `fit_epidlnm()`.
#'
#' Profiles must be ordered chronologically, from the oldest to the most recent
#' observation. The final value is associated internally with lag 0 and the
#' first value with the maximum lag.
#'
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param profiles A numeric vector for one profile in a single-exposure model,
#'   or a named list with one element per fitted exposure. Each exposure element
#'   may be: (1) one numeric vector; (2) a list of numeric vectors; (3) a matrix
#'   or data.frame with one profile per row. Every profile must contain
#'   `lag_max + 1` values. Across exposures, numbers of profiles must be equal,
#'   except that a single profile is recycled across the other exposures.
#' @param re Character. `"population"` excludes random effects;
#'   `"conditional"` includes them where supported.
#' @param id Optional grouping levels. Supply one value to recycle, or one value
#'   per exposure profile.
#' @param allow_new_levels Logical. Allow unseen grouping levels where supported.
#' @param type Character. Prediction scale: `"response"` or `"link"`.
#' @param uncertainty Logical. If `TRUE`, propagate coefficient or posterior
#'   uncertainty.
#' @param output Character. `"summary"` returns the median, standard deviation,
#'   and empirical 95 percent interval. `"samples"` returns individual draws.
#' @param n_samples Positive integer number of draws. At least two are required
#'   when `uncertainty = TRUE`.
#' @param seed Optional finite integer seed.
#'
#' @return A data.frame. When more than one exposure profile is supplied, the
#'   `profile` column identifies profile combinations in their input order.
#'   Deterministic output contains `prediction`; uncertainty summaries also
#'   contain `sd`, `lower`, and `upper`; sample output contains `sample`.
#'
#' @details
#' A list of 100 profiles for each fitted exposure produces 100 matched
#' predictions: profile 1 of every exposure is combined, then profile 2, and so
#' forth. The function does not construct the Cartesian product of profiles.
#' A single profile for one exposure may be recycled when other exposures have
#' multiple profiles.
#'
#' When `uncertainty = FALSE`, central engine estimates are used. For `brms`,
#' deterministic expected predictions are summarized by their median. When
#' `uncertainty = TRUE`, predictions are calculated draw by draw and summarized
#' using the median and empirical 2.5 and 97.5 percent quantiles.
#'
#' @export
predict_outcome <- function(
    fit,
    profiles,
    re = c("population", "conditional"),
    id = NULL,
    allow_new_levels = FALSE,
    type = c("response", "link"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    seed = NULL
) {
  re <- match.arg(re)
  type <- match.arg(type)
  output <- match.arg(output)

  if (is.null(fit)) stop("`fit` cannot be NULL.")
  if (!is.logical(allow_new_levels) || length(allow_new_levels) != 1L ||
      is.na(allow_new_levels)) {
    stop("`allow_new_levels` must be TRUE or FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L ||
      is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples <= 0 ||
      n_samples != as.integer(n_samples)) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)
  if (uncertainty && n_samples < 2L) {
    stop("`n_samples` must be at least 2 when `uncertainty = TRUE`.")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != as.integer(seed)) {
      stop("`seed` must be NULL or one finite integer.")
    }
    set.seed(as.integer(seed))
  }

  specification <- attr(fit, "epiexposure_spec")
  fitted_variables <- attr(fit, "epiexposure_vars")
  fitted_cb_columns <- attr(fit, "epiexposure_cb_cols")
  data_template <- attr(fit, "epiexposure_data_template")
  id_column <- attr(fit, "epiexposure_id_col")
  fitted_family <- attr(fit, "epiexposure_family")
  fitted_family_name <- attr(fit, "epiexposure_family_name")
  fitted_link <- attr(fit, "epiexposure_link")
  basis_objects <- attr(fit, "epiexposure_basis_objects")
  fitted_engine <- attr(fit, "epiexposure_engine")

  if (is.null(specification) || !is.list(specification)) {
    stop("Missing or invalid `epiexposure_spec` metadata in `fit`.")
  }
  if (is.null(fitted_variables) || !length(fitted_variables)) {
    if (identical(fitted_engine, "bdlnm") && is.list(basis_objects) &&
        !is.null(names(basis_objects)) && length(basis_objects)) {
      fitted_variables <- names(basis_objects)
    } else {
      stop("Missing `epiexposure_vars` metadata in `fit`.")
    }
  }
  if (is.null(data_template) || !is.data.frame(data_template) ||
      nrow(data_template) < 1L) {
    stop("Missing or invalid `epiexposure_data_template` metadata in `fit`.")
  }

  resolve_family_name <- function(x) {
    if (!is.null(fitted_family_name) && is.character(fitted_family_name) &&
        length(fitted_family_name) == 1L && nzchar(fitted_family_name)) {
      return(tolower(fitted_family_name))
    }
    raw <- NULL
    if (is.character(x) && length(x) >= 1L) raw <- x[[1]]
    if (is.list(x) && !is.null(x$family)) raw <- x$family[[1]]
    if (is.null(raw)) stop("Missing canonical family metadata in `fit`.")
    z <- tolower(gsub("[[:space:]-]+", "_", as.character(raw)))
    if (z %in% c("beta", "beta_family")) return("beta")
    if (z %in% c("binomial", "bernoulli")) return("binomial")
    if (z == "poisson") return("poisson")
    if (z == "gamma") return("gamma")
    if (z %in% c("gaussian", "normal")) return("gaussian")
    if (z %in% c("negbin", "nbinom", "nbinom1", "nbinom2",
                 "negative_binomial")) return("negative_binomial")
    if (z %in% c("ordinal", "cumulative")) return("ordinal")
    stop("Could not determine a supported canonical family from model metadata.")
  }

  family_name <- resolve_family_name(fitted_family)
  if (identical(family_name, "ordinal")) {
    stop("Ordinal predictions are not yet implemented in `predict_outcome()`.")
  }
  if (is.null(fitted_link) || !is.character(fitted_link) ||
      length(fitted_link) != 1L || is.na(fitted_link) || !nzchar(fitted_link)) {
    if (inherits(fitted_family, "family") && !is.null(fitted_family$link)) {
      fitted_link <- fitted_family$link
    } else if (is.list(fitted_family) && !is.null(fitted_family$link)) {
      fitted_link <- fitted_family$link[[1]]
    } else {
      fitted_link <- switch(
        family_name,
        beta = "logit", binomial = "logit", poisson = "log",
        gamma = "log", gaussian = "identity", negative_binomial = "log",
        stop("Could not determine the model link from metadata.")
      )
    }
  }
  fitted_link <- tolower(as.character(fitted_link[[1]]))

  get_linkinv <- function(link_name) {
    switch(
      tolower(link_name), identity = identity, log = exp,
      logit = stats::plogis, probit = stats::pnorm,
      cloglog = function(x) 1 - exp(-exp(x)), inverse = function(x) 1 / x,
      stop("Unsupported inverse-link function for link '", link_name, "'.")
    )
  }
  linkinv <- get_linkinv(fitted_link)
  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    x <- x[is.finite(x)]
    if (!length(x)) return(rep(NA_real_, length(probs)))
    stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)
  }
  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(0)
    stats::sd(x)
  }
  sort_cb_names <- function(x) {
    if (!length(x)) return(x)
    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    missing <- is.na(idx)
    idx[missing] <- seq_along(x)[missing]
    x[order(idx)]
  }
  apply_linkinv_matrix <- function(x) {
    d <- dim(x)
    matrix(linkinv(as.vector(x)), nrow = d[1], ncol = d[2])
  }
  ensure_intercept <- function(x, coefficient_names) {
    if ("(Intercept)" %in% coefficient_names &&
        !"(Intercept)" %in% names(x)) x[["(Intercept)"]] <- 1
    x
  }
  check_vcov_names <- function(vcov_matrix, coefficient_names) {
    if (is.null(rownames(vcov_matrix)) || is.null(colnames(vcov_matrix))) {
      stop("The coefficient covariance matrix must contain row and column names.")
    }
    missing <- setdiff(coefficient_names,
                       intersect(rownames(vcov_matrix), colnames(vcov_matrix)))
    if (length(missing)) {
      stop("Coefficients absent from the covariance matrix: ",
           paste(missing, collapse = ", "), ".")
    }
    invisible(TRUE)
  }
  is_gamm_object <- function(model) {
    is.list(model) && !is.null(model$gam) && inherits(model$gam, "gam")
  }

  # Standardize each exposure as a list of numeric profiles. This is the main
  # extension that allows one or many profiles without breaking single-profile
  # calls.
  as_profile_list <- function(x, variable) {
    if (is.numeric(x) && is.null(dim(x))) return(list(as.numeric(x)))
    if (is.matrix(x) || is.data.frame(x)) {
      x <- as.matrix(x)
      if (!is.numeric(x)) {
        stop("Matrix/data.frame profiles for '", variable, "' must be numeric.")
      }
      return(lapply(seq_len(nrow(x)), function(i) as.numeric(x[i, ])))
    }
    if (is.list(x)) {
      if (!length(x)) stop("No profiles were supplied for '", variable, "'.")
      return(lapply(seq_along(x), function(i) {
        value <- x[[i]]
        if (is.data.frame(value) || is.matrix(value)) {
          if (nrow(as.matrix(value)) != 1L) {
            stop("Nested profile ", i, " for '", variable,
                 "' must contain exactly one row.")
          }
          value <- as.numeric(as.matrix(value)[1, ])
        }
        if (!is.numeric(value) || is.list(value)) {
          stop("Profile ", i, " for '", variable,
               "' must be a numeric vector.")
        }
        as.numeric(value)
      }))
    }
    stop("Profiles for '", variable,
         "' must be a numeric vector, list, matrix, or data.frame.")
  }

  if (is.numeric(profiles) && is.null(dim(profiles))) {
    if (length(fitted_variables) != 1L) {
      stop("A numeric `profiles` vector can only be used with one fitted exposure.")
    }
    profiles <- stats::setNames(list(profiles), fitted_variables)
  } else if (!is.list(profiles)) {
    stop("`profiles` must be a numeric vector or a named list by exposure.")
  }
  if (!length(profiles) || is.null(names(profiles)) || anyNA(names(profiles)) ||
      any(names(profiles) == "") || anyDuplicated(names(profiles))) {
    stop("`profiles` must have unique non-empty exposure-variable names.")
  }
  extra_variables <- setdiff(names(profiles), fitted_variables)
  missing_variables <- setdiff(fitted_variables, names(profiles))
  if (length(extra_variables)) {
    stop("Profile variables not present in the fitted model: ",
         paste(extra_variables, collapse = ", "), ".")
  }
  if (length(missing_variables)) {
    stop("Profiles are missing fitted exposure variables: ",
         paste(missing_variables, collapse = ", "), ".")
  }
  missing_specification <- setdiff(names(profiles), names(specification))
  if (length(missing_specification)) {
    stop("Exposure specifications are missing for: ",
         paste(missing_specification, collapse = ", "), ".")
  }

  profile_sets <- lapply(names(profiles), function(variable) {
    as_profile_list(profiles[[variable]], variable)
  })
  names(profile_sets) <- names(profiles)
  profile_counts <- vapply(profile_sets, length, integer(1))
  n_profiles <- max(profile_counts)
  incompatible <- profile_counts != 1L & profile_counts != n_profiles
  if (any(incompatible)) {
    stop(
      "Exposure variables must provide either one profile or the common maximum ",
      "number of profiles. Received: ",
      paste(paste0(names(profile_counts), "=", profile_counts), collapse = ", "), "."
    )
  }
  for (variable in names(profile_sets)) {
    if (length(profile_sets[[variable]]) == 1L && n_profiles > 1L) {
      profile_sets[[variable]] <- rep(profile_sets[[variable]], n_profiles)
    }
  }

  build_cb_row <- function(profile_values, variable_specification, variable,
                           profile_index) {
    maximum_lag <- as.integer(max(variable_specification$lag_max))
    expected_length <- maximum_lag + 1L
    if (!is.numeric(profile_values) || any(!is.finite(profile_values))) {
      stop("Profile ", profile_index, " for variable '", variable,
           "' must contain finite numeric values.")
    }
    if (length(profile_values) != expected_length) {
      stop("Profile ", profile_index, " for variable '", variable,
           "' must have lag_max + 1 values. Expected ", expected_length,
           " but received ", length(profile_values), ".")
    }
    cb <- dlnm::crossbasis(
      profile_values, lag = maximum_lag,
      argvar = variable_specification$argvar,
      arglag = variable_specification$arglag
    )
    row <- as.numeric(cb[expected_length, , drop = TRUE])
    if (any(!is.finite(row))) {
      stop("The reconstructed cross-basis row is non-finite for variable '",
           variable, "', profile ", profile_index, ".")
    }
    row
  }

  profile_rows <- lapply(names(profile_sets), function(variable) {
    do.call(rbind, lapply(seq_len(n_profiles), function(i) {
      build_cb_row(profile_sets[[variable]][[i]], specification[[variable]],
                   variable, i)
    }))
  })
  names(profile_rows) <- names(profile_sets)

  derive_cb_columns <- function(variable, expected_columns,
                                coefficient_names = NULL) {
    candidates <- character(0)
    if (!is.null(fitted_cb_columns) && length(fitted_cb_columns)) {
      candidates <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"), fitted_cb_columns, value = TRUE
      ))
    }
    if (!length(candidates)) {
      candidates <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"), names(data_template), value = TRUE
      ))
    }
    if (!length(candidates) && !is.null(coefficient_names)) {
      candidates <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"), coefficient_names, value = TRUE
      ))
    }
    if (!length(candidates) && identical(fitted_engine, "bdlnm")) {
      stop("Could not map the stored bdlnm basis for variable '", variable,
           "' to named posterior coefficients.")
    }
    if (length(candidates) != expected_columns) {
      stop("Could not align cross-basis columns for variable '", variable,
           "'. Expected ", expected_columns, " columns but found ",
           length(candidates), ".")
    }
    candidates
  }

  bdlnm_coefficient_names <- NULL
  if (inherits(fit, "bdlnm")) {
    if (!is.null(fit$coefficients) && !is.null(rownames(fit$coefficients))) {
      bdlnm_coefficient_names <- rownames(fit$coefficients)
    } else if (!is.null(fit$coefficients.summary)) {
      bdlnm_coefficient_names <- rownames(fit$coefficients.summary)
    }
  }
  cb_map <- lapply(names(profile_rows), function(variable) {
    derive_cb_columns(variable, ncol(profile_rows[[variable]]),
                      bdlnm_coefficient_names)
  })
  names(cb_map) <- names(profile_rows)
  fitted_cb_columns <- unname(unlist(cb_map, use.names = FALSE))

  newdata <- data_template[rep(1L, n_profiles), , drop = FALSE]
  rownames(newdata) <- NULL
  if (!is.null(id)) {
    if (is.null(id_column) || !is.character(id_column) || !nzchar(id_column)) {
      stop("The fitted model does not contain grouping-column metadata.")
    }
    if (!id_column %in% names(newdata)) {
      stop("Grouping column '", id_column,
           "' is absent from `epiexposure_data_template`.")
    }
    if (!length(id) || anyNA(id)) {
      stop("`id` must contain at least one non-missing grouping level.")
    }
    if (length(id) == 1L) id <- rep(id, n_profiles)
    if (length(id) != n_profiles) {
      stop("`id` must have length 1 or equal the number of profiles (",
           n_profiles, ").")
    }
    newdata[[id_column]] <- id
  }
  for (variable in names(profile_rows)) {
    columns <- cb_map[[variable]]
    replacement <- profile_rows[[variable]]
    colnames(replacement) <- columns
    newdata[columns] <- as.data.frame(replacement)
  }

  expected_fixed_terms <- function(coefficient_names, input_data) {
    terms <- unique(c(if ("(Intercept)" %in% coefficient_names) "(Intercept)",
                      fitted_cb_columns))
    missing_coefficients <- setdiff(terms, coefficient_names)
    missing_data <- setdiff(terms, names(input_data))
    if (length(missing_coefficients)) {
      stop("Expected coefficients missing from the model: ",
           paste(missing_coefficients, collapse = ", "), ".")
    }
    if (length(missing_data)) {
      stop("Expected fixed-effect columns missing from newdata: ",
           paste(missing_data, collapse = ", "), ".")
    }
    terms
  }
  manual_fixed_prediction <- function(beta, input_data, requested_type) {
    if (is.null(names(beta))) stop("The fixed-effect coefficient vector has no names.")
    input_data <- ensure_intercept(input_data, names(beta))
    terms <- expected_fixed_terms(names(beta), input_data)
    eta <- as.numeric(as.matrix(input_data[, terms, drop = FALSE]) %*%
                        as.numeric(beta[terms]))
    if (requested_type == "link") return(eta)
    linkinv(eta)
  }
  point_predict <- function(model, input_data, requested_re, requested_type) {
    population <- identical(requested_re, "population")
    if (inherits(model, "glmmTMB")) return(as.numeric(stats::predict(
      model, newdata = input_data, type = requested_type,
      re.form = if (population) NA else NULL,
      allow.new.levels = allow_new_levels
    )))
    if (inherits(model, "merMod")) return(as.numeric(stats::predict(
      model, newdata = input_data,
      type = if (requested_type == "link") "link" else "response",
      re.form = if (population) NA else NULL,
      allow.new.levels = allow_new_levels
    )))
    if (inherits(model, "brmsfit")) {
      values <- if (requested_type == "link") {
        brms::posterior_linpred(model, newdata = input_data,
                                re_formula = if (population) NA else NULL,
                                allow_new_levels = allow_new_levels)
      } else {
        brms::posterior_epred(model, newdata = input_data,
                              re_formula = if (population) NA else NULL,
                              allow_new_levels = allow_new_levels)
      }
      return(apply(as.matrix(values), 2, stats::median, na.rm = TRUE))
    }
    if (is_gamm_object(model)) return(as.numeric(stats::predict(
      model$gam, newdata = input_data,
      type = if (requested_type == "link") "link" else "response")))
    if (inherits(model, "gam")) return(as.numeric(stats::predict(
      model, newdata = input_data,
      type = if (requested_type == "link") "link" else "response")))
    if (inherits(model, "lme")) return(as.numeric(nlme::predict.lme(
      model, newdata = input_data, level = if (population) 0 else 1)))
    if (inherits(model, "gls")) {
      return(as.numeric(nlme::predict.gls(model, newdata = input_data)))
    }
    if (inherits(model, "HLfit")) {
      direct <- tryCatch(as.numeric(stats::predict(
        model, newdata = input_data,
        type = if (requested_type == "link") "link" else "response",
        re.form = if (population) NA else NULL)), error = function(e) NULL)
      if (!is.null(direct) && all(is.finite(direct))) return(direct)
      return(manual_fixed_prediction(spaMM::fixef(model), input_data,
                                     requested_type))
    }
    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      if (is.null(names(beta))) names(beta) <- rownames(model$summary.fixed)
      if (!population) warning(
        "INLA conditional prediction is not implemented; using fixed effects.",
        call. = FALSE)
      return(manual_fixed_prediction(beta, input_data, requested_type))
    }
    if (inherits(model, "bdlnm")) {
      beta <- model$coefficients.summary[, "mean"]
      if (is.null(names(beta))) names(beta) <- rownames(model$coefficients.summary)
      return(manual_fixed_prediction(beta, input_data, requested_type))
    }
    direct <- tryCatch(as.numeric(stats::predict(
      model, newdata = input_data,
      type = if (requested_type == "link") "link" else "response")),
      error = function(e) NULL)
    if (!is.null(direct) && all(is.finite(direct))) return(direct)
    manual_fixed_prediction(stats::coef(model), input_data, requested_type)
  }

  add_profile_columns <- function(out, profile_values) {
    if (n_profiles > 1L) out$profile <- profile_values
    if (!is.null(id)) out[[id_column]] <- id[profile_values]
    leading <- c(if (!is.null(id)) id_column,
                 if (n_profiles > 1L) "profile")
    out[, c(leading, setdiff(names(out), leading)), drop = FALSE]
  }
  format_samples <- function(draw_matrix) {
    draw_matrix <- as.matrix(draw_matrix)
    rows <- lapply(seq_len(ncol(draw_matrix)), function(j) {
      x <- data.frame(sample = seq_len(nrow(draw_matrix)),
                      prediction = as.numeric(draw_matrix[, j]),
                      stringsAsFactors = FALSE)
      add_profile_columns(x, rep(j, nrow(draw_matrix)))
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
  format_summary <- function(draw_matrix) {
    draw_matrix <- as.matrix(draw_matrix)
    q <- t(apply(draw_matrix, 2, safe_quantile))
    out <- data.frame(
      prediction = apply(draw_matrix, 2, stats::median, na.rm = TRUE),
      sd = apply(draw_matrix, 2, safe_sd), lower = q[, 1], upper = q[, 2],
      stringsAsFactors = FALSE
    )
    add_profile_columns(out, seq_len(n_profiles))
  }

  if (!uncertainty) {
    out <- data.frame(prediction = point_predict(fit, newdata, re, type),
                      stringsAsFactors = FALSE)
    return(add_profile_columns(out, seq_len(n_profiles)))
  }
  if (inherits(fit, "brmsfit")) {
    draws <- if (type == "link") {
      brms::posterior_linpred(fit, newdata = newdata,
                              re_formula = if (re == "population") NA else NULL,
                              allow_new_levels = allow_new_levels, ndraws = n_samples)
    } else {
      brms::posterior_epred(fit, newdata = newdata,
                            re_formula = if (re == "population") NA else NULL,
                            allow_new_levels = allow_new_levels, ndraws = n_samples)
    }
    draws <- as.matrix(draws)
    if (output == "samples") return(format_samples(draws))
    return(format_summary(draws))
  }
  if (inherits(fit, "inla")) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for INLA uncertainty.")
    }
    samples <- tryCatch(INLA::inla.posterior.sample(n = n_samples, result = fit),
                        error = function(e) NULL)
    if (is.null(samples)) stop(
      "INLA posterior samples could not be drawn. Fit with ",
      "`control.compute = list(config = TRUE)`."
    )
    beta_names <- rownames(fit$summary.fixed)
    if (is.null(beta_names)) beta_names <- names(fit$summary.fixed$mean)
    inla_data <- ensure_intercept(newdata, beta_names)
    terms <- expected_fixed_terms(beta_names, inla_data)
    beta_draws <- do.call(rbind, lapply(samples, function(x) {
      latent <- x$latent
      names(latent) <- gsub(":1$", "", names(latent))
      values <- latent[terms]
      if (anyNA(values)) stop("Could not match INLA draws to fixed-effect terms.")
      as.numeric(values)
    }))
    draws <- beta_draws %*% t(as.matrix(inla_data[, terms, drop = FALSE]))
    if (type != "link") draws <- apply_linkinv_matrix(draws)
    if (re != "population") warning(
      "INLA conditional uncertainty is not implemented; using fixed effects.",
      call. = FALSE)
    if (output == "samples") return(format_samples(draws))
    return(format_summary(draws))
  }
  if (inherits(fit, "bdlnm")) {
    coefficient_draws <- fit$coefficients
    if (is.null(dim(coefficient_draws))) {
      coefficient_draws <- matrix(coefficient_draws, ncol = 1L)
    }
    if (is.null(rownames(coefficient_draws))) {
      stop("The bdlnm coefficient-draw matrix has no row names.")
    }
    bdlnm_data <- ensure_intercept(newdata, rownames(coefficient_draws))
    terms <- expected_fixed_terms(rownames(coefficient_draws), bdlnm_data)
    coefficient_draws <- coefficient_draws[terms, , drop = FALSE]
    if (ncol(coefficient_draws) > n_samples) {
      keep <- sample(seq_len(ncol(coefficient_draws)), n_samples)
      coefficient_draws <- coefficient_draws[, keep, drop = FALSE]
    }
    draws <- t(as.matrix(bdlnm_data[, terms, drop = FALSE]) %*%
                 coefficient_draws)
    if (type != "link") draws <- apply_linkinv_matrix(draws)
    if (output == "samples") return(format_samples(draws))
    return(format_summary(draws))
  }

  extract_coef_vcov <- function(model) {
    if (inherits(model, "glmmTMB")) return(list(
      beta = glmmTMB::fixef(model)$cond,
      vcov = as.matrix(stats::vcov(model)$cond)))
    if (inherits(model, "merMod")) return(list(
      beta = lme4::fixef(model), vcov = as.matrix(stats::vcov(model))))
    if (is_gamm_object(model)) return(list(
      beta = stats::coef(model$gam), vcov = as.matrix(stats::vcov(model$gam))))
    if (inherits(model, "lme")) return(list(
      beta = nlme::fixef(model), vcov = as.matrix(stats::vcov(model))))
    if (inherits(model, "gls") || inherits(model, "gam")) return(list(
      beta = stats::coef(model), vcov = as.matrix(stats::vcov(model))))
    if (inherits(model, "HLfit")) return(list(
      beta = spaMM::fixef(model), vcov = as.matrix(stats::vcov(model))))
    beta <- stats::coef(model)
    V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
    if (!is.numeric(beta) || is.null(V)) {
      stop("Could not extract coefficients and covariance matrix from the model.")
    }
    list(beta = beta, vcov = V)
  }
  info <- extract_coef_vcov(fit)
  beta_hat <- info$beta
  covariance_hat <- info$vcov
  fixed_data <- ensure_intercept(newdata, names(beta_hat))
  fixed_terms <- expected_fixed_terms(names(beta_hat), fixed_data)
  check_vcov_names(covariance_hat, fixed_terms)
  beta_hat <- beta_hat[fixed_terms]
  covariance_hat <- covariance_hat[fixed_terms, fixed_terms, drop = FALSE]
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for frequentist uncertainty.")
  }
  beta_draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat,
                              Sigma = covariance_hat)
  if (is.null(dim(beta_draws))) beta_draws <- matrix(beta_draws, nrow = 1L)
  draws <- beta_draws %*% t(as.matrix(fixed_data[, fixed_terms, drop = FALSE]))
  if (re != "population" &&
      (inherits(fit, "glmmTMB") || inherits(fit, "merMod") ||
       inherits(fit, "lme"))) {
    warning("Frequentist mixed-model uncertainty reflects fixed effects only.",
            call. = FALSE)
  }
  if (type != "link") draws <- apply_linkinv_matrix(draws)
  if (output == "samples") return(format_samples(draws))
  format_summary(draws)
}
