#' Predict outcomes under user-defined chronological exposure profiles
#'
#' Predicts the outcome associated with one or more user-defined exposure
#' profiles using the DLNM specification and model metadata stored by
#' `fit_epidlnm()`.
#'
#' Every profile must be supplied in chronological order, from the earliest
#' observation to the most recent observation. Profiles are not reversed
#' internally. When the final cross-basis row is reconstructed, the last profile
#' value is associated with lag 0 and the first value with lag max.
#'
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param profiles Numeric vector for a single-exposure model or a named list of
#'   numeric vectors for one or more exposures. Each profile must be chronological
#'   and have length `lag_max + 1`, using the lag stored for that exposure in the
#'   fitted-model specification.
#' @param re Character. Prediction level: `"population"` excludes random
#'   effects; `"conditional"` includes them where supported.
#' @param id Optional vector of grouping levels. When supplied, one prediction is
#'   generated for each value. The fitted model must contain the grouping-column
#'   metadata stored by `fit_epidlnm()`.
#' @param allow_new_levels Logical. Allow unseen grouping levels where supported.
#' @param type Character. `"response"`, `"link"`, or `"conditional"`.
#'   `"conditional"` requests a conditional prediction on the response scale.
#' @param uncertainty Logical. If `TRUE`, propagate coefficient uncertainty.
#' @param output Character. `"summary"` returns median, standard deviation, and
#'   empirical 95 percent intervals. `"samples"` returns individual draws.
#' @param n_samples Positive integer number of coefficient or posterior draws.
#'
#' @return A data.frame. Deterministic output contains `prediction`. Summary
#'   output additionally contains `sd`, `lower`, and `upper`. Sample output
#'   contains `sample` and `prediction`. The grouping column is included when
#'   `id` is supplied.
#'
#' @details
#' Bayesian engines use posterior draws when available. Frequentist engines use
#' draws from the asymptotic multivariate normal distribution of the fixed-effect
#' coefficients. For frequentist mixed models, uncertainty reflects fixed
#' effects only.
#'
#' Cross-basis columns are aligned using `epiexposure_cb_cols`, with
#' `epiexposure_data_template` as a fallback. Coefficient and covariance-matrix
#' names are checked before simulation to prevent silent term omission.
#'
#' @export
predict_outcome <- function(
    fit,
    profiles,
    re = c("population", "conditional"),
    id = NULL,
    allow_new_levels = FALSE,
    type = c("response", "link", "conditional"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
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
      !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  specification <- attr(fit, "epiexposure_spec")
  fitted_variables <- attr(fit, "epiexposure_vars")
  fitted_cb_columns <- attr(fit, "epiexposure_cb_cols")
  data_template <- attr(fit, "epiexposure_data_template")
  id_column <- attr(fit, "epiexposure_id_col")
  fitted_family <- attr(fit, "epiexposure_family")

  if (is.null(specification) || !is.list(specification)) {
    stop("Missing or invalid `epiexposure_spec` metadata in `fit`.")
  }
  if (is.null(fitted_variables) || !length(fitted_variables)) {
    stop("Missing `epiexposure_vars` metadata in `fit`.")
  }
  if (is.null(fitted_cb_columns) || !length(fitted_cb_columns)) {
    stop("Missing `epiexposure_cb_cols` metadata in `fit`.")
  }
  if (is.null(data_template) || !is.data.frame(data_template) ||
      nrow(data_template) < 1L) {
    stop("Missing or invalid `epiexposure_data_template` metadata in `fit`.")
  }
  if (is.null(fitted_family)) {
    stop("Missing `epiexposure_family` metadata in `fit`.")
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

  get_cb_names <- function(variable, expected_columns) {
    matched <- sort_cb_names(grep(
      paste0("^cb_", variable, "_"),
      fitted_cb_columns,
      value = TRUE
    ))

    if (!length(matched)) {
      matched <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"),
        names(data_template),
        value = TRUE
      ))
    }

    if (length(matched) != expected_columns) {
      stop(
        "Could not align cross-basis columns for variable '", variable,
        "'. Expected ", expected_columns, " columns but found ",
        length(matched), "."
      )
    }

    missing_template <- setdiff(matched, names(data_template))
    if (length(missing_template)) {
      stop(
        "Cross-basis columns missing from `epiexposure_data_template`: ",
        paste(missing_template, collapse = ", "), "."
      )
    }

    matched
  }

  build_cb_row <- function(profile_values, variable_specification, variable) {
    maximum_lag <- as.integer(max(variable_specification$lag_max))
    expected_length <- maximum_lag + 1L

    if (!is.numeric(profile_values) || any(!is.finite(profile_values))) {
      stop("Profile for variable '", variable, "' must contain finite numeric values.")
    }
    if (length(profile_values) != expected_length) {
      stop(
        "Profile length for variable '", variable,
        "' must equal lag_max + 1. Expected ", expected_length,
        " values but received ", length(profile_values), "."
      )
    }

    cross_basis <- dlnm::crossbasis(
      as.numeric(profile_values),
      lag = maximum_lag,
      argvar = variable_specification$argvar,
      arglag = variable_specification$arglag
    )

    cross_basis_row <- as.numeric(
      cross_basis[expected_length, , drop = TRUE]
    )

    if (any(!is.finite(cross_basis_row))) {
      stop(
        "The reconstructed cross-basis row contains non-finite values for variable '",
        variable, "'."
      )
    }

    cross_basis_row
  }

  get_linkinv <- function(family_object) {
    if (inherits(family_object, "family") && !is.null(family_object$linkinv)) {
      return(family_object$linkinv)
    }
    if (is.character(family_object)) {
      family_name <- tolower(family_object[1])
      if (family_name %in% c("beta", "binomial")) return(stats::plogis)
      if (family_name %in% c(
        "poisson", "gamma", "negbin", "negative_binomial"
      )) return(exp)
      if (family_name == "gaussian") return(identity)
    }
    stop("Could not determine the inverse-link function from model metadata.")
  }

  apply_linkinv_matrix <- function(matrix_object, linkinv) {
    dimensions <- dim(matrix_object)
    transformed <- linkinv(as.vector(matrix_object))
    matrix(transformed, nrow = dimensions[1], ncol = dimensions[2])
  }

  ensure_intercept <- function(input_data, coefficient_names) {
    if ("(Intercept)" %in% coefficient_names &&
        !"(Intercept)" %in% names(input_data)) {
      input_data[["(Intercept)"]] <- 1
    }
    input_data
  }

  check_vcov_names <- function(covariance_matrix, coefficient_names) {
    if (is.null(rownames(covariance_matrix)) ||
        is.null(colnames(covariance_matrix))) {
      stop("The coefficient covariance matrix must contain row and column names.")
    }
    missing_names <- setdiff(
      coefficient_names,
      intersect(rownames(covariance_matrix), colnames(covariance_matrix))
    )
    if (length(missing_names)) {
      stop(
        "Coefficients absent from the covariance matrix: ",
        paste(missing_names, collapse = ", "), "."
      )
    }
    invisible(TRUE)
  }

  is_gamm_object <- function(model) {
    is.list(model) && !is.null(model$gam) && inherits(model$gam, "gam")
  }

  # Normalize and validate the chronological profiles.
  if (is.numeric(profiles)) {
    if (length(fitted_variables) != 1L) {
      stop("A numeric `profiles` vector can only be used with one fitted exposure.")
    }
    profiles <- stats::setNames(list(as.numeric(profiles)), fitted_variables)
  } else if (!is.list(profiles)) {
    stop("`profiles` must be a numeric vector or a named list of numeric vectors.")
  }

  if (is.null(names(profiles)) || anyNA(names(profiles)) ||
      any(names(profiles) == "") || anyDuplicated(names(profiles))) {
    stop("`profiles` must have unique non-empty variable names.")
  }

  profile_variables <- names(profiles)
  extra_variables <- setdiff(profile_variables, fitted_variables)
  missing_variables <- setdiff(fitted_variables, profile_variables)

  if (length(extra_variables)) {
    stop("Profile variables not present in the fitted model: ",
         paste(extra_variables, collapse = ", "), ".")
  }
  if (length(missing_variables)) {
    stop("Profiles are missing fitted exposure variables: ",
         paste(missing_variables, collapse = ", "), ".")
  }

  missing_specification <- setdiff(profile_variables, names(specification))
  if (length(missing_specification)) {
    stop("Exposure specifications are missing for: ",
         paste(missing_specification, collapse = ", "), ".")
  }

  for (variable in profile_variables) {
    if (!is.numeric(profiles[[variable]]) ||
        any(!is.finite(profiles[[variable]]))) {
      stop("Profile for variable '", variable, "' must contain finite numeric values.")
    }
    profiles[[variable]] <- as.numeric(profiles[[variable]])
  }

  # Build newdata from the fitted-model template.
  newdata <- data_template[1, , drop = FALSE]

  if (!is.null(id)) {
    if (is.null(id_column) || !nzchar(id_column)) {
      stop("The fitted model does not contain grouping-column metadata.")
    }
    if (!id_column %in% names(newdata)) {
      stop("Grouping column '", id_column,
           "' is absent from `epiexposure_data_template`.")
    }
    if (!length(id) || anyNA(id)) {
      stop("`id` must contain at least one non-missing grouping level.")
    }
    newdata <- newdata[rep(1L, length(id)), , drop = FALSE]
    newdata[[id_column]] <- id
  }

  for (variable in profile_variables) {
    cross_basis_row <- build_cb_row(
      profile_values = profiles[[variable]],
      variable_specification = specification[[variable]],
      variable = variable
    )
    cross_basis_columns <- get_cb_names(
      variable = variable,
      expected_columns = length(cross_basis_row)
    )

    replacement <- matrix(
      rep(cross_basis_row, each = nrow(newdata)),
      nrow = nrow(newdata),
      ncol = length(cross_basis_row),
      byrow = FALSE
    )
    colnames(replacement) <- cross_basis_columns
    newdata[cross_basis_columns] <- as.data.frame(replacement)
  }

  expected_fixed_terms <- function(coefficient_names, input_data) {
    terms <- unique(c(
      if ("(Intercept)" %in% coefficient_names) "(Intercept)",
      fitted_cb_columns
    ))
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

  manual_fixed_prediction <- function(
    coefficient_vector,
    input_data,
    requested_type
  ) {
    if (is.null(names(coefficient_vector))) {
      stop("The fixed-effect coefficient vector has no names.")
    }
    input_data <- ensure_intercept(input_data, names(coefficient_vector))
    terms <- expected_fixed_terms(names(coefficient_vector), input_data)
    eta <- as.numeric(
      as.matrix(input_data[, terms, drop = FALSE]) %*%
        as.numeric(coefficient_vector[terms])
    )
    if (requested_type == "link") return(eta)
    get_linkinv(fitted_family)(eta)
  }

  point_predict <- function(model, input_data, requested_re, requested_type) {
    population <- identical(requested_re, "population")
    type_use <- if (requested_type == "conditional") "response" else requested_type

    if (inherits(model, "glmmTMB")) {
      return(as.numeric(stats::predict(
        model,
        newdata = input_data,
        type = type_use,
        re.form = if (population) NA else NULL,
        allow.new.levels = allow_new_levels
      )))
    }

    if (inherits(model, "merMod")) {
      return(as.numeric(stats::predict(
        model,
        newdata = input_data,
        type = if (type_use == "link") "link" else "response",
        re.form = if (population) NA else NULL,
        allow.new.levels = allow_new_levels
      )))
    }

    if (inherits(model, "brmsfit")) {
      re_formula <- if (population) NA else NULL
      posterior_values <- if (type_use == "link") {
        brms::posterior_linpred(
          model,
          newdata = input_data,
          re_formula = re_formula,
          allow_new_levels = allow_new_levels
        )
      } else {
        brms::posterior_epred(
          model,
          newdata = input_data,
          re_formula = re_formula,
          allow_new_levels = allow_new_levels
        )
      }
      return(apply(posterior_values, 2, stats::median, na.rm = TRUE))
    }

    if (is_gamm_object(model)) {
      return(as.numeric(stats::predict(
        model$gam,
        newdata = input_data,
        type = if (type_use == "link") "link" else "response"
      )))
    }

    if (inherits(model, "gam")) {
      return(as.numeric(stats::predict(
        model,
        newdata = input_data,
        type = if (type_use == "link") "link" else "response"
      )))
    }

    if (inherits(model, "lme")) {
      prediction <- as.numeric(nlme::predict.lme(
        model,
        newdata = input_data,
        level = if (population) 0 else 1
      ))
      if (type_use == "link") return(prediction)
      return(prediction)
    }

    if (inherits(model, "gls")) {
      return(as.numeric(nlme::predict.gls(model, newdata = input_data)))
    }

    if (inherits(model, "HLfit")) {
      direct <- tryCatch(
        as.numeric(stats::predict(
          model,
          newdata = input_data,
          type = if (type_use == "link") "link" else "response",
          re.form = if (population) NA else NULL
        )),
        error = function(e) NULL
      )
      if (!is.null(direct) && all(is.finite(direct))) return(direct)
      return(manual_fixed_prediction(spaMM::fixef(model), input_data, type_use))
    }

    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      if (is.null(names(beta))) names(beta) <- rownames(model$summary.fixed)
      if (!population) {
        warning("INLA conditional prediction is not implemented; using fixed effects.",
                call. = FALSE)
      }
      return(manual_fixed_prediction(beta, input_data, type_use))
    }

    if (inherits(model, "bdlnm")) {
      if (is.null(model$coefficients.summary)) {
        stop("The bdlnm object does not contain coefficient summaries.")
      }
      beta <- model$coefficients.summary[, "mean"]
      return(manual_fixed_prediction(beta, input_data, type_use))
    }

    direct <- tryCatch(
      as.numeric(stats::predict(
        model,
        newdata = input_data,
        type = if (type_use == "link") "link" else "response"
      )),
      error = function(e) NULL
    )
    if (!is.null(direct) && all(is.finite(direct))) return(direct)

    beta <- stats::coef(model)
    manual_fixed_prediction(beta, input_data, type_use)
  }

  format_samples <- function(draw_matrix) {
    draw_matrix <- as.matrix(draw_matrix)
    output_list <- lapply(seq_len(ncol(draw_matrix)), function(column_index) {
      current <- data.frame(
        sample = seq_len(nrow(draw_matrix)),
        prediction = as.numeric(draw_matrix[, column_index]),
        stringsAsFactors = FALSE
      )
      if (!is.null(id)) current[[id_column]] <- id[column_index]
      current
    })
    output_data <- do.call(rbind, output_list)
    rownames(output_data) <- NULL
    if (!is.null(id)) {
      output_data <- output_data[, c(id_column, "sample", "prediction"), drop = FALSE]
    }
    output_data
  }

  format_summary <- function(draw_matrix) {
    draw_matrix <- as.matrix(draw_matrix)
    quantiles <- t(apply(draw_matrix, 2, safe_quantile))
    output_data <- data.frame(
      prediction = apply(draw_matrix, 2, stats::median, na.rm = TRUE),
      sd = apply(draw_matrix, 2, safe_sd),
      lower = quantiles[, 1],
      upper = quantiles[, 2],
      stringsAsFactors = FALSE
    )
    if (!is.null(id)) {
      output_data[[id_column]] <- id
      output_data <- output_data[
        ,
        c(id_column, "prediction", "sd", "lower", "upper"),
        drop = FALSE
      ]
    }
    output_data
  }

  if (!uncertainty) {
    prediction <- point_predict(fit, newdata, re, type)
    output_data <- data.frame(
      prediction = as.numeric(prediction),
      stringsAsFactors = FALSE
    )
    if (!is.null(id)) {
      output_data[[id_column]] <- id
      output_data <- output_data[, c(id_column, "prediction"), drop = FALSE]
    }
    return(output_data)
  }

  type_use <- if (type == "conditional") "response" else type

  if (inherits(fit, "brmsfit")) {
    posterior_draws <- if (type_use == "link") {
      brms::posterior_linpred(
        fit,
        newdata = newdata,
        re_formula = if (re == "population") NA else NULL,
        allow_new_levels = allow_new_levels,
        ndraws = n_samples
      )
    } else {
      brms::posterior_epred(
        fit,
        newdata = newdata,
        re_formula = if (re == "population") NA else NULL,
        allow_new_levels = allow_new_levels,
        ndraws = n_samples
      )
    }
    posterior_draws <- as.matrix(posterior_draws)
    if (output == "samples") return(format_samples(posterior_draws))
    return(format_summary(posterior_draws))
  }

  if (inherits(fit, "inla")) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for INLA uncertainty.")
    }
    posterior_samples <- tryCatch(
      INLA::inla.posterior.sample(n = n_samples, result = fit),
      error = function(e) NULL
    )
    if (is.null(posterior_samples)) {
      stop("INLA posterior samples could not be drawn; fit with config = TRUE.")
    }

    beta_names <- rownames(fit$summary.fixed)
    if (is.null(beta_names)) beta_names <- names(fit$summary.fixed$mean)
    if (is.null(beta_names)) stop("Could not determine INLA fixed-effect names.")
    inla_data <- ensure_intercept(newdata, beta_names)
    terms <- expected_fixed_terms(beta_names, inla_data)

    beta_draws <- do.call(rbind, lapply(posterior_samples, function(sample_object) {
      latent <- sample_object$latent
      names(latent) <- gsub(":1$", "", names(latent))
      values <- latent[terms]
      if (anyNA(values)) stop("Could not match INLA draws to fixed-effect terms.")
      as.numeric(values)
    }))
    colnames(beta_draws) <- terms
    eta_draws <- beta_draws %*% t(as.matrix(inla_data[, terms, drop = FALSE]))
    if (type_use != "link") {
      eta_draws <- apply_linkinv_matrix(eta_draws, get_linkinv(fitted_family))
    }
    if (re != "population") {
      warning("INLA conditional uncertainty is not implemented; using fixed effects.",
              call. = FALSE)
    }
    if (output == "samples") return(format_samples(eta_draws))
    return(format_summary(eta_draws))
  }

  if (inherits(fit, "bdlnm")) {
    if (is.null(fit$coefficients)) {
      stop("The bdlnm object does not contain posterior coefficient draws.")
    }
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
      coefficient_draws <- coefficient_draws[
        ,
        sample(seq_len(ncol(coefficient_draws)), n_samples),
        drop = FALSE
      ]
    }
    eta_draws <- t(
      as.matrix(bdlnm_data[, terms, drop = FALSE]) %*% coefficient_draws
    )
    if (type_use != "link") {
      eta_draws <- apply_linkinv_matrix(eta_draws, get_linkinv(fitted_family))
    }
    if (output == "samples") return(format_samples(eta_draws))
    return(format_summary(eta_draws))
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
    if (is_gamm_object(model)) return(list(
      beta = stats::coef(model$gam),
      vcov = as.matrix(stats::vcov(model$gam))
    ))
    if (inherits(model, "lme")) return(list(
      beta = nlme::fixef(model),
      vcov = as.matrix(stats::vcov(model))
    ))
    if (inherits(model, "gls") || inherits(model, "gam")) return(list(
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
    beta <- stats::coef(model)
    covariance_matrix <- tryCatch(
      as.matrix(stats::vcov(model)),
      error = function(e) NULL
    )
    if (!is.numeric(beta) || is.null(covariance_matrix)) {
      stop("Could not extract coefficients and covariance matrix from the model.")
    }
    list(beta = beta, vcov = covariance_matrix)
  }

  coefficient_information <- extract_coef_vcov(fit)
  beta_hat <- coefficient_information$beta
  covariance_hat <- coefficient_information$vcov

  if (is.null(names(beta_hat))) {
    stop("The fitted coefficient vector has no names.")
  }

  fixed_data <- ensure_intercept(newdata, names(beta_hat))
  fixed_terms <- expected_fixed_terms(names(beta_hat), fixed_data)
  check_vcov_names(covariance_hat, fixed_terms)

  beta_hat <- beta_hat[fixed_terms]
  covariance_hat <- covariance_hat[fixed_terms, fixed_terms, drop = FALSE]

  if (any(!is.finite(beta_hat)) || any(!is.finite(covariance_hat))) {
    stop(
      "Non-finite coefficient or covariance values were detected. ",
      "Check model convergence before requesting uncertainty."
    )
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for frequentist uncertainty.")
  }

  beta_draws <- MASS::mvrnorm(
    n = n_samples,
    mu = beta_hat,
    Sigma = covariance_hat
  )
  if (is.null(dim(beta_draws))) {
    beta_draws <- matrix(beta_draws, nrow = 1L)
  }
  colnames(beta_draws) <- fixed_terms

  draw_sd <- apply(beta_draws, 2, stats::sd)
  if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
    warning(
      "Near-zero coefficient-draw variability was detected. Prediction intervals may collapse.",
      call. = FALSE
    )
  }

  prediction_draws <- beta_draws %*%
    t(as.matrix(fixed_data[, fixed_terms, drop = FALSE]))

  if (re != "population" &&
      (inherits(fit, "glmmTMB") || inherits(fit, "merMod") ||
       inherits(fit, "lme"))) {
    warning(
      "Frequentist mixed-model uncertainty currently reflects fixed effects only.",
      call. = FALSE
    )
  }

  if (type_use != "link") {
    prediction_draws <- apply_linkinv_matrix(
      prediction_draws,
      get_linkinv(fitted_family)
    )
  }

  if (output == "samples") return(format_samples(prediction_draws))
  format_summary(prediction_draws)
}
