#' Compute Exposure Cumulative Impact (ECI - the exposure profile with the fitted model coefficients)
#'
#' This function supports both deterministic estimation and uncertainty
#' propagation. When `uncertainty = TRUE`, the weighted ECI is recomputed
#' across simulated or posterior draws of the model coefficients.
#'
#' If `output = "summary"`, the central estimate is computed as the median of
#' the simulated ECI values, while interval limits are obtained from empirical
#' quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** when `uncertainty = TRUE` and `output = "summary"`,
#' `ECI_weighted` represents the *central estimate*, computed as the median
#' of the simulated distribution.
#'
#' @param profile Numeric vector of exposure values representing a lag profile.
#'   Its length must be equal to `lag_max + 1` for the exposure variable stored
#'   in the fitted model.
#' @param fit Fitted model from `fit_epidlnm()`. Required to compute
#'   `ECI_weighted`. If `NULL`, only `ECI_raw` is returned.
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Integer. Number of samples used for uncertainty quantification.
#'
#' @return A data.frame.
#'
#' - If `fit = NULL`, returns:
#'   - `ECI_raw`
#'   - `ECI_weighted = NA`
#'
#' - If `uncertainty = FALSE`, returns:
#'   - `ECI_raw`
#'   - `ECI_weighted`
#'
#' - If `uncertainty = TRUE` and `output = "summary"`, returns:
#'   - `ECI_raw`
#'   - `ECI_weighted` (median-based central estimate)
#'   - `sd`
#'   - `lower`
#'   - `upper`
#'
#' - If `uncertainty = TRUE` and `output = "samples"`, returns:
#'   - `sample`
#'   - `ECI_raw`
#'   - `ECI_weighted`
#'
#' @details
#' Uncertainty is propagated using model-consistent sampling:
#' - Bayesian models (e.g., `brms`, `INLA`, `bdlnm`) use posterior draws
#' - Frequentist models use simulation from the asymptotic coefficient distribution
#'
#' The use of the median as the central estimate improves robustness under
#' asymmetric or non-normal simulated ECI distributions.
#'
#' @export
compute_eci <- function(
    profile,
    fit = NULL,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  output <- match.arg(output)

  # -----------------------
  # Validation
  # -----------------------
  if (!is.numeric(profile)) {
    stop("`profile` must be a numeric vector.")
  }

  if (any(!is.finite(profile))) {
    stop("`profile` must contain only finite values.")
  }

  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }

  if (!is.numeric(n_samples) || length(n_samples) != 1L || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }

  n_samples <- as.integer(n_samples)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  # -----------------------
  # 1) ECI RAW
  # -----------------------
  eci_raw <- sum(profile)

  # -----------------------
  # No model → return RAW only
  # -----------------------
  if (is.null(fit)) {
    return(data.frame(
      ECI_raw = eci_raw,
      ECI_weighted = NA_real_
    ))
  }

  # -----------------------
  # Extract model metadata
  # -----------------------
  spec <- attr(fit, "epiexposure_spec")
  vars <- attr(fit, "epiexposure_vars")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")

  if (is.null(spec)) stop("`fit` does not contain `epiexposure_spec`.")
  if (length(vars) != 1) {
    stop("`compute_eci()` currently supports a single exposure variable.")
  }

  var <- vars[1]
  spec_v <- spec[[var]]

  if (is.null(spec_v)) {
    stop("Missing spec for variable: ", var)
  }

  lag_max <- as.integer(spec_v$lag_max)

  if (length(profile) != lag_max + 1L) {
    stop("`profile` length must be equal to lag_max + 1.")
  }

  # -----------------------
  # Rebuild crossbasis and extract final row
  # -----------------------
  cb <- dlnm::crossbasis(
    profile,
    lag = lag_max,
    argvar = spec_v$argvar,
    arglag = spec_v$arglag
  )

  cb_row <- as.numeric(cb[lag_max + 1L, ])

  # -----------------------
  # Helper: deterministic coefficients
  # -----------------------
  extract_coef <- function(model) {

    if (inherits(model, "glmmTMB")) {
      return(glmmTMB::fixef(model)$cond)
    }

    if (inherits(model, "merMod")) {
      return(lme4::fixef(model))
    }

    if (inherits(model, "lme")) {
      return(nlme::fixef(model))
    }

    if (inherits(model, "gls")) {
      return(stats::coef(model))
    }

    if (inherits(model, "gam")) {
      return(stats::coef(model))
    }

    if (inherits(model, "HLfit")) {
      return(spaMM::fixef(model))
    }

    if (inherits(model, "brmsfit")) {
      b <- brms::fixef(model)
      beta <- b[, "Estimate"]
      names(beta) <- rownames(b)
      return(beta)
    }

    if (inherits(model, "inla")) {
      return(fit$summary.fixed$mean)
    }

    if (inherits(model, "bdlnm")) {
      return(model$coefficients.summary[, "mean"])
    }

    return(stats::coef(model))
  }

  # -----------------------
  # Helper: coefficient names for cb terms
  # -----------------------
  get_cb_names <- function(coef_names, row_names = NULL) {
    cb_names <- grep(paste0("^cb_", var, "_"), coef_names, value = TRUE)

    if (length(cb_names) == 0 && !is.null(row_names)) {
      cb_names <- intersect(row_names, paste0("cb_", var, "_", seq_along(cb_row)))
    }

    if (length(cb_names) == 0 && !is.null(cb_cols_fit)) {
      cb_names <- intersect(cb_cols_fit, coef_names %||% row_names %||% character(0))
    }

    cb_names
  }

  # -----------------------
  # Deterministic (no uncertainty)
  # -----------------------
  if (!uncertainty) {

    beta <- extract_coef(fit)

    cb_names <- get_cb_names(names(beta))

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between basis columns and model coefficients for weighted ECI.")
    }

    beta_cb <- beta[cb_names]
    eci_weighted <- sum(cb_row * beta_cb)

    return(data.frame(
      ECI_raw = eci_raw,
      ECI_weighted = eci_weighted
    ))
  }

  # -----------------------
  # Uncertainty: coefficient draws
  # -----------------------
  if (inherits(fit, "brmsfit")) {

    # REAL posterior draws
    beta_draws <- as.matrix(brms::as_draws_matrix(fit))

    cb_names <- grep(paste0("^b_cb_", var, "_"), colnames(beta_draws), value = TRUE)
    if (length(cb_names) == 0) {
      cb_names <- grep(paste0("^b_.*", var), colnames(beta_draws), value = TRUE)
    }

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between brms posterior draws and crossbasis structure.")
    }

    if (nrow(beta_draws) > n_samples) {
      set.seed(1)
      keep <- sample(seq_len(nrow(beta_draws)), n_samples)
      beta_draws <- beta_draws[keep, cb_names, drop = FALSE]
    } else {
      beta_draws <- beta_draws[, cb_names, drop = FALSE]
    }

    eci_draws <- as.numeric(beta_draws %*% cb_row)

  } else if (inherits(fit, "inla")) {

    # REAL posterior draws
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for INLA uncertainty quantification.")
    }

    posterior <- tryCatch(
      INLA::inla.posterior.sample(n = n_samples, result = fit),
      error = function(e) NULL
    )

    if (is.null(posterior)) {
      stop("INLA posterior samples could not be drawn. Ensure the model was fitted with control.compute = list(config = TRUE).")
    }

    cb_names <- paste0("cb_", var, "_", seq_along(cb_row))

    beta_draws <- do.call(cbind, lapply(posterior, function(s) {
      latent <- s$latent
      names(latent) <- gsub(":1$", "", names(latent))
      latent[cb_names]
    }))

    if (nrow(beta_draws) != length(cb_row)) {
      stop("Mismatch between INLA posterior draws and crossbasis structure.")
    }

    eci_draws <- as.numeric(t(beta_draws) %*% cb_row)

  } else if (inherits(fit, "bdlnm")) {

    # REAL posterior draws
    beta_draws <- fit$coefficients

    if (is.null(dim(beta_draws))) {
      beta_draws <- matrix(beta_draws, ncol = 1)
    }

    cb_names <- intersect(rownames(beta_draws), paste0("cb_", var, "_", seq_along(cb_row)))
    if (length(cb_names) == 0) {
      cb_names <- intersect(cb_cols_fit %||% character(0), rownames(beta_draws))
    }

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between bdlnm posterior draws and crossbasis structure.")
    }

    if (ncol(beta_draws) > n_samples) {
      set.seed(1)
      keep <- sample(seq_len(ncol(beta_draws)), n_samples)
      beta_draws <- beta_draws[cb_names, keep, drop = FALSE]
    } else {
      beta_draws <- beta_draws[cb_names, , drop = FALSE]
    }

    eci_draws <- as.numeric(t(beta_draws) %*% cb_row)

  } else {

    # Frequentist: normal approximation
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

      if (inherits(model, "gam")) {
        return(list(
          beta = stats::coef(model),
          vcov = as.matrix(stats::vcov(model))
        ))
      }

      if (inherits(model, "HLfit")) {
        b <- spaMM::fixef(model)
        V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
        if (is.null(V)) stop("Could not extract vcov from spaMM model.")
        return(list(beta = b, vcov = V))
      }

      return(list(
        beta = stats::coef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    cv <- extract_coef_vcov(fit)

    beta_hat <- cv$beta
    V_hat <- cv$vcov

    cb_names <- get_cb_names(names(beta_hat))

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between coefficient vector and crossbasis structure.")
    }

    beta_hat <- beta_hat[cb_names]
    V_hat <- V_hat[cb_names, cb_names, drop = FALSE]

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty approximation.")
    }

    beta_draws <- MASS::mvrnorm(
      n = n_samples,
      mu = beta_hat,
      Sigma = V_hat
    )

    if (is.null(dim(beta_draws))) {
      beta_draws <- matrix(beta_draws, nrow = 1)
    }

    eci_draws <- as.numeric(beta_draws %*% cb_row)
  }

  # -----------------------
  # OUTPUT
  # -----------------------
  if (output == "samples") {
    return(data.frame(
      sample = seq_along(eci_draws),
      ECI_raw = eci_raw,
      ECI_weighted = eci_draws
    ))
  }

  return(data.frame(
    ECI_raw = eci_raw,
    ECI_weighted = stats::median(eci_draws, na.rm = TRUE),
    sd = stats::sd(eci_draws, na.rm = TRUE),
    lower = safe_quantile(eci_draws)[1],
    upper = safe_quantile(eci_draws)[2]
  ))
}
#'
#' Computes two complementary summaries for a lagged exposure profile:
#'
#' 1. `ECI_raw`: the unweighted cumulative exposure (simple sum of profile values)
#' 2. `ECI_weighted`: the DLNM-based cumulative impact, obtained by combining
