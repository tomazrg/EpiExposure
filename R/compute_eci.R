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
#' of the simulated distribution on the requested `scale`.
#'
#' @param profile Numeric vector of exposure values representing a lag profile.
#'   Its length must be equal to `lag_max + 1` for the exposure variable stored
#'   in the fitted model.
#' @param fit Fitted model from `fit_epidlnm()`. Required to compute
#'   `ECI_weighted`. If `NULL`, only `ECI_raw` is returned.
#' @param scale Character. Scale for `ECI_weighted`:
#'   - `"link"`: linear predictor scale (default)
#'   - `"response"`: inverse-link transformed scale
#'   - `"percent"`: relative change scale, computed as `(exp(eta) - 1) * 100`
#' @param reverse Logical. If `TRUE`, reverses the input profile before analysis.
#'   Use this when the supplied profile is ordered from oldest to most recent
#'   (e.g., `dpp` order). If `FALSE`, the function assumes the first value
#'   corresponds to lag 0 (most recent) and the last value to lag max (oldest).
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
#'   - `ECI_weighted` (median-based central estimate, on requested scale)
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
#' For `scale = "response"`, the inverse link is obtained from the fitted model
#' metadata when available.
#'
#' For `scale = "percent"`, the transformation is defined as:
#' \deqn{(exp(\eta) - 1) * 100}
#' which is most directly interpretable for log-linked or relative-effect contexts.
#'
#' @export
compute_eci <- function(
    profile,
    fit = NULL,
    scale = c("link", "response", "percent"),
    reverse = FALSE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  scale <- match.arg(scale)
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

  if (!is.logical(reverse) || length(reverse) != 1L) {
    stop("`reverse` must be TRUE or FALSE.")
  }

  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }

  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }

  n_samples <- as.integer(n_samples)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1) return(0)
    stats::sd(x)
  }

  # ----------------------------------------------------------
  # ✅ helper: matching robusto de nomes para draws do brms
  # ----------------------------------------------------------
  match_brms_draw_names <- function(cb_names_ref, draw_colnames) {

    bnames <- paste0("b_", cb_names_ref)

    if (all(bnames %in% draw_colnames)) {
      return(bnames)
    }

    # fallback defensivo
    raw_match <- cb_names_ref[cb_names_ref %in% draw_colnames]
    if (length(raw_match) == length(cb_names_ref)) {
      return(raw_match)
    }

    stop("Could not match brms posterior draw names to crossbasis columns.")
  }

  # ----------------------------------------------------------
  # Orientation handling (message always shown)
  # ----------------------------------------------------------
  if (reverse) {
    profile <- rev(profile)
    message(
      "Profile reversed: interpreting input as oldest \u2192 most recent. ",
      "After reversal, first value = lag 0 (most recent); last value = lag max (oldest)."
    )
  } else {
    message(
      "Assuming profile is lag-ordered: first value = lag 0 (most recent); ",
      "last value = lag max (oldest)."
    )
  }

  # ----------------------------------------------------------
  # 1) ECI RAW
  # ----------------------------------------------------------
  eci_raw <- sum(profile)

  # ----------------------------------------------------------
  # No model -> return RAW only
  # ----------------------------------------------------------
  if (is.null(fit)) {
    return(data.frame(
      ECI_raw = eci_raw,
      ECI_weighted = NA_real_,
      scale = scale
    ))
  }

  # ----------------------------------------------------------
  # Extract model metadata
  # ----------------------------------------------------------
  spec <- attr(fit, "epiexposure_spec")
  vars <- attr(fit, "epiexposure_vars")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")
  family_fit <- attr(fit, "epiexposure_family")

  if (is.null(spec)) {
    stop("`fit` does not contain `epiexposure_spec`.")
  }

  if (length(vars) != 1) {
    stop("`compute_eci()` currently supports a single exposure variable.")
  }

  var <- vars[1]
  spec_v <- spec[[var]]

  if (is.null(spec_v)) {
    stop("Missing spec for variable: ", var)
  }

  # ----------------------------------------------------------
  # Resolve lag dimension
  # ----------------------------------------------------------
  lag_max_use <- as.integer(max(spec_v$lag_max))

  if (length(profile) != lag_max_use + 1L) {
    stop("`profile` length must be equal to lag_max + 1.")
  }

  # ----------------------------------------------------------
  # Rebuild crossbasis and extract final row
  # ----------------------------------------------------------
  cb <- dlnm::crossbasis(
    profile,
    lag = lag_max_use,
    argvar = spec_v$argvar,
    arglag = spec_v$arglag
  )

  cb_row <- as.numeric(cb[lag_max_use + 1L, ])

  # ----------------------------------------------------------
  # Helper: order coefficient names robustly
  # ----------------------------------------------------------
  sort_cb_names <- function(x) {
    if (length(x) == 0) return(x)

    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    idx[is.na(idx)] <- seq_along(x)

    x[order(idx)]
  }

  # ----------------------------------------------------------
  # Helper: inverse link
  # ----------------------------------------------------------
  get_linkinv <- function(family_fit) {

    if (is.character(family_fit)) {
      if (family_fit %in% c("beta", "binomial")) return(plogis)
      if (family_fit %in% c("poisson", "gamma")) return(exp)
      if (family_fit %in% c("gaussian")) return(identity)
    }

    if (inherits(family_fit, "family") && !is.null(family_fit$linkinv)) {
      return(family_fit$linkinv)
    }

    identity
  }

  transform_scale <- function(x, scale, linkinv) {
    if (scale == "link") return(x)
    if (scale == "response") return(linkinv(x))
    if (scale == "percent") return((exp(x) - 1) * 100)
    stop("Unsupported scale.")
  }

  # ----------------------------------------------------------
  # Helper: deterministic coefficients + vcov
  # ----------------------------------------------------------
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
      V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
      if (is.null(V)) stop("Could not extract vcov from spaMM model.")
      return(list(
        beta = spaMM::fixef(model),
        vcov = V
      ))
    }

    if (inherits(model, "brmsfit")) {
      b <- brms::fixef(model)
      beta <- b[, "Estimate"]
      names(beta) <- rownames(b)
      return(list(
        beta = beta,
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      V <- diag(model$summary.fixed$sd^2)
      rownames(V) <- names(beta)
      colnames(V) <- names(beta)
      return(list(beta = beta, vcov = V))
    }

    if (inherits(model, "bdlnm")) {
      beta <- model$coefficients.summary[, "mean"]
      V <- stats::cov(t(model$coefficients))
      return(list(beta = beta, vcov = V))
    }

    return(list(
      beta = stats::coef(model),
      vcov = as.matrix(stats::vcov(model))
    ))
  }

  # ----------------------------------------------------------
  # Helper: coefficient names for cb terms
  # ----------------------------------------------------------
  get_cb_names <- function(coef_names, cb_cols_fit = NULL) {

    cb_names <- grep(paste0("^cb_", var, "_"), coef_names, value = TRUE)

    if (length(cb_names) == 0 && !is.null(cb_cols_fit)) {
      cb_names <- intersect(cb_cols_fit, coef_names)
    }

    cb_names <- sort_cb_names(cb_names)

    cb_names
  }

  # ----------------------------------------------------------
  # Build stable coefficient reference
  # ----------------------------------------------------------
  cv <- extract_coef_vcov(fit)
  beta_full <- cv$beta

  cb_names_ref <- get_cb_names(names(beta_full), cb_cols_fit)

  if (length(cb_names_ref) != length(cb_row)) {
    stop("Mismatch between coefficient vector and crossbasis structure.")
  }

  linkinv <- get_linkinv(family_fit)

  # ----------------------------------------------------------
  # Deterministic (no uncertainty)
  # ----------------------------------------------------------
  if (!uncertainty) {

    beta_cb <- beta_full[cb_names_ref]
    eta_weighted <- sum(cb_row * beta_cb)

    eci_weighted <- transform_scale(
      x = eta_weighted,
      scale = scale,
      linkinv = linkinv
    )

    return(data.frame(
      ECI_raw = eci_raw,
      ECI_weighted = eci_weighted,
      scale = scale
    ))
  }

  # ----------------------------------------------------------
  # Helper: posterior / simulated draws by engine
  # returns matrix n_draws x p
  # ----------------------------------------------------------
  extract_beta_draws <- function(model, cb_names_ref, n_samples) {

    # ---------- brms: REAL posterior draws ----------
    if (inherits(model, "brmsfit")) {

      draws <- as.matrix(brms::as_draws_matrix(model))
      draw_names <- match_brms_draw_names(cb_names_ref, colnames(draws))

      draws <- draws[, draw_names, drop = FALSE]

      if (nrow(draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(nrow(draws)), n_samples)
        draws <- draws[keep, , drop = FALSE]
      }

      colnames(draws) <- cb_names_ref
      return(draws)
    }

    # ---------- INLA: REAL posterior draws ----------
    if (inherits(model, "inla")) {

      if (!requireNamespace("INLA", quietly = TRUE)) {
        stop("Package 'INLA' is required for INLA uncertainty quantification.")
      }

      posterior <- tryCatch(
        INLA::inla.posterior.sample(n = n_samples, result = model),
        error = function(e) NULL
      )

      if (is.null(posterior)) {
        stop("INLA posterior samples could not be drawn. Ensure the model was fitted with control.compute = list(config = TRUE).")
      }

      beta_draws <- do.call(rbind, lapply(posterior, function(s) {
        latent <- s$latent
        names(latent) <- gsub(":1$", "", names(latent))
        vals <- latent[cb_names_ref]
        as.numeric(vals)
      }))

      colnames(beta_draws) <- cb_names_ref
      return(beta_draws)
    }

    # ---------- bdlnm: REAL posterior draws ----------
    if (inherits(model, "bdlnm")) {

      beta_draws <- model$coefficients

      if (is.null(dim(beta_draws))) {
        beta_draws <- matrix(beta_draws, ncol = 1)
      }

      if (!all(cb_names_ref %in% rownames(beta_draws))) {
        stop("Mismatch between bdlnm posterior draws and crossbasis structure.")
      }

      beta_draws <- beta_draws[cb_names_ref, , drop = FALSE]

      if (ncol(beta_draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(ncol(beta_draws)), n_samples)
        beta_draws <- beta_draws[, keep, drop = FALSE]
      }

      out <- t(beta_draws)
      colnames(out) <- cb_names_ref
      return(out)
    }

    # ---------- Frequentist: normal approximation ----------
    beta_hat <- cv$beta[cb_names_ref]
    V_hat <- cv$vcov[cb_names_ref, cb_names_ref, drop = FALSE]

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

    colnames(beta_draws) <- cb_names_ref
    return(beta_draws)
  }

  # ----------------------------------------------------------
  # Uncertainty
  # ----------------------------------------------------------
  beta_draws <- extract_beta_draws(
    model = fit,
    cb_names_ref = cb_names_ref,
    n_samples = n_samples
  )

  draw_sd <- apply(beta_draws, 2, stats::sd)
  if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
    warning(
      "Near-zero variability detected in coefficient draws for variable '", var,
      "'. ECI intervals may collapse to a single value."
    )
  }

  # eta draws
  eta_draws <- as.numeric(beta_draws %*% cb_row)

  # transformação na escala solicitada draw-by-draw
  eci_draws <- transform_scale(
    x = eta_draws,
    scale = scale,
    linkinv = linkinv
  )

  # ----------------------------------------------------------
  # OUTPUT
  # ----------------------------------------------------------
  if (output == "samples") {
    return(data.frame(
      sample = seq_along(eci_draws),
      ECI_raw = eci_raw,
      ECI_weighted = eci_draws,
      scale = scale
    ))
  }

  return(data.frame(
    ECI_raw = eci_raw,
    ECI_weighted = stats::median(eci_draws, na.rm = TRUE),
    sd = safe_sd(eci_draws),
    lower = safe_quantile(eci_draws)[1],
    upper = safe_quantile(eci_draws)[2],
    scale = scale
  ))
}
