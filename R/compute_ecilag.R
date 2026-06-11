#' Compute lag-specific decomposition of Exposure Cumulative Impact (ECI)
#'
#' Decomposes the weighted Exposure Cumulative Impact (ECI) into lag-specific
#' contributions using numerical derivatives of the DLNM linear predictor.
#'
#' This function returns:
#' - `ECI_raw`: the unweighted cumulative exposure
#' - `ECI_weighted`: the weighted cumulative impact on the model scale
#' - `by_lag`: lag-specific contribution summaries
#'
#' The lag-specific contribution is computed as:
#'
#' \deqn{
#' \eta(x) = \sum_k \beta_k B_k(x)
#' }
#'
#' \deqn{
#' w_l(x) \approx \frac{\eta(x + \varepsilon e_l) - \eta(x)}{\varepsilon}
#' }
#'
#' \deqn{
#' contribution_l = x_l \times w_l
#' }
#'
#' When `center = TRUE`, a central finite difference is used; otherwise a
#' forward difference is applied.
#'
#' This function supports both deterministic estimation and uncertainty
#' propagation. When `uncertainty = TRUE`, lag-specific contributions are
#' recomputed across simulated or posterior draws of the model coefficients.
#'
#' If `output = "summary"`, the central estimate is computed as the median of
#' the simulated distributions, and interval limits are derived from empirical
#' quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** when `uncertainty = TRUE`, all central estimates returned in
#' `ECI_weighted` and `by_lag` are based on the median of the simulated
#' distribution.
#'
#' @param profile Numeric vector, length = `lag_max + 1`, ordered as lag
#'   `0, 1, ..., lag_max`.
#' @param fit Fitted model returned by `fit_epidlnm()`. The fitted model must
#'   contain `epiexposure_spec` metadata.
#' @param eps Small numeric perturbation used for finite differences
#'   (default = `1e-6`).
#' @param center Logical. If `TRUE`, uses central difference; otherwise uses
#'   forward difference.
#' @param absolute Logical. If `TRUE`, also reports absolute contributions
#'   and absolute weights.
#' @param reverse Logical. If `TRUE`, reverses the input profile before analysis.
#'   Use this when the supplied profile is ordered from oldest to most recent
#'   (e.g., `dpp` order). If `FALSE`, the function assumes the first value
#'   corresponds to lag 0 (most recent) and the last value to lag max (oldest).
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Integer. Number of samples used for uncertainty propagation.
#'
#' @return A list with:
#' - `ECI_raw`
#' - `ECI_weighted`
#' - `by_lag` (always returned as a summary-style table)
#'
#' If `uncertainty = TRUE`, the list additionally includes:
#' - `ECI_weighted_sd`
#' - `ECI_weighted_lower`
#' - `ECI_weighted_upper`
#'
#' If `uncertainty = TRUE` and `output = "samples"`, the list also contains:
#' - `by_lag_samples`
#' - `ECI_weighted_samples`
#'
#' @details
#' Uncertainty is propagated using model-consistent sampling:
#' - Bayesian models (e.g., `brms`, `INLA`, `bdlnm`) use posterior draws
#' - Frequentist models use simulation from the asymptotic coefficient distribution
#'
#' For summary outputs under uncertainty, the median is used instead of the
#' mean to provide a more robust central estimate under asymmetric
#' distributions.
#'
#' @export
compute_ecilag <- function(
    profile,
    fit,
    eps = 1e-6,
    center = TRUE,
    absolute = TRUE,
    reverse = FALSE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  output <- match.arg(output)

  # -----------------------
  # Validations
  # -----------------------
  if (!is.numeric(profile) || any(!is.finite(profile))) {
    stop("`profile` must be a finite numeric vector.")
  }
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a positive numeric scalar.")
  }
  if (!is.logical(center) || length(center) != 1) {
    stop("`center` must be TRUE/FALSE.")
  }
  if (!is.logical(absolute) || length(absolute) != 1) {
    stop("`absolute` must be TRUE/FALSE.")
  }
  if (!is.logical(reverse) || length(reverse) != 1) {
    stop("`reverse` must be TRUE/FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1) {
    stop("`uncertainty` must be TRUE/FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1 || !is.finite(n_samples) || n_samples <= 0) {
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

  spec <- attr(fit, "epiexposure_spec")
  vars <- attr(fit, "epiexposure_vars")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")

  if (is.null(spec) || !is.list(spec)) {
    stop("`fit` is missing `epiexposure_spec`.")
  }
  if (is.null(vars) || length(vars) != 1) {
    stop("`compute_ecilag()` currently supports a single exposure variable fit.")
  }

  var <- vars[1]
  spec_v <- spec[[var]]
  if (is.null(spec_v)) {
    stop("Missing spec for variable: ", var)
  }

  # =========================================================
  # ✅ AJUSTE 1 — DIMENSÃO TEMPORAL
  # =========================================================
  lag_max_use <- as.integer(max(spec_v$lag_max))

  if (length(profile) != lag_max_use + 1L) {
    stop("`profile` length must be equal to lag_max + 1.")
  }

  ECI_raw <- sum(profile)

  # -----------------------
  # Rebuild crossbasis row from profile
  # -----------------------
  cb_from_profile <- function(p) {
    cb <- dlnm::crossbasis(
      p,
      lag = lag_max_use,
      argvar = spec_v$argvar,
      arglag = spec_v$arglag
    )
    as.numeric(cb[lag_max_use + 1L, ])
  }

  cb_row0 <- cb_from_profile(profile)

  # =========================================================
  # ✅ AJUSTE 4 — ORDEM DOS COEFICIENTES
  # =========================================================
  sort_cb_names <- function(x) {
    if (length(x) == 0) return(x)
    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    idx[is.na(idx)] <- seq_along(x)
    x[order(idx)]
  }

  # -----------------------
  # Match cb coefficient names robustly
  # -----------------------
  get_cb_names <- function(names_vec, var, p, cb_cols_fit = NULL) {

    nm <- grep(paste0("^cb_", var, "_"), names_vec, value = TRUE)

    # brms fixed effects / draws
    if (length(nm) == 0) {
      nm <- grep(paste0("^b_cb_", var, "_"), names_vec, value = TRUE)
    }

    if (length(nm) == 0 && !is.null(cb_cols_fit)) {
      nm <- intersect(cb_cols_fit, names_vec)
    }

    if (length(nm) == 0 && !is.null(cb_cols_fit)) {
      nm <- intersect(paste0("b_", cb_cols_fit), names_vec)
    }

    nm <- sort_cb_names(nm)

    if (length(nm) != p) {
      stop("Could not match coefficient names to crossbasis dimension.")
    }

    nm
  }

  # -----------------------
  # Extract coefficients / vcov
  # -----------------------
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

    if (inherits(model, "glm") || inherits(model, "gam") || inherits(model, "gls") || inherits(model, "lme")) {
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

    if (inherits(model, "brmsfit")) {
      fe <- brms::fixef(model)
      beta <- fe[, "Estimate"]
      names(beta) <- rownames(fe)
      return(list(beta = beta, vcov = as.matrix(stats::vcov(model))))
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

    stop("Unsupported model class for lag-specific ECI decomposition.")
  }

  cv <- extract_coef_vcov(fit)
  beta_full <- cv$beta
  vcov_full <- cv$vcov

  cb_names_ref <- get_cb_names(
    names_vec = names(beta_full),
    var = var,
    p = length(cb_row0),
    cb_cols_fit = cb_cols_fit
  )

  # -----------------------
  # eta(profile | beta_cb)
  # -----------------------
  eta_from_profile_beta <- function(p, beta_cb) {
    cb_row <- cb_from_profile(p)
    if (length(cb_row) != length(beta_cb)) {
      stop("Mismatch between crossbasis dimension and cb coefficients.")
    }
    sum(cb_row * beta_cb)
  }

  # -----------------------
  # Deterministic case
  # -----------------------
  if (!uncertainty) {

    beta_cb <- beta_full[cb_names_ref]

    eta0 <- eta_from_profile_beta(profile, beta_cb)
    ECI_weighted <- eta0

    lags <- 0:lag_max_use
    w <- numeric(length(lags))

    if (center) {
      for (i in seq_along(lags)) {
        p_plus <- profile
        p_minus <- profile
        p_plus[i]  <- p_plus[i]  + eps
        p_minus[i] <- p_minus[i] - eps
        w[i] <- (eta_from_profile_beta(p_plus, beta_cb) -
                   eta_from_profile_beta(p_minus, beta_cb)) / (2 * eps)
      }
    } else {
      for (i in seq_along(lags)) {
        p_plus <- profile
        p_plus[i] <- p_plus[i] + eps
        w[i] <- (eta_from_profile_beta(p_plus, beta_cb) - eta0) / eps
      }
    }

    contribution <- profile * w

    by_lag <- data.frame(
      var = var,
      lag = lags,
      exposure = profile,
      weight = w,
      contribution = contribution,
      stringsAsFactors = FALSE
    )

    if (absolute) {
      by_lag$abs_contribution <- abs(by_lag$contribution)
      by_lag$abs_weight <- abs(by_lag$weight)
    }

    denom <- sum(abs(by_lag$contribution), na.rm = TRUE)
    by_lag$percent_contribution <- if (denom > 0) 100 * abs(by_lag$contribution) / denom else NA_real_

    return(list(
      ECI_raw = ECI_raw,
      ECI_weighted = ECI_weighted,
      by_lag = by_lag
    ))
  }

  # -----------------------
  # Extract posterior / simulated draws
  # returns matrix n_samples x p
  # -----------------------
  extract_beta_cb_draws <- function(model, cb_names_ref, n_samples) {

    # Frequentist: asymptotic approximation
    if (!inherits(model, c("brmsfit", "inla", "bdlnm"))) {

      beta_hat <- beta_full[cb_names_ref]
      V_hat <- vcov_full[cb_names_ref, cb_names_ref, drop = FALSE]

      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for frequentist uncertainty approximation.")
      }

      draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat, Sigma = V_hat)

      if (is.null(dim(draws))) {
        draws <- matrix(draws, nrow = 1)
      }

      colnames(draws) <- cb_names_ref
      return(draws)
    }

    # brms: REAL posterior draws
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

    # INLA: REAL posterior draws
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

      draws <- do.call(rbind, lapply(posterior, function(s) {
        latent <- s$latent
        names(latent) <- gsub(":1$", "", names(latent))
        as.numeric(latent[cb_names_ref])
      }))

      colnames(draws) <- cb_names_ref
      return(draws)
    }

    # bdlnm: REAL posterior draws
    if (inherits(model, "bdlnm")) {
      beta_draws <- model$coefficients

      if (is.null(dim(beta_draws))) {
        beta_draws <- matrix(beta_draws, ncol = 1)
      }

      if (!all(cb_names_ref %in% rownames(beta_draws))) {
        stop("Could not match bdlnm draws to cb dimension.")
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

    stop("Unsupported model class for uncertainty in `compute_ecilag()`.")
  }

  # -----------------------
  # Uncertainty case
  # -----------------------
  beta_draws <- extract_beta_cb_draws(
    model = fit,
    cb_names_ref = cb_names_ref,
    n_samples = n_samples
  )

  draw_sd <- apply(beta_draws, 2, stats::sd)
  if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
    warning(
      "Near-zero variability detected in coefficient draws for variable '", var,
      "'. ECI-lag intervals may collapse to a single value."
    )
  }

  n_draws <- nrow(beta_draws)
  lags <- 0:lag_max_use

  eta_draws <- numeric(n_draws)

  weight_mat <- matrix(NA_real_, nrow = length(lags), ncol = n_draws)
  contrib_mat <- matrix(NA_real_, nrow = length(lags), ncol = n_draws)
  pct_mat <- matrix(NA_real_, nrow = length(lags), ncol = n_draws)

  if (absolute) {
    abs_contrib_mat <- matrix(NA_real_, nrow = length(lags), ncol = n_draws)
    abs_weight_mat <- matrix(NA_real_, nrow = length(lags), ncol = n_draws)
  }

  samples_list <- vector("list", n_draws)

  for (s in seq_len(n_draws)) {

    beta_cb <- beta_draws[s, ]
    eta0 <- eta_from_profile_beta(profile, beta_cb)
    eta_draws[s] <- eta0

    w <- numeric(length(lags))

    if (center) {
      for (i in seq_along(lags)) {
        p_plus <- profile
        p_minus <- profile
        p_plus[i]  <- p_plus[i]  + eps
        p_minus[i] <- p_minus[i] - eps
        w[i] <- (eta_from_profile_beta(p_plus, beta_cb) -
                   eta_from_profile_beta(p_minus, beta_cb)) / (2 * eps)
      }
    } else {
      for (i in seq_along(lags)) {
        p_plus <- profile
        p_plus[i] <- p_plus[i] + eps
        w[i] <- (eta_from_profile_beta(p_plus, beta_cb) - eta0) / eps
      }
    }

    contribution <- profile * w
    denom <- sum(abs(contribution), na.rm = TRUE)
    pct_contrib <- if (denom > 0) 100 * abs(contribution) / denom else rep(NA_real_, length(contribution))

    weight_mat[, s] <- w
    contrib_mat[, s] <- contribution
    pct_mat[, s] <- pct_contrib

    if (absolute) {
      abs_contrib_mat[, s] <- abs(contribution)
      abs_weight_mat[, s] <- abs(w)
    }

    if (output == "samples") {
      df_s <- data.frame(
        sample = s,
        var = var,
        lag = lags,
        exposure = profile,
        weight = w,
        contribution = contribution,
        percent_contribution = pct_contrib,
        stringsAsFactors = FALSE
      )

      if (absolute) {
        df_s$abs_contribution <- abs(df_s$contribution)
        df_s$abs_weight <- abs(df_s$weight)
      }

      samples_list[[s]] <- df_s
    }
  }

  # -----------------------
  # Summary directly by matrix
  # -----------------------
  by_lag_summary <- data.frame(
    var = var,
    lag = lags,
    exposure = profile,

    weight = apply(weight_mat, 1, stats::median, na.rm = TRUE),
    weight_sd = apply(weight_mat, 1, safe_sd),
    weight_lower = apply(weight_mat, 1, function(z) safe_quantile(z)[1]),
    weight_upper = apply(weight_mat, 1, function(z) safe_quantile(z)[2]),

    contribution = apply(contrib_mat, 1, stats::median, na.rm = TRUE),
    contribution_sd = apply(contrib_mat, 1, safe_sd),
    contribution_lower = apply(contrib_mat, 1, function(z) safe_quantile(z)[1]),
    contribution_upper = apply(contrib_mat, 1, function(z) safe_quantile(z)[2]),

    percent_contribution = apply(pct_mat, 1, stats::median, na.rm = TRUE),
    percent_contribution_sd = apply(pct_mat, 1, safe_sd),
    percent_contribution_lower = apply(pct_mat, 1, function(z) safe_quantile(z)[1]),
    percent_contribution_upper = apply(pct_mat, 1, function(z) safe_quantile(z)[2]),

    stringsAsFactors = FALSE
  )

  if (absolute) {
    by_lag_summary$abs_contribution <- apply(abs_contrib_mat, 1, stats::median, na.rm = TRUE)
    by_lag_summary$abs_contribution_sd <- apply(abs_contrib_mat, 1, safe_sd)
    by_lag_summary$abs_contribution_lower <- apply(abs_contrib_mat, 1, function(z) safe_quantile(z)[1])
    by_lag_summary$abs_contribution_upper <- apply(abs_contrib_mat, 1, function(z) safe_quantile(z)[2])

    by_lag_summary$abs_weight <- apply(abs_weight_mat, 1, stats::median, na.rm = TRUE)
    by_lag_summary$abs_weight_sd <- apply(abs_weight_mat, 1, safe_sd)
    by_lag_summary$abs_weight_lower <- apply(abs_weight_mat, 1, function(z) safe_quantile(z)[1])
    by_lag_summary$abs_weight_upper <- apply(abs_weight_mat, 1, function(z) safe_quantile(z)[2])
  }

  res <- list(
    ECI_raw = ECI_raw,
    ECI_weighted = stats::median(eta_draws, na.rm = TRUE),
    ECI_weighted_sd = safe_sd(eta_draws),
    ECI_weighted_lower = safe_quantile(eta_draws)[1],
    ECI_weighted_upper = safe_quantile(eta_draws)[2],
    by_lag = by_lag_summary
  )

  if (output == "samples") {
    res$by_lag_samples <- do.call(rbind, samples_list)
    res$ECI_weighted_samples <- data.frame(
      sample = seq_len(length(eta_draws)),
      ECI_weighted = eta_draws
    )
  }

  return(res)
}
