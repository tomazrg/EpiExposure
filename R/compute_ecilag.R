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

  lag_max <- as.integer(spec_v$lag_max)
  if (length(profile) != lag_max + 1L) {
    stop("`profile` length must be equal to lag_max + 1.")
  }

  ECI_raw <- sum(profile)

  # -----------------------
  # Rebuild crossbasis row from profile
  # -----------------------
  cb_from_profile <- function(p) {
    cb <- dlnm::crossbasis(
      p,
      lag = lag_max,
      argvar = spec_v$argvar,
      arglag = spec_v$arglag
    )
    as.numeric(cb[lag_max + 1L, ])
  }

  cb_row0 <- cb_from_profile(profile)

  # -----------------------
  # Match cb coefficient names robustly
  # -----------------------
  get_cb_names <- function(coef_names = NULL, row_names = NULL, p, var, cb_cols_fit = NULL) {

    nm <- character(0)

    if (!is.null(coef_names)) {
      nm <- grep(paste0("^cb_", var, "_"), coef_names, value = TRUE)
    }

    if (length(nm) == 0 && !is.null(row_names)) {
      nm <- intersect(row_names, paste0("cb_", var, "_", seq_len(p)))
    }

    if (length(nm) == 0 && !is.null(cb_cols_fit)) {
      nm <- intersect(cb_cols_fit, coef_names %||% row_names %||% character(0))
    }

    nm
  }

  # -----------------------
  # Extract cb coefficients (deterministic)
  # -----------------------
  extract_beta_cb <- function(model, var, p, cb_cols_fit = NULL) {

    if (inherits(model, "glmmTMB")) {
      beta <- glmmTMB::fixef(model)$cond
      nm <- get_cb_names(coef_names = names(beta), p = p, var = var, cb_cols_fit = cb_cols_fit)
      if (length(nm) != p) stop("Mismatch between glmmTMB coefficients and cb dimension.")
      return(beta[nm])
    }

    if (inherits(model, "merMod")) {
      beta <- lme4::fixef(model)
      nm <- get_cb_names(coef_names = names(beta), p = p, var = var, cb_cols_fit = cb_cols_fit)
      if (length(nm) != p) stop("Mismatch between merMod coefficients and cb dimension.")
      return(beta[nm])
    }

    if (inherits(model, "glm") || inherits(model, "gam") || inherits(model, "gls") || inherits(model, "lme")) {
      beta <- stats::coef(model)
      nm <- get_cb_names(coef_names = names(beta), p = p, var = var, cb_cols_fit = cb_cols_fit)
      if (length(nm) != p) stop("Mismatch between coefficients and cb dimension.")
      return(beta[nm])
    }

    if (inherits(model, "HLfit")) {
      beta <- spaMM::fixef(model)
      nm <- get_cb_names(coef_names = names(beta), p = p, var = var, cb_cols_fit = cb_cols_fit)
      if (length(nm) != p) stop("Mismatch between spaMM coefficients and cb dimension.")
      return(beta[nm])
    }

    if (inherits(model, "brmsfit")) {
      fe <- brms::fixef(model)
      beta <- fe[, "Estimate"]
      names(beta) <- rownames(fe)
      nm <- grep(paste0("^b_cb_", var, "_"), names(beta), value = TRUE)
      if (length(nm) == 0) nm <- grep(paste0("^b_.*", var), names(beta), value = TRUE)
      if (length(nm) != p) stop("Mismatch between brms coefficients and cb dimension.")
      return(beta[nm])
    }

    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      nm <- get_cb_names(coef_names = names(beta), p = p, var = var, cb_cols_fit = cb_cols_fit)
      if (length(nm) != p) stop("Mismatch between INLA coefficients and cb dimension.")
      return(beta[nm])
    }

    if (inherits(model, "bdlnm")) {
      beta <- model$coefficients.summary[, "mean"]
      nm <- get_cb_names(coef_names = names(beta), p = p, var = var, cb_cols_fit = cb_cols_fit)
      if (length(nm) == 0) {
        nm <- get_cb_names(row_names = names(beta), p = p, var = var, cb_cols_fit = cb_cols_fit)
      }
      if (length(nm) != p) stop("Mismatch between bdlnm coefficients and cb dimension.")
      return(beta[nm])
    }

    stop("Unsupported model class for lag-specific ECI decomposition.")
  }

  # -----------------------
  # Extract cb posterior / simulated draws
  # returns matrix n_samples x p
  # -----------------------
  extract_beta_cb_draws <- function(model, var, p, n_samples, cb_cols_fit = NULL) {

    # Frequentist: asymptotic approximation
    if (!inherits(model, c("brmsfit", "inla", "bdlnm"))) {

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

      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for frequentist uncertainty approximation.")
      }

      cv <- extract_coef_vcov(model)
      beta_hat <- cv$beta
      V_hat <- cv$vcov

      nm <- get_cb_names(coef_names = names(beta_hat), p = p, var = var, cb_cols_fit = cb_cols_fit)
      if (length(nm) != p) stop("Mismatch between frequentist coefficients and cb dimension.")

      beta_hat <- beta_hat[nm]
      V_hat <- V_hat[nm, nm, drop = FALSE]

      draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat, Sigma = V_hat)
      if (is.null(dim(draws))) {
        draws <- matrix(draws, nrow = 1)
      }
      colnames(draws) <- nm
      return(draws)
    }

    # brms: REAL posterior draws
    if (inherits(model, "brmsfit")) {
      draws <- as.matrix(brms::as_draws_matrix(model))
      nm <- grep(paste0("^b_cb_", var, "_"), colnames(draws), value = TRUE)
      if (length(nm) == 0) nm <- grep(paste0("^b_.*", var), colnames(draws), value = TRUE)
      if (length(nm) != p) stop("Could not match brms draws to cb dimension.")

      if (nrow(draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(nrow(draws)), n_samples)
        draws <- draws[keep, nm, drop = FALSE]
      } else {
        draws <- draws[, nm, drop = FALSE]
      }
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

      cb_names <- paste0("cb_", var, "_", seq_len(p))

      draws <- do.call(rbind, lapply(posterior, function(s) {
        latent <- s$latent
        names(latent) <- gsub(":1$", "", names(latent))
        as.numeric(latent[cb_names])
      }))
      colnames(draws) <- cb_names
      return(draws)
    }

    # bdlnm: REAL posterior draws
    if (inherits(model, "bdlnm")) {
      beta_draws <- model$coefficients

      if (is.null(dim(beta_draws))) {
        beta_draws <- matrix(beta_draws, ncol = 1)
      }

      nm <- get_cb_names(
        row_names = rownames(beta_draws),
        p = p,
        var = var,
        cb_cols_fit = cb_cols_fit
      )

      if (length(nm) != p) stop("Could not match bdlnm draws to cb dimension.")

      if (ncol(beta_draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(ncol(beta_draws)), n_samples)
        beta_draws <- beta_draws[nm, keep, drop = FALSE]
      } else {
        beta_draws <- beta_draws[nm, , drop = FALSE]
      }

      return(t(beta_draws))
    }

    stop("Unsupported model class for uncertainty in `compute_ecilag()`.")
  }

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

    beta_cb <- extract_beta_cb(fit, var, length(cb_row0), cb_cols_fit = cb_cols_fit)

    eta0 <- eta_from_profile_beta(profile, beta_cb)
    ECI_weighted <- eta0

    lags <- 0:lag_max
    w <- numeric(length(lags))

    if (center) {
      for (i in seq_along(lags)) {
        p_plus <- profile
        p_minus <- profile
        p_plus[i]  <- p_plus[i]  + eps
        p_minus[i] <- p_minus[i] - eps
        w[i] <- (eta_from_profile_beta(p_plus, beta_cb) - eta_from_profile_beta(p_minus, beta_cb)) / (2 * eps)
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
  # Uncertainty case
  # -----------------------
  beta_draws <- extract_beta_cb_draws(
    model = fit,
    var = var,
    p = length(cb_row0),
    n_samples = n_samples,
    cb_cols_fit = cb_cols_fit
  )

  n_draws <- nrow(beta_draws)
  lags <- 0:lag_max

  eta_draws <- numeric(n_draws)
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
        w[i] <- (eta_from_profile_beta(p_plus, beta_cb) - eta_from_profile_beta(p_minus, beta_cb)) / (2 * eps)
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

  by_lag_samples <- do.call(rbind, samples_list)

  by_lag_summary <- by_lag_samples |>
    dplyr::group_by(var, lag, exposure) |>
    dplyr::summarise(
      weight = stats::median(weight, na.rm = TRUE),
      weight_sd = stats::sd(weight, na.rm = TRUE),
      weight_lower = safe_quantile(weight)[1],
      weight_upper = safe_quantile(weight)[2],

      contribution = stats::median(contribution, na.rm = TRUE),
      contribution_sd = stats::sd(contribution, na.rm = TRUE),
      contribution_lower = safe_quantile(contribution)[1],
      contribution_upper = safe_quantile(contribution)[2],

      percent_contribution = stats::median(percent_contribution, na.rm = TRUE),
      percent_contribution_sd = stats::sd(percent_contribution, na.rm = TRUE),
      percent_contribution_lower = safe_quantile(percent_contribution)[1],
      percent_contribution_upper = safe_quantile(percent_contribution)[2],

      abs_contribution = if ("abs_contribution" %in% names(by_lag_samples)) stats::median(abs_contribution, na.rm = TRUE) else NA_real_,
      abs_contribution_sd = if ("abs_contribution" %in% names(by_lag_samples)) stats::sd(abs_contribution, na.rm = TRUE) else NA_real_,
      abs_contribution_lower = if ("abs_contribution" %in% names(by_lag_samples)) safe_quantile(abs_contribution)[1] else NA_real_,
      abs_contribution_upper = if ("abs_contribution" %in% names(by_lag_samples)) safe_quantile(abs_contribution)[2] else NA_real_,

      abs_weight = if ("abs_weight" %in% names(by_lag_samples)) stats::median(abs_weight, na.rm = TRUE) else NA_real_,
      abs_weight_sd = if ("abs_weight" %in% names(by_lag_samples)) stats::sd(abs_weight, na.rm = TRUE) else NA_real_,
      abs_weight_lower = if ("abs_weight" %in% names(by_lag_samples)) safe_quantile(abs_weight)[1] else NA_real_,
      abs_weight_upper = if ("abs_weight" %in% names(by_lag_samples)) safe_quantile(abs_weight)[2] else NA_real_,

      .groups = "drop"
    )

  res <- list(
    ECI_raw = ECI_raw,
    ECI_weighted = stats::median(eta_draws, na.rm = TRUE),
    ECI_weighted_sd = stats::sd(eta_draws, na.rm = TRUE),
    ECI_weighted_lower = safe_quantile(eta_draws)[1],
    ECI_weighted_upper = safe_quantile(eta_draws)[2],
    by_lag = by_lag_summary
  )

  if (output == "samples") {
    res$by_lag_samples <- by_lag_samples
    res$ECI_weighted_samples <- data.frame(
      sample = seq_len(length(eta_draws)),
      ECI_weighted = eta_draws
    )
  }

  return(res)
}
