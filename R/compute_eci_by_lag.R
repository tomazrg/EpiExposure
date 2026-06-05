#' Compute ECI decomposition by lag (ECI-by-lag)
#'
#' Returns ECI_raw and ECI_weighted, plus a per-lag decomposition of the
#' weighted ECI using numerical derivatives on the DLNM linear predictor.
#'
#' Idea:
#'   eta(x) = sum_k beta_k * B_k(x)   (DLNM linear predictor contribution)
#'   w_l(x) ≈ [eta(x + eps*e_l) - eta(x)] / eps
#'   contribution_l = x_l * w_l
#'
#' @param profile Numeric vector, length = lag_max + 1, ordered as lag 0..lag_max
#' @param fit Fitted model returned by fit_epidlnm() (stores epiexposure_spec attrs)
#' @param eps Small numeric perturbation for finite differences (default 1e-6)
#' @param center If TRUE, uses central difference; else forward difference
#' @param absolute If TRUE, report absolute contributions as well
#' @param uncertainty Logical; if TRUE, quantify uncertainty
#' @param output "summary" or "samples"
#' @param n_samples Number of samples used for uncertainty quantification
#'
#' @return list with:
#'   - ECI_raw
#'   - ECI_weighted
#'   - by_lag (always summary-style; backward compatible)
#'   - by_lag_samples (only if uncertainty = TRUE and output = "samples")
#' @export
compute_eci_by_lag <- function(
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
    stop("profile must be a finite numeric vector.")
  }
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive numeric scalar.")
  }
  if (!is.logical(center) || length(center) != 1) {
    stop("center must be TRUE/FALSE.")
  }
  if (!is.logical(absolute) || length(absolute) != 1) {
    stop("absolute must be TRUE/FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1) {
    stop("uncertainty must be TRUE/FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1 || !is.finite(n_samples) || n_samples <= 0) {
    stop("n_samples must be a positive integer.")
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
    stop("fit is missing epiexposure_spec attribute.")
  }
  if (is.null(vars) || length(vars) != 1) {
    stop("compute_eci_by_lag currently supports a single exposure variable fit.")
  }

  var <- vars[1]
  spec_v <- spec[[var]]
  if (is.null(spec_v)) {
    stop("Missing spec for variable: ", var)
  }

  lag_max <- as.integer(spec_v$lag_max)
  if (length(profile) != lag_max + 1L) {
    stop("Profile length must be lag_max + 1.")
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
  # Extract cb coefficients (deterministic)
  # -----------------------
  extract_beta_cb <- function(model, var, p, cb_cols_fit = NULL) {

    # glmmTMB
    if (inherits(model, "glmmTMB")) {
      beta <- glmmTMB::fixef(model)$cond
      nm <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)
      if (length(nm) != p) stop("Mismatch between glmmTMB coefficients and cb dimension.")
      return(beta[nm])
    }

    # merMod
    if (inherits(model, "merMod")) {
      beta <- lme4::fixef(model)
      nm <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)
      if (length(nm) != p) stop("Mismatch between merMod coefficients and cb dimension.")
      return(beta[nm])
    }

    # glm / gam / gls / lme / default
    if (inherits(model, "glm") || inherits(model, "gam") || inherits(model, "gls") || inherits(model, "lme")) {
      beta <- stats::coef(model)
      nm <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)
      if (length(nm) != p) stop("Mismatch between coefficients and cb dimension.")
      return(beta[nm])
    }

    # spaMM
    if (inherits(model, "HLfit")) {
      beta <- stats::coef(model)
      nm <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)
      if (length(nm) != p) stop("Mismatch between spaMM coefficients and cb dimension.")
      return(beta[nm])
    }

    # brms
    if (inherits(model, "brmsfit")) {
      fe <- brms::fixef(model)
      beta <- fe[, "Estimate"]
      names(beta) <- rownames(fe)
      nm <- grep(paste0("^b_cb_", var, "_"), names(beta), value = TRUE)
      if (length(nm) == 0) nm <- grep(paste0("^b_.*", var), names(beta), value = TRUE)
      if (length(nm) != p) stop("Mismatch between brms coefficients and cb dimension.")
      return(beta[nm])
    }

    # INLA
    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      nm <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)
      if (length(nm) != p) stop("Mismatch between INLA coefficients and cb dimension.")
      return(beta[nm])
    }

    # bdlnm
    if (inherits(model, "bdlnm")) {
      beta <- model$coefficients.summary[, "mean"]
      nm <- intersect(names(beta), paste0("cb_", var, "_", seq_len(p)))
      if (length(nm) == 0) {
        nm <- intersect(cb_cols_fit %||% character(0), names(beta))
      }
      if (length(nm) != p) stop("Mismatch between bdlnm coefficients and cb dimension.")
      return(beta[nm])
    }

    stop("Unsupported model class for ECI-by-lag decomposition.")
  }

  # -----------------------
  # Extract cb posterior / simulated draws
  # returns matrix n_samples x p
  # -----------------------
  extract_beta_cb_draws <- function(model, var, p, n_samples, cb_cols_fit = NULL) {

    # frequentist normal approximation
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
          b <- stats::coef(model)
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

      nm <- grep(paste0("^cb_", var, "_"), names(beta_hat), value = TRUE)
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

    # brms posterior draws
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

    # INLA posterior draws
    if (inherits(model, "inla")) {
      if (!requireNamespace("INLA", quietly = TRUE)) {
        stop("Package 'INLA' is required for INLA uncertainty quantification.")
      }

      posterior <- INLA::inla.posterior.sample(n = n_samples, result = model)
      cb_names <- paste0("cb_", var, "_", seq_len(p))

      draws <- do.call(rbind, lapply(posterior, function(s) {
        latent <- s$latent
        names(latent) <- gsub(":1$", "", names(latent))
        as.numeric(latent[cb_names])
      }))
      colnames(draws) <- cb_names
      return(draws)
    }

    # bdlnm posterior draws
    if (inherits(model, "bdlnm")) {
      beta_draws <- model$coefficients

      if (is.null(dim(beta_draws))) {
        beta_draws <- matrix(beta_draws, ncol = 1)
      }

      nm <- intersect(rownames(beta_draws), paste0("cb_", var, "_", seq_len(p)))
      if (length(nm) == 0) {
        nm <- intersect(cb_cols_fit %||% character(0), rownames(beta_draws))
      }
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

    stop("Unsupported model class for uncertainty in compute_eci_by_lag().")
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

  # store results
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

  # summary preserving old column names for downstream compatibility
  by_lag_summary <- by_lag_samples |>
    dplyr::group_by(var, lag, exposure) |>
    dplyr::summarise(
      weight = mean(weight, na.rm = TRUE),
      weight_sd = stats::sd(weight, na.rm = TRUE),
      weight_lower = safe_quantile(weight)[1],
      weight_upper = safe_quantile(weight)[2],

      contribution = mean(contribution, na.rm = TRUE),
      contribution_sd = stats::sd(contribution, na.rm = TRUE),
      contribution_lower = safe_quantile(contribution)[1],
      contribution_upper = safe_quantile(contribution)[2],

      percent_contribution = mean(percent_contribution, na.rm = TRUE),
      percent_contribution_sd = stats::sd(percent_contribution, na.rm = TRUE),
      percent_contribution_lower = safe_quantile(percent_contribution)[1],
      percent_contribution_upper = safe_quantile(percent_contribution)[2],

      abs_contribution = if ("abs_contribution" %in% names(by_lag_samples)) mean(abs_contribution, na.rm = TRUE) else NA_real_,
      abs_contribution_sd = if ("abs_contribution" %in% names(by_lag_samples)) stats::sd(abs_contribution, na.rm = TRUE) else NA_real_,
      abs_contribution_lower = if ("abs_contribution" %in% names(by_lag_samples)) safe_quantile(abs_contribution)[1] else NA_real_,
      abs_contribution_upper = if ("abs_contribution" %in% names(by_lag_samples)) safe_quantile(abs_contribution)[2] else NA_real_,

      abs_weight = if ("abs_weight" %in% names(by_lag_samples)) mean(abs_weight, na.rm = TRUE) else NA_real_,
      abs_weight_sd = if ("abs_weight" %in% names(by_lag_samples)) stats::sd(abs_weight, na.rm = TRUE) else NA_real_,
      abs_weight_lower = if ("abs_weight" %in% names(by_lag_samples)) safe_quantile(abs_weight)[1] else NA_real_,
      abs_weight_upper = if ("abs_weight" %in% names(by_lag_samples)) safe_quantile(abs_weight)[2] else NA_real_,

      .groups = "drop"
    )

  res <- list(
    ECI_raw = ECI_raw,
    ECI_weighted = mean(eta_draws, na.rm = TRUE),
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
