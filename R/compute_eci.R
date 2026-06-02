#' Compute Exposure Cumulative Impact (ECI)
#'
#' Computes:
#' 1) ECI_raw (sum of exposure values; unweighted)
#' 2) ECI_weighted (DLNM-based, using model coefficients)
#'
#' @param profile Numeric vector of exposure values (lag profile)
#' @param fit Fitted model from fit_epidlnm() (required for weighted ECI)
#' @param uncertainty Logical; if TRUE, quantify uncertainty
#' @param output "summary" or "samples"
#' @param n_samples Number of samples used for uncertainty quantification
#'
#' @return data.frame
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
    stop("profile must be a numeric vector.")
  }

  if (any(!is.finite(profile))) {
    stop("profile must contain only finite values.")
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

  is_bayesian_model <- function(model) {
    inherits(model, "brmsfit") || inherits(model, "inla") || inherits(model, "bdlnm")
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
  family_fit <- attr(fit, "epiexposure_family")

  if (is.null(spec)) stop("fit does not contain epiexposure_spec.")
  if (length(vars) != 1) {
    stop("compute_eci currently supports a single exposure variable.")
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

  # -----------------------
  # Rebuild crossbasis
  # -----------------------
  cb <- dlnm::crossbasis(
    profile,
    lag = lag_max,
    argvar = spec_v$argvar,
    arglag = spec_v$arglag
  )

  cb_row <- as.numeric(cb[lag_max + 1L, ])

  # -----------------------
  # Helper: linkinv
  # -----------------------
  get_linkinv <- function(model, family_fit) {

    if (!is.null(model$family) && !is.null(model$family$linkinv)) {
      return(model$family$linkinv)
    }

    if (is.character(family_fit)) {
      if (family_fit %in% c("binomial", "beta")) {
        return(stats::binomial(link = "logit")$linkinv)
      }
      if (family_fit %in% c("poisson", "gamma")) {
        return(stats::poisson(link = "log")$linkinv)
      }
      if (family_fit %in% c("gaussian")) {
        return(stats::gaussian()$linkinv)
      }
    }

    return(function(x) x)
  }

  # -----------------------
  # Deterministic (no uncertainty)
  # -----------------------
  if (!uncertainty) {

    beta <- NULL

    if (inherits(fit, "glmmTMB")) {
      beta <- glmmTMB::fixef(fit)$cond
    } else if (inherits(fit, "merMod")) {
      beta <- lme4::fixef(fit)
    } else if (inherits(fit, "glm") || inherits(fit, "gam")) {
      beta <- stats::coef(fit)
    } else if (inherits(fit, "gls")) {
      beta <- stats::coef(fit)
    } else if (inherits(fit, "brmsfit")) {
      b <- brms::fixef(fit)
      beta <- b[, "Estimate"]
      names(beta) <- rownames(b)
    } else if (inherits(fit, "inla")) {
      beta <- fit$summary.fixed$mean
    } else {
      stop("Unsupported model class for weighted ECI.")
    }

    cb_names <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between cb coefficients and basis.")
    }

    beta_cb <- beta[cb_names]

    eci_weighted <- sum(cb_row * beta_cb)

    return(data.frame(
      ECI_raw = eci_raw,
      ECI_weighted = eci_weighted
    ))
  }

  # -----------------------
  # ✅ UNCERTAINTY
  # -----------------------

  # -----------------------
  # Bayesian: brms
  # -----------------------
  if (inherits(fit, "brmsfit")) {

    beta_draws <- as.matrix(brms::as_draws_matrix(fit))

    cb_names <- grep(paste0("^b_cb_", var, "_"), colnames(beta_draws), value = TRUE)

    if (length(cb_names) == 0) {
      cb_names <- grep(paste0("^b_.*", var), colnames(beta_draws), value = TRUE)
    }

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between bdlnm/brms coefficient names and cb structure.")
    }

    beta_draws <- beta_draws[, cb_names, drop = FALSE]

    eci_draws <- as.numeric(beta_draws %*% cb_row)

    # -----------------------
    # Bayesian: INLA
    # -----------------------
  } else if (inherits(fit, "inla")) {

    posterior <- INLA::inla.posterior.sample(n = n_samples, result = fit)

    cb_names <- paste0("cb_", var, "_", seq_along(cb_row))

    beta_draws <- do.call(cbind, lapply(posterior, function(s) {
      latent <- s$latent
      latent[cb_names]
    }))

    eci_draws <- as.numeric(t(beta_draws) %*% cb_row)

    # -----------------------
    # Bayesian: bdlnm
    # -----------------------
  } else if (inherits(fit, "bdlnm")) {

    beta_draws <- fit$coefficients

    if (is.null(dim(beta_draws))) {
      beta_draws <- matrix(beta_draws, ncol = 1)
    }

    cb_names <- intersect(rownames(beta_draws), paste0("cb_", var, "_", seq_along(cb_row)))

    if (length(cb_names) == 0) {
      cb_names <- intersect(cb_cols_fit %||% character(0), rownames(beta_draws))
    }

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between bdlnm coefficients and cb structure.")
    }

    beta_draws <- beta_draws[cb_names, , drop = FALSE]

    if (ncol(beta_draws) > n_samples) {
      set.seed(1)
      keep <- sample(seq_len(ncol(beta_draws)), n_samples)
      beta_draws <- beta_draws[, keep, drop = FALSE]
    }

    eci_draws <- as.numeric(t(beta_draws) %*% cb_row)

    # -----------------------
    # Frequentist → normal approximation
    # -----------------------
  } else {

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

      return(list(
        beta = stats::coef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    cv <- extract_coef_vcov(fit)

    beta_hat <- cv$beta
    V_hat <- cv$vcov

    cb_names <- grep(paste0("^cb_", var, "_"), names(beta_hat), value = TRUE)

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between coefficient vector and cb.")
    }

    beta_hat <- beta_hat[cb_names]
    V_hat <- V_hat[cb_names, cb_names, drop = FALSE]

    beta_draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat, Sigma = V_hat)

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
    ECI_weighted = mean(eci_draws),
    sd = stats::sd(eci_draws),
    lower = safe_quantile(eci_draws)[1],
    upper = safe_quantile(eci_draws)[2]
  ))
}
