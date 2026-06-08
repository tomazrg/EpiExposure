#' Recenter DLNM effects using a new reference value
#'
#' Recomputes the DLNM exposure-lag-response surface using a different
#' centering (reference) value, without refitting the model.
#'
#' When `uncertainty = TRUE`, the effect surface is recomputed across
#' simulated/posterior samples of the model coefficients. The central
#' estimate is obtained as the median of the simulated effects, while
#' uncertainty intervals are derived from empirical quantiles
#' (default: 2.5% and 97.5%).
#'
#' **Important:** although the output element is named `mean` for backward
#' compatibility, it represents the *central estimate*, computed as the
#' median when uncertainty is propagated.
#'
#' @param fit Fitted model from `fit_epidlnm()`.
#' @param wx_long Long-format weather data.
#' @param var Exposure variable (e.g. `"tmax"`).
#' @param lag_max Maximum lag.
#' @param df_var Degrees of freedom (exposure).
#' @param df_lag Degrees of freedom (lag).
#' @param fun_var Basis (`"ns"`, `"bs"`, `"poly"`, `"lin"`).
#' @param fun_lag Basis (`"ns"`, `"ps"`, `"lin"`).
#' @param ref New reference definition (centering value).
#' @param probs Quantiles used to define the exposure grid.
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty using
#' simulated or posterior draws of the model coefficients.
#' @param output Character. `"summary"` returns aggregated surfaces;
#' `"samples"` returns all simulated surfaces.
#' @param n_samples Integer. Number of samples used for uncertainty propagation.
#'
#' @return
#' - If `uncertainty = FALSE`: a `crosspred` object.
#'
#' - If `uncertainty = TRUE` and `output = "summary"`:
#'   a list with:
#'   - `mean`: central estimate surface (median-based)
#'   - `lower`: lower interval surface (quantile-based)
#'   - `upper`: upper interval surface (quantile-based)
#'
#' - If `uncertainty = TRUE` and `output = "samples"`:
#'   a list with:
#'   - `mean`: central estimate surface (median-based)
#'   - `lower`: lower interval surface
#'   - `upper`: upper interval surface
#'   - `samples`: list of `crosspred` objects for each simulation
#'
#' @details
#' Uncertainty is propagated using model-consistent sampling:
#' - Bayesian models (e.g., `brms`, `INLA`, `bdlnm`) use posterior draws
#' - Frequentist models use a normal approximation of the coefficient distribution
#'
#' The use of the median as the central estimate improves robustness to
#' asymmetry and non-normality in DLNM effect distributions, which commonly
#' arise from nonlinear exposure-lag-response relationships.
#'
#' @export
recenter_effects <- function(
    fit,
    wx_long,
    var,
    lag_max,
    df_var = 4,
    df_lag = 4,
    fun_var = "ns",
    fun_lag = "ns",
    ref = list(method = "median", value = NULL),
    probs = seq(0.05, 0.95, by = 0.01),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  output <- match.arg(output)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  if (!is.data.frame(wx_long)) stop("`wx_long` must be a data.frame.")
  if (!all(c("epi_id", "dpp") %in% names(wx_long))) {
    stop("`wx_long` must contain at least 'epi_id' and 'dpp'.")
  }
  if (!var %in% names(wx_long)) stop("`var` not found in `wx_long`.")
  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  fit_spec <- attr(fit, "epiexposure_spec")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")

  # ------------------------------------------------------------
  # 1) Resolve basis spec (prefer fit metadata)
  # ------------------------------------------------------------
  if (!is.null(fit_spec) && !is.null(fit_spec[[var]])) {
    spec_v <- fit_spec[[var]]
    lag_max_use <- as.integer(spec_v$lag_max)
    argvar <- spec_v$argvar
    arglag <- spec_v$arglag
  } else {
    lag_max_use <- lag_max

    argvar <- switch(
      fun_var,
      ns   = list(fun = "ns", df = df_var),
      bs   = list(fun = "bs", df = df_var),
      poly = list(fun = "poly", degree = df_var),
      lin  = list(fun = "lin")
    )

    if (!is.null(argvar$fun) && argvar$fun != "lin") {
      argvar$intercept <- FALSE
    }

    arglag <- switch(
      fun_lag,
      ns  = list(fun = "ns", df = df_lag),
      ps  = list(fun = "ps", df = df_lag),
      lin = list(fun = "lin")
    )
  }

  # ------------------------------------------------------------
  # 2) Helper: pooled series
  # ------------------------------------------------------------
  build_pooled_series <- function(dat, var, sep_n) {

    ids <- unique(dat$epi_id)
    out <- vector("list", length(ids))

    for (i in seq_along(ids)) {
      v <- dat |>
        dplyr::filter(epi_id == ids[i]) |>
        dplyr::arrange(dpp) |>
        dplyr::pull(.data[[var]])

      out[[i]] <- c(v, rep(NA_real_, sep_n))
    }

    unlist(out)
  }

  x_pool <- build_pooled_series(wx_long, var, lag_max_use)

  # ------------------------------------------------------------
  # 3) Rebuild crossbasis
  # ------------------------------------------------------------
  cb <- dlnm::crossbasis(
    x_pool,
    lag    = lag_max_use,
    argvar = argvar,
    arglag = arglag
  )

  p <- ncol(cb)

  # ------------------------------------------------------------
  # 4) Exposure grid and new center
  # ------------------------------------------------------------
  x_all <- wx_long[[var]]

  at_vals <- sort(unique(
    as.numeric(stats::quantile(x_all, probs = probs, na.rm = TRUE))
  ))

  cen <- switch(
    ref$method,
    median     = stats::median(x_all, na.rm = TRUE),
    percentile = stats::quantile(x_all, ref$value, na.rm = TRUE),
    fixed      = ref$value,
    stop("Invalid ref$method")
  )

  cen <- as.numeric(cen)

  # ------------------------------------------------------------
  # 5) Deterministic mean coef/vcov extractor
  # ------------------------------------------------------------
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
      return(list(
        beta = spaMM::fixef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
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

    return(list(
      beta = stats::coef(model),
      vcov = as.matrix(stats::vcov(model))
    ))
  }

  # ------------------------------------------------------------
  # 6) Posterior/simulated draws extractor
  # returns matrix n_draws x p
  # ------------------------------------------------------------
  extract_beta_draws <- function(model, var, p, n_samples, cb_cols_fit = NULL) {

    # ---------- brms: REAL posterior draws ----------
    if (inherits(model, "brmsfit")) {
      draws <- as.matrix(brms::as_draws_matrix(model))

      nm <- grep(paste0("^b_cb_", var, "_"), colnames(draws), value = TRUE)
      if (length(nm) == 0) {
        nm <- grep(paste0("^b_.*", var), colnames(draws), value = TRUE)
      }

      if (length(nm) != p) {
        stop("Could not match brms posterior draws to crossbasis columns.")
      }

      if (nrow(draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(nrow(draws)), n_samples)
        draws <- draws[keep, nm, drop = FALSE]
      } else {
        draws <- draws[, nm, drop = FALSE]
      }

      return(draws)
    }

    # ---------- INLA: REAL posterior draws ----------
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

    # ---------- bdlnm: REAL posterior draws ----------
    if (inherits(model, "bdlnm")) {
      beta_draws <- model$coefficients
      if (is.null(dim(beta_draws))) {
        beta_draws <- matrix(beta_draws, ncol = 1)
      }

      cb_names <- intersect(
        rownames(beta_draws),
        paste0("cb_", var, "_", seq_len(p))
      )

      if (length(cb_names) == 0) {
        cb_names <- grep(paste0("^cb_", var, "_"), rownames(beta_draws), value = TRUE)
      }
      if (length(cb_names) == 0 && !is.null(cb_cols_fit)) {
        cb_names <- intersect(cb_cols_fit, rownames(beta_draws))
      }

      if (length(cb_names) != p) {
        stop("Could not match bdlnm posterior draws to crossbasis columns.")
      }

      if (ncol(beta_draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(ncol(beta_draws)), n_samples)
        beta_draws <- beta_draws[cb_names, keep, drop = FALSE]
      } else {
        beta_draws <- beta_draws[cb_names, , drop = FALSE]
      }

      return(t(beta_draws))
    }

    # ---------- Frequentist: normal approximation ----------
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty.")
    }

    cv <- extract_coef_vcov(model)
    beta <- cv$beta
    V    <- cv$vcov

    idx <- grepl(paste0("^cb_", var, "_"), names(beta))
    beta_sub <- beta[idx]
    V_sub    <- V[idx, idx, drop = FALSE]

    if (length(beta_sub) != p) {
      stop("Mismatch between coefficient vector and crossbasis columns.")
    }

    draws <- MASS::mvrnorm(
      n = n_samples,
      mu = beta_sub,
      Sigma = V_sub
    )

    if (is.null(dim(draws))) {
      draws <- matrix(draws, nrow = 1)
    }

    return(draws)
  }

  # ------------------------------------------------------------
  # 7) Deterministic coefficients for the base return
  # ------------------------------------------------------------
  cv <- extract_coef_vcov(fit)
  beta <- cv$beta
  V    <- cv$vcov

  idx <- grepl(paste0("^cb_", var, "_"), names(beta))
  beta_sub <- beta[idx]
  V_sub    <- V[idx, idx, drop = FALSE]

  if (length(beta_sub) == 0) {
    stop("No DLNM terms found for variable: ", var)
  }
  if (length(beta_sub) != p) {
    stop("Mismatch between selected coefficients and crossbasis dimension.")
  }

  # ------------------------------------------------------------
  # 8) Deterministic crosspred (backward compatible)
  # ------------------------------------------------------------
  base_cp <- dlnm::crosspred(
    cb,
    coef  = beta_sub,
    vcov  = V_sub,
    at    = at_vals,
    cen   = cen,
    bylag = 1
  )

  if (!uncertainty) {
    return(base_cp)
  }

  # ------------------------------------------------------------
  # 9) Draw-based uncertainty
  # ------------------------------------------------------------
  beta_draws <- extract_beta_draws(
    model = fit,
    var = var,
    p = p,
    n_samples = n_samples,
    cb_cols_fit = cb_cols_fit
  )

  n_draws <- nrow(beta_draws)
  cp_samples <- vector("list", n_draws)

  for (i in seq_len(n_draws)) {
    cp_samples[[i]] <- dlnm::crosspred(
      cb,
      coef  = beta_draws[i, ],
      vcov  = NULL,
      at    = at_vals,
      cen   = cen,
      bylag = 1
    )
  }

  # ------------------------------------------------------------
  # 10) Build summary surfaces
  # ------------------------------------------------------------
  # Extract arrays of matfit and allfit from samples
  matfit_arr <- simplify2array(lapply(cp_samples, function(x) x$matfit))
  allfit_arr <- simplify2array(lapply(cp_samples, function(x) x$allfit))

  # matfit: [at x lag x sample]
  matfit_mean  <- apply(matfit_arr, c(1, 2), stats::median, na.rm = TRUE)
  matfit_lower <- apply(matfit_arr, c(1, 2), stats::quantile, probs = 0.025, na.rm = TRUE)
  matfit_upper <- apply(matfit_arr, c(1, 2), stats::quantile, probs = 0.975, na.rm = TRUE)

  # allfit usually [at x sample]
  if (length(dim(allfit_arr)) == 2) {
    allfit_mean  <- apply(allfit_arr, 1, stats::median, na.rm = TRUE)
    allfit_lower <- apply(allfit_arr, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
    allfit_upper <- apply(allfit_arr, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
  } else {
    allfit_mean  <- base_cp$allfit
    allfit_lower <- base_cp$allfit
    allfit_upper <- base_cp$allfit
  }

  cp_mean  <- base_cp
  cp_lower <- base_cp
  cp_upper <- base_cp

  cp_mean$matfit  <- matfit_mean
  cp_lower$matfit <- matfit_lower
  cp_upper$matfit <- matfit_upper

  if (!is.null(base_cp$allfit)) {
    cp_mean$allfit  <- allfit_mean
    cp_lower$allfit <- allfit_lower
    cp_upper$allfit <- allfit_upper
  }

  # try to keep RR fields coherent if present
  if (!is.null(base_cp$matRRfit)) {
    cp_mean$matRRfit  <- exp(matfit_mean)
    cp_lower$matRRfit <- exp(matfit_lower)
    cp_upper$matRRfit <- exp(matfit_upper)
  }
  if (!is.null(base_cp$allRRfit)) {
    cp_mean$allRRfit  <- exp(allfit_mean)
    cp_lower$allRRfit <- exp(allfit_lower)
    cp_upper$allRRfit <- exp(allfit_upper)
  }

  if (output == "samples") {
    return(list(
      mean = cp_mean,
      lower = cp_lower,
      upper = cp_upper,
      samples = cp_samples
    ))
  }

  return(list(
    mean = cp_mean,
    lower = cp_lower,
    upper = cp_upper
  ))
}
