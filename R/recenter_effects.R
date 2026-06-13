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
#' @param lag_max Optional maximum lag. Ignored if `fit` contains `epiexposure_spec`.
#' @param df_var Optional degrees of freedom (exposure). Ignored if `fit`
#'   contains `epiexposure_spec`.
#' @param df_lag Optional degrees of freedom (lag). Ignored if `fit`
#'   contains `epiexposure_spec`.
#' @param fun_var Optional basis (`"ns"`, `"bs"`, `"poly"`, `"lin"`). Ignored
#'   if `fit` contains `epiexposure_spec`.
#' @param fun_lag Optional basis (`"ns"`, `"ps"`, `"lin"`). Ignored
#'   if `fit` contains `epiexposure_spec`.
#' @param ref New reference definition (centering value).
#' @param probs Quantiles used to define the exposure grid.
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty using
#'   simulated or posterior draws of the model coefficients.
#' @param output Character. `"summary"` returns aggregated surfaces;
#'   `"samples"` returns all simulated surfaces.
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
    lag_max = NULL,
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

  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(0)
    stats::sd(x)
  }

  # ------------------------------------------------------------
  # validations
  # ------------------------------------------------------------
  if (is.null(fit)) stop("`fit` cannot be NULL.")
  if (!is.data.frame(wx_long)) stop("`wx_long` must be a data.frame.")
  if (!all(c("epi_id", "dpp") %in% names(wx_long))) {
    stop("`wx_long` must contain at least 'epi_id' and 'dpp'.")
  }
  if (!is.character(var) || length(var) != 1L) {
    stop("`var` must be a single character string.")
  }
  if (!var %in% names(wx_long)) stop("`var` not found in `wx_long`.")
  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  fit_spec     <- attr(fit, "epiexposure_spec")
  cb_cols_fit  <- attr(fit, "epiexposure_cb_cols")
  dat_template <- attr(fit, "epiexposure_dat_template")

  # ------------------------------------------------------------
  # 1) Resolve basis spec (prefer fit metadata)
  # ------------------------------------------------------------
  if (!is.null(fit_spec) && !is.null(fit_spec[[var]])) {
    spec_v <- fit_spec[[var]]

    if (is.null(spec_v$lag_max)) {
      stop("Missing `lag_max` inside `epiexposure_spec` for variable: ", var)
    }

    # ✅ dimensão temporal
    lag_max_use <- as.integer(max(spec_v$lag_max))
    argvar <- spec_v$argvar
    arglag <- spec_v$arglag

  } else {

    if (is.null(lag_max)) {
      stop("Model does not contain `epiexposure_spec` for variable '", var,
           "'. Please provide `lag_max`.")
    }

    lag_max_use <- as.integer(max(lag_max))

    argvar <- switch(
      fun_var,
      ns   = list(fun = "ns", df = df_var),
      bs   = list(fun = "bs", df = df_var),
      poly = list(fun = "poly", degree = df_var),
      lin  = list(fun = "lin"),
      stop("Unsupported fun_var: ", fun_var)
    )

    if (!is.null(argvar$fun) && argvar$fun != "lin") {
      argvar$intercept <- FALSE
    }

    arglag <- switch(
      fun_lag,
      ns  = list(fun = "ns", df = df_lag),
      ps  = list(fun = "ps", df = df_lag),
      lin = list(fun = "lin"),
      stop("Unsupported fun_lag: ", fun_lag)
    )
  }

  # ------------------------------------------------------------
  # ✅ Check temporal coverage
  # ------------------------------------------------------------
  .check_lag_coverage <- function(dat, lag_max) {
    n_required <- lag_max + 1L

    bad_ids <- dat |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(n_days = dplyr::n_distinct(dpp), .groups = "drop") |>
      dplyr::filter(n_days < n_required)

    if (nrow(bad_ids) > 0) {
      stop(
        paste0(
          "Some epidemics do not have enough temporal coverage for lag_max.\n",
          "Required days per epi_id: ", n_required, "\n",
          "Example problematic epi_id: ",
          paste(head(bad_ids$epi_id, 5), collapse = ", ")
        )
      )
    }
  }

  .check_lag_coverage(wx_long, lag_max_use)

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
  # 5) Helpers: coefficient matching
  # ------------------------------------------------------------
  sort_cb_names <- function(x) {
    if (length(x) == 0) return(x)
    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    idx[is.na(idx)] <- seq_along(x)
    x[order(idx)]
  }

  match_brms_draw_names <- function(cb_names_ref, draw_colnames) {
    bnames <- paste0("b_", cb_names_ref)

    if (all(bnames %in% draw_colnames)) {
      return(bnames)
    }

    raw_match <- cb_names_ref[cb_names_ref %in% draw_colnames]
    if (length(raw_match) == length(cb_names_ref)) {
      return(raw_match)
    }

    stop("Could not match brms posterior draw names to crossbasis columns.")
  }

  # ------------------------------------------------------------
  # 6) Deterministic mean coef/vcov extractor
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
      V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
      if (is.null(V)) stop("Could not extract vcov from spaMM model.")
      return(list(
        beta = spaMM::fixef(model),
        vcov = V
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
  # 7) Select coefficient names robustly
  # ------------------------------------------------------------
  cv <- extract_coef_vcov(fit)
  beta <- cv$beta
  V    <- cv$vcov

  cb_cols_var <- NULL
  if (!is.null(cb_cols_fit)) {
    cb_cols_var <- grep(paste0("^cb_", var, "_"), cb_cols_fit, value = TRUE)
    cb_cols_var <- sort_cb_names(cb_cols_var)
  } else if (!is.null(dat_template)) {
    cb_cols_var <- grep(paste0("^cb_", var, "_"), names(dat_template), value = TRUE)
    cb_cols_var <- sort_cb_names(cb_cols_var)
  }

  cb_names_ref <- if (!is.null(cb_cols_var) && length(cb_cols_var) > 0) {
    cb_cols_var[cb_cols_var %in% names(beta)]
  } else {
    sort_cb_names(grep(paste0("^cb_", var, "_"), names(beta), value = TRUE))
  }

  if (length(cb_names_ref) == 0) {
    stop("No DLNM terms found for variable: ", var)
  }
  if (length(cb_names_ref) != p) {
    stop("Mismatch between selected coefficients and crossbasis dimension.")
  }

  beta_sub <- beta[cb_names_ref]
  V_sub    <- V[cb_names_ref, cb_names_ref, drop = FALSE]

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
  # 9) Posterior/simulated draws extractor
  # returns matrix n_draws x p
  # ------------------------------------------------------------
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

      draws <- do.call(rbind, lapply(posterior, function(s) {
        latent <- s$latent
        names(latent) <- gsub(":1$", "", names(latent))
        as.numeric(latent[cb_names_ref])
      }))

      colnames(draws) <- cb_names_ref
      return(draws)
    }

    # ---------- bdlnm: REAL posterior draws ----------
    if (inherits(model, "bdlnm")) {
      beta_draws <- model$coefficients

      if (is.null(dim(beta_draws))) {
        beta_draws <- matrix(beta_draws, ncol = 1)
      }

      if (!all(cb_names_ref %in% rownames(beta_draws))) {
        stop("Could not match bdlnm posterior draws to crossbasis columns.")
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
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty.")
    }

    draws <- MASS::mvrnorm(
      n = n_samples,
      mu = beta_sub,
      Sigma = V_sub
    )

    if (is.null(dim(draws))) {
      draws <- matrix(draws, nrow = 1)
    }

    colnames(draws) <- cb_names_ref
    return(draws)
  }

  beta_draws <- extract_beta_draws(
    model = fit,
    cb_names_ref = cb_names_ref,
    n_samples = n_samples
  )

  draw_sd <- apply(beta_draws, 2, stats::sd)
  if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
    warning(
      "Near-zero variability detected in coefficient draws for variable '", var,
      "'. Recentered surface intervals may collapse to a single value."
    )
  }

  # ------------------------------------------------------------
  # 10) Draw-based uncertainty
  # ------------------------------------------------------------
  n_draws <- nrow(beta_draws)
  cp_samples <- vector("list", n_draws)

  for (i in seq_len(n_draws)) {
    beta_i <- beta_draws[i, ]

    cp_samples[[i]] <- dlnm::crosspred(
      cb,
      coef  = beta_i,
      # ✅ API compatível com draws fixos
      vcov  = diag(0, length(beta_i)),
      at    = at_vals,
      cen   = cen,
      bylag = 1
    )
  }

  if (output == "samples") {
    # ainda assim retornamos mean/lower/upper para conveniência
  }

  # ------------------------------------------------------------
  # 11) Build summary surfaces
  # ------------------------------------------------------------
  n_at  <- length(at_vals)
  n_lag <- ncol(cp_samples[[1]]$matfit)

  matfit_arr <- array(NA_real_, dim = c(n_at, n_lag, n_draws))
  for (i in seq_len(n_draws)) {
    matfit_arr[, , i] <- cp_samples[[i]]$matfit
  }

  center_fit <- apply(matfit_arr, c(1, 2), stats::median, na.rm = TRUE)
  lower_fit  <- apply(matfit_arr, c(1, 2), function(z) safe_quantile(z)[1])
  upper_fit  <- apply(matfit_arr, c(1, 2), function(z) safe_quantile(z)[2])

  cp_mean  <- base_cp
  cp_lower <- base_cp
  cp_upper <- base_cp

  cp_mean$matfit   <- center_fit
  cp_lower$matfit  <- lower_fit
  cp_upper$matfit  <- upper_fit

  if (!is.null(base_cp$allfit)) {
    allfit_mat <- matrix(NA_real_, nrow = length(base_cp$allfit), ncol = n_draws)
    for (i in seq_len(n_draws)) {
      allfit_mat[, i] <- cp_samples[[i]]$allfit
    }

    cp_mean$allfit  <- apply(allfit_mat, 1, stats::median, na.rm = TRUE)
    cp_lower$allfit <- apply(allfit_mat, 1, function(z) safe_quantile(z)[1])
    cp_upper$allfit <- apply(allfit_mat, 1, function(z) safe_quantile(z)[2])
  }

  # manter campos RR coerentes, se existirem
  if (!is.null(base_cp$matRRfit)) {
    cp_mean$matRRfit  <- exp(cp_mean$matfit)
    cp_lower$matRRfit <- exp(cp_lower$matfit)
    cp_upper$matRRfit <- exp(cp_upper$matfit)
  }

  if (!is.null(base_cp$allRRfit) && !is.null(cp_mean$allfit)) {
    cp_mean$allRRfit  <- exp(cp_mean$allfit)
    cp_lower$allRRfit <- exp(cp_lower$allfit)
    cp_upper$allRRfit <- exp(cp_upper$allfit)
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
