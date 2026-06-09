#' Predict full DLNM exposure-lag-response surface
#'
#' Computes the full DLNM exposure-lag-response surface using the fitted model,
#' without refitting. The surface is evaluated over a grid of exposure values
#' and lags defined by the model specification or user input.
#'
#' This function supports both deterministic predictions and uncertainty
#' propagation. When `uncertainty = TRUE`, the surface is recomputed across
#' simulated or posterior draws of the model coefficients.
#'
#' If `output = "summary"`, the central surface is computed as the median
#' of the simulated surfaces, while interval limits are obtained from
#' empirical quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** although the returned object uses the name `mean` for backward
#' compatibility, it represents the *central estimate*, computed as the median
#' when uncertainty is propagated.
#'
#' @param fit Fitted model object.
#' @param wx_long Long-format weather data.
#' @param var Exposure variable name.
#' @param lag_max Optional maximum lag. Ignored if fit contains `epiexposure_spec`.
#' @param df_var Optional degrees of freedom (exposure). Ignored if fit contains metadata.
#' @param df_lag Optional degrees of freedom (lag). Ignored if fit contains metadata.
#' @param fun_var Optional basis function for exposure.
#' @param fun_lag Optional basis function for lag.
#' @param ref Reference exposure definition.
#' @param probs Quantiles used for the exposure grid.
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Number of simulations.
#'
#' @return
#' - If `uncertainty = FALSE`: a `crosspred` object.
#'
#' - If `uncertainty = TRUE` and `output = "summary"`:
#'   a list with:
#'   - `mean`: central surface (median-based)
#'   - `lower`: lower surface (quantile-based)
#'   - `upper`: upper surface (quantile-based)
#'
#' - If `uncertainty = TRUE` and `output = "samples"`:
#'   a list of `crosspred` objects (one per simulation).
#'
#' @details
#' Uncertainty is propagated using model-consistent sampling:
#' - Bayesian models (e.g., `brms`, `INLA`, `bdlnm`) use posterior draws
#' - Frequentist models use simulation from the asymptotic coefficient distribution
#'
#' Using the median as the central estimate improves robustness to asymmetric
#' distributions commonly observed in DLNM surfaces.
#'
#' @export
predict_surface <- function(
    fit,
    wx_long,
    var,
    lag_max = NULL,
    df_var = NULL,
    df_lag = NULL,
    fun_var = NULL,
    fun_lag = NULL,
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

  # ----------------------------------------------------------
  # VALIDATIONS
  # ----------------------------------------------------------
  if (is.null(fit)) stop("`fit` cannot be NULL.")
  if (!is.data.frame(wx_long)) stop("`wx_long` must be a data.frame.")
  if (!all(c("epi_id", "dpp") %in% names(wx_long))) {
    stop("`wx_long` must contain 'epi_id' and 'dpp'.")
  }
  if (!is.character(var) || length(var) != 1L) {
    stop("`var` must be a single character string.")
  }
  if (!var %in% names(wx_long)) {
    stop("`var` not found in `wx_long`.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE/FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  fit_spec     <- attr(fit, "epiexposure_spec")
  cb_cols_fit  <- attr(fit, "epiexposure_cb_cols")
  dat_template <- attr(fit, "epiexposure_dat_template")

  # ----------------------------------------------------------
  # RESOLVE BASIS SPEC
  # ----------------------------------------------------------
  if (!is.null(fit_spec) && !is.null(fit_spec[[var]])) {

    spec_v <- fit_spec[[var]]

    if (is.null(spec_v$lag_max)) {
      stop("Missing `lag_max` inside `epiexposure_spec` for variable: ", var)
    }

    lag_max_use <- as.integer(max(spec_v$lag_max))
    argvar <- spec_v$argvar
    arglag <- spec_v$arglag

  } else {

    if (is.null(lag_max)) {
      stop("`lag_max` must be provided when fit does not contain `epiexposure_spec`.")
    }

    lag_max_use <- as.integer(lag_max)

    fun_var_use <- fun_var %||% "ns"
    fun_lag_use <- fun_lag %||% "ns"
    df_var_use  <- df_var %||% 4
    df_lag_use  <- df_lag %||% 4

    argvar <- switch(
      fun_var_use,
      ns   = list(fun = "ns", df = df_var_use),
      bs   = list(fun = "bs", df = df_var_use),
      poly = list(fun = "poly", degree = df_var_use),
      lin  = list(fun = "lin"),
      stop("Unsupported fun_var: ", fun_var_use)
    )

    if (!is.null(argvar$fun) && argvar$fun != "lin") {
      argvar$intercept <- FALSE
    }

    arglag <- switch(
      fun_lag_use,
      ns  = list(fun = "ns", df = df_lag_use),
      ps  = list(fun = "ps", df = df_lag_use),
      lin = list(fun = "lin"),
      stop("Unsupported fun_lag: ", fun_lag_use)
    )
  }

  # ----------------------------------------------------------
  # BUILD POOLED SERIES
  # ----------------------------------------------------------
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

  cb <- dlnm::crossbasis(
    x_pool,
    lag    = lag_max_use,
    argvar = argvar,
    arglag = arglag
  )

  # ----------------------------------------------------------
  # EXPOSURE GRID / REFERENCE
  # ----------------------------------------------------------
  x_all <- wx_long[[var]]

  at_vals <- sort(unique(
    as.numeric(stats::quantile(x_all, probs = probs, na.rm = TRUE))
  ))

  cen <- switch(
    ref$method,
    median     = stats::median(x_all, na.rm = TRUE),
    percentile = stats::quantile(x_all, ref$value, na.rm = TRUE),
    fixed      = ref$value,
    stop("Invalid ref$method. Use 'median', 'percentile', or 'fixed'.")
  )
  cen <- as.numeric(cen)

  # ----------------------------------------------------------
  # IDENTIFY CB COLUMN NAMES FOR THIS VARIABLE
  # ----------------------------------------------------------
  cb_cols_var <- NULL

  if (!is.null(cb_cols_fit)) {
    cb_cols_var <- grep(paste0("^cb_", var, "_"), cb_cols_fit, value = TRUE)
  } else if (!is.null(dat_template)) {
    cb_cols_var <- grep(paste0("^cb_", var, "_"), names(dat_template), value = TRUE)
  }

  # fallback by expected names
  if (is.null(cb_cols_var) || length(cb_cols_var) == 0) {
    cb_cols_var <- paste0("cb_", var, "_", seq_len(ncol(cb)))
  }

  # ----------------------------------------------------------
  # EXTRACT COEF / VCOV BY ENGINE
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
      fe <- brms::fixef(model)
      beta <- fe[, "Estimate"]
      names(beta) <- rownames(fe)
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
      if (is.null(model$coefficients.summary) || is.null(model$coefficients)) {
        stop("bdlnm object does not contain coefficient summaries / posterior samples.")
      }
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
  # EXTRACT DRAWS BY ENGINE
  # ----------------------------------------------------------
  extract_beta_draws <- function(model, n_samples, var, cb_ncol, cb_cols_var = NULL) {

    # ---------- brms: REAL posterior draws ----------
    if (inherits(model, "brmsfit")) {

      draws <- as.matrix(brms::as_draws_matrix(model))

      nm <- if (!is.null(cb_cols_var)) {
        paste0("b_", cb_cols_var)
      } else {
        grep(paste0("^b_cb_", var, "_"), colnames(draws), value = TRUE)
      }

      if (length(nm) == 0) {
        nm <- grep(paste0("^b_cb_", var, "_"), colnames(draws), value = TRUE)
      }
      if (length(nm) == 0) {
        nm <- grep(paste0("^b_.*", var), colnames(draws), value = TRUE)
      }

      if (length(nm) != cb_ncol) {
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
        stop("Package 'INLA' is required for INLA uncertainty.")
      }

      posterior <- tryCatch(
        INLA::inla.posterior.sample(n = n_samples, result = model),
        error = function(e) NULL
      )

      if (is.null(posterior)) {
        stop("Could not draw INLA posterior samples. Ensure model was fitted with config = TRUE.")
      }

      cb_names <- cb_cols_var %||% paste0("cb_", var, "_", seq_len(cb_ncol))

      draws <- do.call(rbind, lapply(posterior, function(s) {
        latent <- s$latent
        names(latent) <- gsub(":1$", "", names(latent))
        vals <- latent[cb_names]
        if (any(is.na(vals))) {
          stop("Could not match INLA posterior draw names to crossbasis columns.")
        }
        as.numeric(vals)
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

      cb_names <- intersect(rownames(beta_draws), cb_cols_var %||% character(0))
      if (length(cb_names) == 0) {
        cb_names <- intersect(rownames(beta_draws), paste0("cb_", var, "_", seq_len(cb_ncol)))
      }
      if (length(cb_names) == 0) {
        cb_names <- grep(paste0("^cb_", var, "_"), rownames(beta_draws), value = TRUE)
      }

      if (length(cb_names) != cb_ncol) {
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
    cv <- extract_coef_vcov(model)
    beta_hat <- cv$beta
    V_hat <- cv$vcov

    cb_names <- if (!is.null(cb_cols_var)) {
      intersect(cb_cols_var, names(beta_hat))
    } else {
      grep(paste0("^cb_", var, "_"), names(beta_hat), value = TRUE)
    }

    if (length(cb_names) != cb_ncol) {
      stop("Could not match coefficient names to crossbasis structure.")
    }

    beta_hat <- beta_hat[cb_names]
    V_hat <- V_hat[cb_names, cb_names, drop = FALSE]

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty approximation.")
    }

    draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat, Sigma = V_hat)

    if (is.null(dim(draws))) {
      draws <- matrix(draws, nrow = 1)
    }

    colnames(draws) <- cb_names
    draws
  }

  # ----------------------------------------------------------
  # DETERMINISTIC: EXTRACT BETA / VCOV
  # ----------------------------------------------------------
  cv <- extract_coef_vcov(fit)
  beta_full <- cv$beta
  vcov_full <- cv$vcov

  if (is.null(names(beta_full))) {
    stop("Could not determine coefficient names from fitted model.")
  }
  if (!is.matrix(vcov_full)) {
    stop("Could not extract a valid covariance matrix from fitted model.")
  }

  cb_names <- if (!is.null(cb_cols_var)) {
    intersect(cb_cols_var, names(beta_full))
  } else {
    grep(paste0("^cb_", var, "_"), names(beta_full), value = TRUE)
  }

  if (length(cb_names) == 0) {
    cb_names <- grep(paste0("^cb_", var, "_"), names(beta_full), value = TRUE)
  }

  if (length(cb_names) != ncol(cb)) {
    stop(
      "Could not match coefficient names to the crossbasis structure for variable '", var, "'. ",
      "Expected ", ncol(cb), " coefficients, found ", length(cb_names), "."
    )
  }

  beta_sub <- beta_full[cb_names]
  V_sub <- vcov_full[cb_names, cb_names, drop = FALSE]

  if (!(nrow(V_sub) == length(beta_sub) && ncol(V_sub) == length(beta_sub))) {
    stop("Subset covariance matrix does not match coefficient dimension.")
  }

  # ----------------------------------------------------------
  # BASE CROSSPRED
  # ----------------------------------------------------------
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

  # ----------------------------------------------------------
  # UNCERTAINTY
  # ----------------------------------------------------------
  beta_draws <- extract_beta_draws(
    model = fit,
    n_samples = n_samples,
    var = var,
    cb_ncol = ncol(cb),
    cb_cols_var = cb_cols_var
  )

  samples <- lapply(seq_len(n_samples), function(i) {

    beta_i <- beta_draws[i, ]

    dlnm::crosspred(
      cb,
      coef  = beta_i,
      vcov  = diag(0, length(beta_i)),
      at    = at_vals,
      cen   = cen,
      bylag = 1
    )
  })

  if (output == "samples") {
    return(samples)
  }

  # ----------------------------------------------------------
  # SUMMARY OF SURFACES
  # ----------------------------------------------------------
  mats <- lapply(samples, function(x) x$matfit)

  n_r <- nrow(mats[[1]])
  n_c <- ncol(mats[[1]])
  n_s <- length(mats)

  matfit_arr <- array(NA_real_, dim = c(n_r, n_c, n_s))
  for (i in seq_len(n_s)) {
    matfit_arr[, , i] <- mats[[i]]
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

  list(
    mean  = cp_mean,
    lower = cp_lower,
    upper = cp_upper
  )
}
