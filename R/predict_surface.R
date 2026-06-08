#' Predict full DLNM exposure-lag-response surface
#' Computes the full DLNM exposure–lag–response surface using the fitted model,#'
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

  if (is.null(fit)) stop("`fit` cannot be NULL.")
  if (!is.data.frame(wx_long)) stop("`wx_long` must be a data.frame.")
  if (!all(c("epi_id", "dpp") %in% names(wx_long))) {
    stop("`wx_long` must contain 'epi_id' and 'dpp'.")
  }
  if (!var %in% names(wx_long)) {
    stop("`var` not found in `wx_long`.")
  }
  if (!is.logical(uncertainty)) stop("`uncertainty` must be TRUE/FALSE.")
  if (!is.numeric(n_samples) || n_samples <= 0) stop("`n_samples` must be positive.")
  n_samples <- as.integer(n_samples)

  fit_spec <- attr(fit, "epiexposure_spec")

  if (!is.null(fit_spec) && !is.null(fit_spec[[var]])) {
    spec_v <- fit_spec[[var]]
    lag_max_use <- as.integer(spec_v$lag_max)
    argvar <- spec_v$argvar
    arglag <- spec_v$arglag
  } else {
    lag_max_use <- lag_max

    argvar <- list(fun = fun_var %||% "ns", df = df_var %||% 4)
    arglag <- list(fun = fun_lag %||% "ns", df = df_lag %||% 4)
  }

  x_pool <- unlist(lapply(split(wx_long[[var]], wx_long$epi_id),
                          function(v) c(v, rep(NA, lag_max_use))))

  cb <- dlnm::crossbasis(x_pool, lag = lag_max_use, argvar = argvar, arglag = arglag)

  x_all <- wx_long[[var]]

  at_vals <- sort(unique(stats::quantile(x_all, probs = probs, na.rm = TRUE)))

  cen <- switch(
    ref$method,
    median = stats::median(x_all, na.rm = TRUE),
    percentile = stats::quantile(x_all, ref$value, na.rm = TRUE),
    fixed = ref$value
  )

  coef <- stats::coef(fit)
  vcov <- stats::vcov(fit)

  idx <- grepl(paste0("^cb_", var, "_"), names(coef))
  beta_sub <- coef[idx]
  V_sub <- vcov[idx, idx, drop = FALSE]

  base_cp <- dlnm::crosspred(cb, coef = beta_sub, vcov = V_sub,
                             at = at_vals, cen = cen, bylag = 1)

  if (!uncertainty) return(base_cp)

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' required.")
  }

  beta_draws <- MASS::mvrnorm(n = n_samples, mu = beta_sub, Sigma = V_sub)

  samples <- lapply(seq_len(n_samples), function(i) {
    dlnm::crosspred(cb,
                    coef = beta_draws[i, ],
                    vcov = NULL,
                    at = at_vals,
                    cen = cen,
                    bylag = 1)
  })

  if (output == "samples") return(samples)

  matfit_arr <- simplify2array(lapply(samples, function(x) x$matfit))

  center_fit <- apply(matfit_arr, c(1,2), stats::median, na.rm = TRUE)
  lower_fit  <- apply(matfit_arr, c(1,2), safe_quantile)[1,,]
  upper_fit  <- apply(matfit_arr, c(1,2), safe_quantile)[2,,]

  cp_mean  <- base_cp
  cp_lower <- base_cp
  cp_upper <- base_cp

  cp_mean$matfit  <- center_fit
  cp_lower$matfit <- lower_fit
  cp_upper$matfit <- upper_fit

  list(
    mean  = cp_mean,
    lower = cp_lower,
    upper = cp_upper
  )
}
