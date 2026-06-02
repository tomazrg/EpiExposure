#' Predict full DLNM exposure-lag-response surface
#'
#' @param fit Fitted model object
#' @param wx_long Long-format weather data
#' @param var Exposure variable name
#' @param lag_max Maximum lag
#' @param df_var Degrees of freedom (exposure)
#' @param df_lag Degrees of freedom (lag)
#' @param fun_var Basis function
#' @param fun_lag Basis function
#' @param ref Reference exposure definition
#' @param probs Quantiles for exposure grid
#' @param uncertainty Logical
#' @param output "summary" or "samples"
#' @param n_samples Number of simulations
#'
#' @return crosspred object OR list
#' @export
predict_surface <- function(
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

  # ----------------------------------------------------------
  # pooled series
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

  SEPARATOR <- lag_max
  x_pool <- build_pooled_series(wx_long, var, SEPARATOR)

  # ----------------------------------------------------------
  # crossbasis
  # ----------------------------------------------------------
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

  cb <- dlnm::crossbasis(
    x_pool,
    lag    = lag_max,
    argvar = argvar,
    arglag = arglag
  )

  # ----------------------------------------------------------
  # extract beta / vcov
  # ----------------------------------------------------------
  extract_coef_vcov <- function(model) {

    if (inherits(model, "glmmTMB")) {
      return(list(
        beta = glmmTMB::fixef(model)$cond,
        vcov = vcov(model)$cond
      ))
    }

    if (inherits(model, "brmsfit")) {
      fe <- brms::fixef(model)
      beta <- fe[, "Estimate"]
      V <- as.matrix(stats::vcov(model))
      names(beta) <- rownames(fe)
      return(list(beta = beta, vcov = V))
    }

    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      V <- diag(model$summary.fixed$sd^2)
      rownames(V) <- names(beta)
      colnames(V) <- names(beta)
      return(list(beta = beta, vcov = V))
    }

    if (inherits(model, "bdlnm")) {
      return(list(
        beta = model$coefficients.summary[, "mean"],
        vcov = stats::cov(t(model$coefficients))
      ))
    }

    return(list(
      beta = coef(model),
      vcov = vcov(model)
    ))
  }

  bv <- extract_coef_vcov(fit)
  beta <- bv$beta
  V    <- bv$vcov

  idx <- grepl(paste0("^cb_", var, "_"), names(beta))
  beta_sub <- beta[idx]
  V_sub    <- V[idx, idx, drop = FALSE]

  # ----------------------------------------------------------
  # grid
  # ----------------------------------------------------------
  x_all <- wx_long[[var]]

  at_vals <- sort(unique(
    as.numeric(quantile(x_all, probs = probs, na.rm = TRUE))
  ))

  cen <- switch(
    ref$method,
    median     = median(x_all, na.rm = TRUE),
    percentile = quantile(x_all, ref$value, na.rm = TRUE),
    fixed      = ref$value,
    stop("Invalid ref$method")
  )

  cen <- as.numeric(cen)

  # ----------------------------------------------------------
  # BASE (no uncertainty)
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

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' required for uncertainty.")
  }

  beta_draws <- MASS::mvrnorm(
    n = n_samples,
    mu = beta_sub,
    Sigma = V_sub
  )

  samples <- vector("list", n_samples)

  for (i in seq_len(n_samples)) {
    samples[[i]] <- dlnm::crosspred(
      cb,
      coef  = beta_draws[i, ],
      vcov  = NULL,
      at    = at_vals,
      cen   = cen,
      bylag = 1
    )
  }

  if (output == "samples") {
    return(samples)
  }

  # ----------------------------------------------------------
  # SUMMARY
  # ----------------------------------------------------------

  get_array <- function(samples, slot) {
    simplify2array(lapply(samples, function(x) x[[slot]]))
  }

  matfit_arr <- get_array(samples, "matfit")

  mean_fit <- apply(matfit_arr, c(1, 2), mean)
  lower_fit <- apply(matfit_arr, c(1, 2), quantile, probs = 0.025)
  upper_fit <- apply(matfit_arr, c(1, 2), quantile, probs = 0.975)

  cp_mean  <- base_cp
  cp_lower <- base_cp
  cp_upper <- base_cp

  cp_mean$matfit  <- mean_fit
  cp_lower$matfit <- lower_fit
  cp_upper$matfit <- upper_fit

  return(list(
    mean  = cp_mean,
    lower = cp_lower,
    upper = cp_upper
  ))
}
