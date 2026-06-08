#' Reduce DLNM effects to one dimension (article-consistent)
#'
#' Reduces a DLNM exposure-lag-response surface into a one-dimensional summary
#' using `dlnm::crossreduce()`, with optional uncertainty propagation.
#'
#' When `uncertainty = TRUE`, summaries are computed from simulated/posterior
#' samples of the model coefficients. The central estimate is obtained as the
#' median of the simulated effects, while uncertainty intervals are derived
#' from empirical quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** although the output element is named `mean` for backward
#' compatibility, it represents the *central estimate*, computed as the median
#' when uncertainty is propagated.
#'
#' @param fit Fitted model object.
#' @param wx_long Long-format weather data.
#' @param var Exposure variable (e.g. `"tmax"`).
#' @param lag_max Optional maximum lag. Ignored if `fit` contains `epiexposure_spec`.
#' @param df_var Optional degrees of freedom (exposure). Ignored if `fit` contains `epiexposure_spec`.
#' @param df_lag Optional degrees of freedom (lag). Ignored if `fit` contains `epiexposure_spec`.
#' @param fun_var Optional basis (`"ns"`, `"bs"`, `"poly"`, `"lin"`). Ignored if `fit` contains `epiexposure_spec`.
#' @param fun_lag Optional basis (`"ns"`, `"ps"`, `"lin"`). Ignored if `fit` contains `epiexposure_spec`.
#' @param type Reduction type: `"overall"`, `"lag"`, or `"var"`.
#' @param value Required for `type = "lag"` or `type = "var"`.
#' @param scale Output scale: `"link"`, `"response"`, or `"percent"`.
#' @param uncertainty Logical. If `TRUE`, propagate uncertainty using simulated
#' or posterior draws of the model coefficients.
#' @param output Character. `"summary"` returns aggregated estimates; `"samples"`
#' returns all simulated values.
#' @param n_samples Integer. Number of samples used for uncertainty propagation.
#'
#' @return A data.frame containing reduced effects. When `uncertainty = TRUE`
#' and `output = "summary"`, the result includes:
#'   - central estimate (median; stored in `eta`)
#'   - `eta_sd`: standard deviation of simulated values
#'   - `low` / `high`: empirical interval limits (quantiles)
#'
#' @details
#' Uncertainty is propagated using model-consistent sampling:
#' - Bayesian models (e.g., `brms`, `INLA`, `bdlnm`) use posterior draws
#' - Frequentist models use a normal approximation of the coefficient distribution
#'
#' The use of the median as the central estimate improves robustness under
#' non-normal or asymmetric effect distributions, which commonly arise in
#' DLNM applications.
#'
#' @export
reduce_effects <- function(
    fit,
    wx_long,
    var,
    lag_max = NULL,
    df_var = NULL,
    df_lag = NULL,
    fun_var = NULL,
    fun_lag = NULL,
    type = c("overall", "lag", "var"),
    value = NULL,
    scale = c("percent", "response", "link"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  type  <- match.arg(type)
  scale <- match.arg(scale)
  output <- match.arg(output)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  # ----------------------------------------------------------
  # basic validation
  # ----------------------------------------------------------
  if (is.null(fit)) stop("`fit` cannot be NULL.")
  if (!is.data.frame(wx_long)) stop("`wx_long` must be a data.frame.")
  if (!all(c("epi_id", "dpp") %in% names(wx_long))) {
    stop("`wx_long` must contain at least 'epi_id' and 'dpp'.")
  }
  if (!var %in% names(wx_long)) {
    stop("`var` not found in `wx_long`.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  fit_spec <- attr(fit, "epiexposure_spec")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")

  # ----------------------------------------------------------
  # 1) Resolve basis spec (PRIORITIZE FIT)
  # ----------------------------------------------------------
  if (!is.null(fit_spec) && !is.null(fit_spec[[var]])) {

    spec_v <- fit_spec[[var]]
    lag_max_use <- as.integer(spec_v$lag_max)
    argvar <- spec_v$argvar
    arglag <- spec_v$arglag

  } else {

    if (is.null(lag_max)) {
      stop("Model does not contain `epiexposure_spec` for variable '", var,
           "'. Please provide `lag_max`.")
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
  # 2) pooled series
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

  # ----------------------------------------------------------
  # 3) crossbasis
  # ----------------------------------------------------------
  cb <- dlnm::crossbasis(
    x_pool,
    lag    = lag_max_use,
    argvar = argvar,
    arglag = arglag
  )

  p <- ncol(cb)

  cb_cols_var <- NULL
  if (!is.null(cb_cols_fit)) {
    cb_cols_var <- grep(paste0("^cb_", var, "_"), cb_cols_fit, value = TRUE)
  }

  # ----------------------------------------------------------
  # 4) coef/vcov extraction
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

  get_cb_coef_names <- function(coef_names, var, p, cb_cols_var = NULL) {

    nm <- if (!is.null(cb_cols_var)) {
      intersect(cb_cols_var, coef_names)
    } else {
      grep(paste0("^cb_", var, "_"), coef_names, value = TRUE)
    }

    if (length(nm) == 0) {
      nm <- grep(paste0("^cb_", var, "_"), coef_names, value = TRUE)
    }

    if (length(nm) != p) {
      stop("Could not match coefficient names to crossbasis columns.")
    }

    nm
  }

  # ----------------------------------------------------------
  # 5) posterior/simulated draws by engine
  # returns matrix n_draws x p
  # ----------------------------------------------------------
  extract_beta_draws <- function(model, var, p, n_samples, cb_cols_var = NULL) {

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

      cb_names <- cb_cols_var %||% paste0("cb_", var, "_", seq_len(p))

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

      cb_names <- intersect(rownames(beta_draws), cb_cols_var %||% character(0))
      if (length(cb_names) == 0) {
        cb_names <- intersect(rownames(beta_draws), paste0("cb_", var, "_", seq_len(p)))
      }
      if (length(cb_names) == 0) {
        cb_names <- grep(paste0("^cb_", var, "_"), rownames(beta_draws), value = TRUE)
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
      stop("Package 'MASS' required for frequentist uncertainty.")
    }

    cv <- extract_coef_vcov(model)
    beta <- cv$beta
    V    <- cv$vcov

    nm <- get_cb_coef_names(names(beta), var, p, cb_cols_var)

    beta_sub <- beta[nm]
    V_sub    <- V[nm, nm, drop = FALSE]

    draws <- MASS::mvrnorm(
      n = n_samples,
      mu = beta_sub,
      Sigma = V_sub
    )

    if (is.null(dim(draws))) {
      draws <- matrix(draws, nrow = 1)
    }

    colnames(draws) <- nm
    return(draws)
  }

  cv <- extract_coef_vcov(fit)
  beta <- cv$beta
  V    <- cv$vcov

  nm <- get_cb_coef_names(names(beta), var, p, cb_cols_var)

  beta_sub <- beta[nm]
  V_sub    <- V[nm, nm, drop = FALSE]

  # ----------------------------------------------------------
  # 6) NO UNCERTAINTY
  # ----------------------------------------------------------
  if (!uncertainty) {

    cr <- dlnm::crossreduce(
      cb,
      coef  = beta_sub,
      vcov  = V_sub,
      type  = type,
      value = value
    )

    df <- data.frame(
      x   = cr$predvar,
      eta = cr$fit,
      low = cr$low,
      high = cr$high
    )

  } else {

    beta_draws <- extract_beta_draws(
      model = fit,
      var = var,
      p = p,
      n_samples = n_samples,
      cb_cols_var = cb_cols_var
    )

    n_draws <- nrow(beta_draws)
    res_list <- vector("list", n_draws)

    for (i in seq_len(n_draws)) {

      cr_i <- dlnm::crossreduce(
        cb,
        coef  = beta_draws[i, ],
        vcov  = NULL,
        type  = type,
        value = value
      )

      res_list[[i]] <- data.frame(
        x = cr_i$predvar,
        eta = cr_i$fit,
        sample = i
      )
    }

    all_draws <- do.call(rbind, res_list)

    if (output == "samples") {
      df <- all_draws
    } else {
      df <- all_draws |>
        dplyr::group_by(x) |>
        dplyr::summarise(
          eta = stats::median(eta, na.rm = TRUE),
          eta_sd = stats::sd(eta, na.rm = TRUE),
          low = safe_quantile(eta)[1],
          high = safe_quantile(eta)[2],
          .groups = "drop"
        )
    }
  }

  # ----------------------------------------------------------
  # 7) scale transformation
  # ----------------------------------------------------------
  if (scale == "link") {

    df$effect   <- df$eta
    df$low_eff  <- if ("low" %in% names(df)) df$low else NA_real_
    df$high_eff <- if ("high" %in% names(df)) df$high else NA_real_

  } else if (scale == "response") {

    df$effect   <- exp(df$eta)
    df$low_eff  <- if ("low" %in% names(df)) exp(df$low) else NA_real_
    df$high_eff <- if ("high" %in% names(df)) exp(df$high) else NA_real_

  } else if (scale == "percent") {

    df$effect   <- (exp(df$eta) - 1) * 100
    df$low_eff  <- if ("low" %in% names(df)) (exp(df$low) - 1) * 100 else NA_real_
    df$high_eff <- if ("high" %in% names(df)) (exp(df$high) - 1) * 100 else NA_real_
  }

  df$type  <- type
  df$value <- value

  return(df)
}
