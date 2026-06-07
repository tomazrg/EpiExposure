#' Reduce DLNM effects to one dimension (article-consistent)
#'
#' @param fit Fitted model object
#' @param wx_long Long-format weather data
#' @param var Exposure variable (ex: "tmax")
#' @param lag_max Maximum lag
#' @param df_var Degrees of freedom (exposure)
#' @param df_lag Degrees of freedom (lag)
#' @param fun_var Basis ("ns","bs","poly","lin")
#' @param fun_lag Basis ("ns","ps","lin")
#' @param type "overall", "lag", "var"
#' @param value Required for type = "lag" or "var"
#' @param scale "link", "response", "percent"
#' @param uncertainty Logical
#' @param output "summary" or "samples"
#' @param n_samples Number of samples
#'
#' @return data.frame
#' @export
reduce_effects <- function(
    fit,
    wx_long,
    var,
    lag_max,
    df_var = 4,
    df_lag = 4,
    fun_var = "ns",
    fun_lag = "ns",
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
  # 1) pooled series
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

  x_pool <- build_pooled_series(wx_long, var, lag_max)

  # ----------------------------------------------------------
  # 2) crossbasis
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
  # 3) deterministic coef + vcov
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

  # ----------------------------------------------------------
  # 4) posterior/simulated draws by engine
  # returns matrix n_draws x p
  # ----------------------------------------------------------
  extract_beta_draws <- function(model, var, p, n_samples) {

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

  cv <- extract_coef_vcov(fit)
  beta <- cv$beta
  V    <- cv$vcov

  idx <- grepl(paste0("^cb_", var, "_"), names(beta))
  if (!any(idx)) stop("No DLNM terms found for variable: ", var)

  beta_sub <- beta[idx]
  V_sub    <- V[idx, idx, drop = FALSE]

  # ----------------------------------------------------------
  # NO UNCERTAINTY
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
      p = ncol(cb),
      n_samples = n_samples
    )

    n_draws <- nrow(beta_draws)
    res_list <- vector("list", n_draws)

    for (i in seq_len(n_draws)) {

      beta_i <- beta_draws[i, ]

      cr_i <- dlnm::crossreduce(
        cb,
        coef  = beta_i,
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
          eta = mean(eta, na.rm = TRUE),
          eta_sd = stats::sd(eta, na.rm = TRUE),
          low = safe_quantile(eta)[1],
          high = safe_quantile(eta)[2],
          .groups = "drop"
        )
    }
  }

  # ----------------------------------------------------------
  # scale transformation
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
