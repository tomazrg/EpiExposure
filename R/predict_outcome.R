#' Predict outcome under user-defined exposure-lag profile(s)
#'
#' This function predicts the outcome (single value per scenario/ID) given
#' exposure-lag profile(s). It retrieves DLNM specification and template data
#' from the fitted model (via attributes).
#'
#' @param fit Fitted model returned by fit_epidlnm()
#' @param profiles Numeric vector (single exposure) OR named list of numeric vectors
#' @param re "population" (no random effect) or "conditional" (include random effect)
#' @param id Optional vector of grouping levels (e.g., epi_id)
#' @param allow_new_levels Allow unseen grouping levels (for mixed models)
#' @param type Prediction scale ("response","link","conditional")
#' @param uncertainty Logical; if TRUE, quantify uncertainty
#' @param output "summary" or "samples"
#' @param n_samples Number of samples used for uncertainty quantification
#'
#' @return data.frame
#'   - If uncertainty = FALSE: same as before, with column `prediction`
#'   - If uncertainty = TRUE and output = "summary": `prediction`, `sd`, `lower`, `upper`
#'   - If uncertainty = TRUE and output = "samples": one row per sample (and id if provided)
#' @export
predict_outcome <- function(
    fit,
    profiles,
    re = c("population", "conditional"),
    id = NULL,
    allow_new_levels = FALSE,
    type = c("response", "link", "conditional"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  re <- match.arg(re)
  type <- match.arg(type)
  output <- match.arg(output)

  # -----------------------
  # ✅ VALIDATIONS
  # -----------------------
  if (is.null(fit)) stop("`fit` cannot be NULL.")

  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }

  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  spec <- attr(fit, "epiexposure_spec")
  nd0  <- attr(fit, "epiexposure_dat_template")
  vars_fit <- attr(fit, "epiexposure_vars")
  id_col <- attr(fit, "epiexposure_id_col")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")
  family_fit <- attr(fit, "epiexposure_family")

  if (is.null(spec) || !is.list(spec)) {
    stop("Missing 'epiexposure_spec' in fit. Ensure fit_epidlnm() stored it.")
  }

  if (is.null(nd0) || !is.data.frame(nd0) || nrow(nd0) < 1) {
    stop("Missing or invalid 'dat_template' in fit.")
  }

  if (is.null(vars_fit) || length(vars_fit) == 0) {
    stop("Missing 'epiexposure_vars' in fit.")
  }

  if (!any(grepl("^cb_", names(nd0)))) {
    stop("dat_template does not contain cb_* columns.")
  }

  # -----------------------
  # Helpers
  # -----------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  is_bayesian_model <- function(model) {
    inherits(model, "brmsfit") || inherits(model, "inla") || inherits(model, "bdlnm")
  }

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  get_cb_cols <- function(dat, v) {
    grep(paste0("^cb_", v, "_"), names(dat), value = TRUE)
  }

  build_cb_row <- function(profile_vec, spec_v) {

    lag_max <- as.integer(spec_v$lag_max)

    if (length(profile_vec) != lag_max + 1L) {
      stop("Profile length mismatch for variable (expected ", lag_max + 1L, ").")
    }

    cb <- dlnm::crossbasis(
      profile_vec,
      lag = lag_max,
      argvar = spec_v$argvar,
      arglag = spec_v$arglag
    )

    as.numeric(cb[lag_max + 1L, ])
  }

  # inverse link helper
  get_linkinv <- function(model, family_fit, default_identity = TRUE) {

    if (!is.null(model$family) && !is.null(model$family$linkinv)) {
      return(model$family$linkinv)
    }

    if (is.character(family_fit)) {
      fam <- family_fit
      if (fam %in% c("binomial", "beta")) {
        return(stats::binomial(link = "logit")$linkinv)
      }
      if (fam %in% c("poisson", "gamma")) {
        return(stats::poisson(link = "log")$linkinv)
      }
      if (fam %in% c("gaussian")) {
        return(stats::gaussian(link = "identity")$linkinv)
      }
    }

    if (inherits(family_fit, "family") && !is.null(family_fit$linkinv)) {
      return(family_fit$linkinv)
    }

    if (default_identity) {
      return(function(x) x)
    }

    NULL
  }

  # -----------------------
  # Normalize profiles
  # -----------------------
  if (is.numeric(profiles)) {
    if (length(vars_fit) != 1) {
      stop("Single numeric profile provided but model has multiple exposure variables.")
    }
    profiles <- setNames(list(profiles), vars_fit)

  } else if (!is.list(profiles) || is.null(names(profiles))) {
    stop("profiles must be a numeric vector or a named list.")
  }

  vars <- names(profiles)

  # -----------------------
  # Validate spec
  # -----------------------
  miss <- setdiff(vars, names(spec))
  if (length(miss) > 0) {
    stop("Missing spec for variable(s): ", paste(miss, collapse = ", "))
  }

  # -----------------------
  # Build newdata
  # -----------------------
  nd <- nd0

  if (!is.null(id)) {

    if (is.null(id_col)) {
      stop("Model has no random-effect column stored (id_col missing).")
    }

    if (!id_col %in% names(nd)) {
      stop("id column not found in dat_template: ", id_col)
    }

    nd <- nd[rep(1, length(id)), , drop = FALSE]
    nd[[id_col]] <- id
  }

  # -----------------------
  # Replace cb columns
  # -----------------------
  for (v in vars) {

    cols <- get_cb_cols(nd, v)

    if (length(cols) == 0) {
      stop("No cb_ columns found for variable: ", v)
    }

    cb_row <- build_cb_row(profiles[[v]], spec[[v]])

    if (length(cb_row) != length(cols)) {
      stop("Crossbasis dimension mismatch for variable: ", v)
    }

    if (nrow(nd) == 1) {
      nd[cols] <- as.list(cb_row)
    } else {
      tmp <- as.data.frame(matrix(
        rep(cb_row, each = nrow(nd)),
        nrow = nrow(nd),
        ncol = length(cb_row),
        byrow = FALSE
      ))
      names(tmp) <- cols
      nd[cols] <- tmp
    }
  }

  # -----------------------
  # Deterministic prediction helper
  # -----------------------
  point_predict <- function(model, newdata, re, allow_new_levels, type) {

    want_population <- identical(re, "population")

    if (inherits(model, "glmmTMB")) {

      re_form <- if (want_population) NA else NULL

      return(as.numeric(predict(
        model,
        newdata = newdata,
        type = type,
        re.form = re_form,
        allow.new.levels = allow_new_levels
      )))
    }

    if (inherits(model, "merMod")) {

      re_form <- if (want_population) NA else NULL

      return(as.numeric(predict(
        model,
        newdata = newdata,
        type = if (type == "link") "link" else "response",
        re.form = re_form,
        allow.new.levels = allow_new_levels
      )))
    }

    if (inherits(model, "brmsfit")) {

      re_formula <- if (want_population) NA else NULL

      pp <- brms::fitted(
        model,
        newdata = newdata,
        re_formula = re_formula,
        summary = TRUE
      )

      return(as.numeric(pp[, "Estimate"]))
    }

    if (inherits(model, "HLfit")) {

      re_form <- if (want_population) NA else NULL

      return(as.numeric(predict(
        model,
        newdata = newdata,
        type = if (type == "link") "link" else "response",
        re.form = re_form
      )))
    }

    if (inherits(model, "lme")) {

      level <- if (want_population) 0 else 1

      return(as.numeric(nlme::predict.lme(model, newdata = newdata, level = level)))
    }

    if (inherits(model, "gls")) {
      return(as.numeric(nlme::predict.gls(model, newdata = newdata)))
    }

    if (inherits(model, "gam")) {
      return(as.numeric(mgcv::predict.gam(
        model,
        newdata = newdata,
        type = if (type == "link") "link" else "response"
      )))
    }

    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      common <- intersect(names(beta), names(newdata))

      if (length(common) == 0) {
        stop("No matching covariates for INLA prediction.")
      }

      eta <- as.numeric(as.matrix(newdata[, common, drop = FALSE]) %*% beta[common])

      if (!want_population) {
        warning("INLA conditional prediction not implemented; using fixed-effects only.")
      }

      if (type == "link") return(eta)

      linkinv <- get_linkinv(model, family_fit)
      return(as.numeric(linkinv(eta)))
    }

    if (inherits(model, "bdlnm")) {

      if (is.null(model$coefficients.summary)) {
        stop("bdlnm object does not contain coefficient summaries.")
      }

      beta_mean <- model$coefficients.summary[, "mean"]
      common <- intersect(names(beta_mean), names(newdata))

      if (length(common) == 0) {
        common <- intersect(cb_cols_fit %||% character(0), names(newdata))
        common <- intersect(common, names(beta_mean))
      }

      if (length(common) == 0) {
        stop("Could not match bdlnm coefficient names to newdata columns.")
      }

      eta <- as.numeric(as.matrix(newdata[, common, drop = FALSE]) %*% beta_mean[common])

      if (type == "link") return(eta)

      linkinv <- get_linkinv(model, family_fit)
      return(as.numeric(linkinv(eta)))
    }

    # glm / lm / fallback
    as.numeric(predict(
      model,
      newdata = newdata,
      type = if (type == "link") "link" else "response"
    ))
  }

  # -----------------------
  # No uncertainty: preserve old behavior
  # -----------------------
  if (!uncertainty) {
    pred <- point_predict(fit, nd, re, allow_new_levels, type)

    out <- data.frame(prediction = as.numeric(pred))

    if (!is.null(id)) {
      out[[id_col]] <- id
      out <- out[, c(id_col, "prediction")]
    }

    return(out)
  }

  # -----------------------
  # Bayesian models
  # -----------------------
  if (inherits(fit, "brmsfit")) {

    want_population <- identical(re, "population")
    re_formula <- if (want_population) NA else NULL

    draws <- brms::fitted(
      fit,
      newdata = nd,
      re_formula = re_formula,
      summary = FALSE,
      ndraws = n_samples
    )

    draws <- as.matrix(draws)

    if (output == "samples") {
      if (ncol(draws) == 1) {
        out <- data.frame(sample = seq_len(nrow(draws)), prediction = as.numeric(draws[, 1]))
      } else {
        out_list <- vector("list", ncol(draws))
        for (j in seq_len(ncol(draws))) {
          tmp <- data.frame(sample = seq_len(nrow(draws)), prediction = as.numeric(draws[, j]))
          if (!is.null(id)) tmp[[id_col]] <- id[j]
          out_list[[j]] <- tmp
        }
        out <- do.call(rbind, out_list)
        if (!is.null(id)) out <- out[, c(id_col, "sample", "prediction")]
      }
      return(out)
    }

    mean_pred <- apply(draws, 2, mean)
    sd_pred   <- apply(draws, 2, stats::sd)
    q_pred    <- t(apply(draws, 2, safe_quantile))

    out <- data.frame(
      prediction = mean_pred,
      sd = sd_pred,
      lower = q_pred[, 1],
      upper = q_pred[, 2]
    )

    if (!is.null(id)) {
      out[[id_col]] <- id
      out <- out[, c(id_col, "prediction", "sd", "lower", "upper")]
    }

    return(out)
  }

  if (inherits(fit, "inla")) {

    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for INLA uncertainty quantification.")
    }

    if (!is.null(id)) {
      warning("INLA conditional uncertainty is not implemented; using fixed-effects only.")
    }

    posterior <- tryCatch(
      INLA::inla.posterior.sample(n = n_samples, result = fit),
      error = function(e) NULL
    )

    if (is.null(posterior)) {
      stop("INLA posterior samples could not be drawn. Ensure the model was fitted with control.compute = list(config = TRUE).")
    }

    beta_names <- rownames(fit$summary.fixed)
    if (is.null(beta_names)) stop("Could not determine fixed-effect names from INLA model.")

    common <- intersect(beta_names, names(nd))
    if (length(common) == 0) stop("No matching covariates for INLA uncertainty prediction.")

    beta_draws <- do.call(cbind, lapply(posterior, function(s) {
      latent <- s$latent
      names(latent) <- gsub(":1$", "", names(latent))
      latent[common]
    }))

    eta_draws <- t(as.matrix(nd[, common, drop = FALSE]) %*% beta_draws)

    if (type != "link") {
      linkinv <- get_linkinv(fit, family_fit)
      eta_draws <- apply(eta_draws, 2, linkinv)
    }

    if (output == "samples") {
      if (ncol(eta_draws) == 1) {
        out <- data.frame(sample = seq_len(nrow(eta_draws)), prediction = as.numeric(eta_draws[, 1]))
      } else {
        out_list <- vector("list", ncol(eta_draws))
        for (j in seq_len(ncol(eta_draws))) {
          tmp <- data.frame(sample = seq_len(nrow(eta_draws)), prediction = as.numeric(eta_draws[, j]))
          if (!is.null(id)) tmp[[id_col]] <- id[j]
          out_list[[j]] <- tmp
        }
        out <- do.call(rbind, out_list)
        if (!is.null(id)) out <- out[, c(id_col, "sample", "prediction")]
      }
      return(out)
    }

    mean_pred <- apply(eta_draws, 2, mean)
    sd_pred   <- apply(eta_draws, 2, stats::sd)
    q_pred    <- t(apply(eta_draws, 2, safe_quantile))

    out <- data.frame(
      prediction = mean_pred,
      sd = sd_pred,
      lower = q_pred[, 1],
      upper = q_pred[, 2]
    )

    if (!is.null(id)) {
      out[[id_col]] <- id
      out <- out[, c(id_col, "prediction", "sd", "lower", "upper")]
    }

    return(out)
  }

  if (inherits(fit, "bdlnm")) {

    if (is.null(fit$coefficients)) {
      stop("bdlnm object does not contain posterior coefficient samples.")
    }

    beta_draws <- fit$coefficients

    if (is.null(dim(beta_draws))) {
      beta_draws <- matrix(beta_draws, ncol = 1)
    }

    all_s <- seq_len(ncol(beta_draws))
    if (length(all_s) > n_samples) {
      set.seed(1)
      keep <- sort(sample(all_s, n_samples))
      beta_draws <- beta_draws[, keep, drop = FALSE]
    }

    common <- intersect(rownames(beta_draws), names(nd))

    if (length(common) == 0) {
      common <- intersect(cb_cols_fit %||% character(0), names(nd))
      common <- intersect(common, rownames(beta_draws))
    }

    if (length(common) == 0) {
      stop("Could not match bdlnm coefficient names to newdata columns.")
    }

    eta_draws <- t(as.matrix(nd[, common, drop = FALSE]) %*% beta_draws[common, , drop = FALSE])

    if (type != "link") {
      linkinv <- get_linkinv(fit, family_fit)
      eta_draws <- apply(eta_draws, 2, linkinv)
    }

    if (output == "samples") {
      if (ncol(eta_draws) == 1) {
        out <- data.frame(sample = seq_len(nrow(eta_draws)), prediction = as.numeric(eta_draws[, 1]))
      } else {
        out_list <- vector("list", ncol(eta_draws))
        for (j in seq_len(ncol(eta_draws))) {
          tmp <- data.frame(sample = seq_len(nrow(eta_draws)), prediction = as.numeric(eta_draws[, j]))
          if (!is.null(id)) tmp[[id_col]] <- id[j]
          out_list[[j]] <- tmp
        }
        out <- do.call(rbind, out_list)
        if (!is.null(id)) out <- out[, c(id_col, "sample", "prediction")]
      }
      return(out)
    }

    mean_pred <- apply(eta_draws, 2, mean)
    sd_pred   <- apply(eta_draws, 2, stats::sd)
    q_pred    <- t(apply(eta_draws, 2, safe_quantile))

    out <- data.frame(
      prediction = mean_pred,
      sd = sd_pred,
      lower = q_pred[, 1],
      upper = q_pred[, 2]
    )

    if (!is.null(id)) {
      out[[id_col]] <- id
      out <- out[, c(id_col, "prediction", "sd", "lower", "upper")]
    }

    return(out)
  }

  # -----------------------
  # Frequentist models: automatic normal approximation
  # -----------------------
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

    # glm / lm fallback
    return(list(
      beta = stats::coef(model),
      vcov = as.matrix(stats::vcov(model))
    ))
  }

  cv <- extract_coef_vcov(fit)
  beta_hat <- cv$beta
  V_hat <- cv$vcov

  common <- intersect(names(beta_hat), names(nd))
  if (length(common) == 0) {
    stop("Could not match model coefficients to newdata columns for uncertainty approximation.")
  }

  beta_hat <- beta_hat[common]
  V_hat <- V_hat[common, common, drop = FALSE]

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for frequentist uncertainty approximation.")
  }

  beta_draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat, Sigma = V_hat)

  if (is.null(dim(beta_draws))) {
    beta_draws <- matrix(beta_draws, nrow = n_samples, byrow = FALSE)
    colnames(beta_draws) <- names(beta_hat)
  }

  X <- as.matrix(nd[, common, drop = FALSE])

  eta_draws <- beta_draws %*% t(X)  # n_samples x n_obs

  if (!identical(re, "population") &&
      (inherits(fit, "glmmTMB") || inherits(fit, "merMod") || inherits(fit, "lme"))) {
    warning("For frequentist mixed models, uncertainty currently reflects fixed-effects approximation only.")
  }

  if (type != "link") {
    linkinv <- get_linkinv(fit, family_fit)
    eta_draws <- apply(eta_draws, 2, linkinv)
  }

  if (output == "samples") {
    if (ncol(eta_draws) == 1) {
      out <- data.frame(sample = seq_len(nrow(eta_draws)), prediction = as.numeric(eta_draws[, 1]))
    } else {
      out_list <- vector("list", ncol(eta_draws))
      for (j in seq_len(ncol(eta_draws))) {
        tmp <- data.frame(sample = seq_len(nrow(eta_draws)), prediction = as.numeric(eta_draws[, j]))
        if (!is.null(id)) tmp[[id_col]] <- id[j]
        out_list[[j]] <- tmp
      }
      out <- do.call(rbind, out_list)
      if (!is.null(id)) out <- out[, c(id_col, "sample", "prediction")]
    }
    return(out)
  }

  mean_pred <- apply(eta_draws, 2, mean)
  sd_pred   <- apply(eta_draws, 2, stats::sd)
  q_pred    <- t(apply(eta_draws, 2, safe_quantile))

  out <- data.frame(
    prediction = mean_pred,
    sd = sd_pred,
    lower = q_pred[, 1],
    upper = q_pred[, 2]
  )

  if (!is.null(id)) {
    out[[id_col]] <- id
    out <- out[, c(id_col, "prediction", "sd", "lower", "upper")]
  }

  return(out)
}
