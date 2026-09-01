#' Reduce DLNM effects to one dimension
#'
#' Reduces a fitted DLNM exposure-lag-response association using
#' `dlnm::crossreduce()`. The cross-basis is reconstructed exclusively from the
#' exposure specification stored by `fit_epidlnm()`.
#'
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param data Long-format exposure data containing `group`, `time`, and `var`.
#' @param var Character scalar naming one fitted exposure variable.
#' @param group Character scalar naming the grouping column.
#' @param time Character scalar naming the chronological time column.
#' @param type Reduction type: `"overall"`, `"lag"`, or `"var"`.
#' @param value Finite numeric value required for `type = "lag"` or
#'   `type = "var"`.
#' @param scale Output scale: `"link"`, `"response"`, or `"percent"`.
#' @param uncertainty Logical. Propagate coefficient uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Positive integer number of draws. At least two are required
#'   when `uncertainty = TRUE`.
#' @param seed Optional finite integer used for reproducible sampling.
#'
#' @return A data.frame containing the reduction coordinate `x`, the reduced
#'   linear-predictor effect `eta`, transformed `effect`, and reduction metadata.
#'   Uncertainty summaries additionally contain `eta_sd`, `low`, `high`,
#'   `low_eff`, and `high_eff`. Sample output contains `sample`.
#'
#' @details
#' When `uncertainty = FALSE`, effects are calculated from the central
#' coefficient estimates supplied by the fitted engine. When
#' `uncertainty = TRUE`, effects are calculated for every coefficient draw and
#' summarized using the median and empirical 2.5% and 97.5% quantiles.
#'
#' For log-link models, `scale = "percent"` is a relative percentage change in
#' the expected outcome. For logit-link models, it is a relative exposure effect
#' on the linear-predictor scale and is not a direct percentage change in the
#' response. Percent scale is unavailable for other links.
#'
#' @export
reduce_effects <- function(
    fit,
    data,
    var,
    group = "epi_id",
    time = "time",
    type = c("overall", "lag", "var"),
    value = NULL,
    scale = c("percent", "response", "link"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    seed = NULL
) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  output <- match.arg(output)

  scalar_name <- function(x) is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(x)
  if (is.null(fit)) stop("`fit` cannot be NULL.")
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!scalar_name(var)) stop("`var` must be one non-empty variable name.")
  if (!scalar_name(group)) stop("`group` must be one non-empty column name.")
  if (!scalar_name(time)) stop("`time` must be one non-empty column name.")
  missing_cols <- setdiff(c(group, time, var), names(data))
  if (length(missing_cols)) stop("Columns missing from `data`: ", paste(missing_cols, collapse = ", "), ".")
  if (anyNA(data[[group]])) stop("The grouping column cannot contain missing values.")
  if (!is.numeric(data[[time]]) || any(!is.finite(data[[time]]))) stop("`time` must contain finite numeric values.")
  if (!is.numeric(data[[var]]) || any(!is.finite(data[[var]]))) stop("Exposure variable '", var, "' must contain finite numeric values.")
  if (!is.logical(uncertainty) || length(uncertainty) != 1L || is.na(uncertainty)) stop("`uncertainty` must be TRUE or FALSE.")
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) ||
      n_samples <= 0 || n_samples != as.integer(n_samples)) stop("`n_samples` must be a positive integer.")
  n_samples <- as.integer(n_samples)
  if (uncertainty && n_samples < 2L) stop("`n_samples` must be at least 2 when `uncertainty = TRUE`.")
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) stop("`seed` must be NULL or one finite integer.")
    set.seed(as.integer(seed))
  }
  if (type %in% c("lag", "var")) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) stop("`value` must be one finite numeric value when `type = 'lag'` or `type = 'var'`.")
  } else if (!is.null(value)) {
    warning("`value` is ignored when `type = 'overall'`.", call. = FALSE)
    value <- NULL
  }

  spec <- attr(fit, "epiexposure_spec")
  fitted_vars <- attr(fit, "epiexposure_vars")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")
  data_template <- attr(fit, "epiexposure_data_template")
  basis_objects <- attr(fit, "epiexposure_basis_objects")
  engine <- attr(fit, "epiexposure_engine")
  family_name <- attr(fit, "epiexposure_family_name")
  link_name <- attr(fit, "epiexposure_link")

  if (is.null(spec) || !is.list(spec) || is.null(spec[[var]])) stop("The fitted model does not contain a valid exposure specification for variable '", var, "'.")
  if (is.null(fitted_vars) || !length(fitted_vars)) {
    if (identical(engine, "bdlnm") && is.list(basis_objects) && length(basis_objects)) fitted_vars <- names(basis_objects)
    else stop("Missing `epiexposure_vars` metadata in `fit`.")
  }
  if (!var %in% fitted_vars) stop("Variable '", var, "' was not fitted by the model.")
  if (is.null(link_name) || !scalar_name(link_name)) stop("Missing or invalid `epiexposure_link` metadata in `fit`.")
  link_name <- tolower(link_name)
  if (identical(family_name, "ordinal")) stop("Ordinal reductions are not yet implemented in `reduce_effects()`.")
  if (scale == "percent" && !link_name %in% c("log", "logit")) stop("`scale = 'percent'` is supported only for log and logit links.")

  get_linkinv <- function(link) switch(
    tolower(link), identity = identity, log = exp, logit = stats::plogis,
    probit = stats::pnorm, cloglog = function(x) 1 - exp(-exp(x)),
    inverse = function(x) 1 / x,
    stop("Unsupported inverse-link function for link '", link, "'.")
  )
  linkinv <- get_linkinv(link_name)
  safe_q <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(c(NA_real_, NA_real_)); stats::quantile(x, c(.025, .975), names = FALSE) }
  safe_sd <- function(x) { x <- x[is.finite(x)]; if (length(x) <= 1L) return(0); stats::sd(x) }
  sort_cb <- function(x) { if (!length(x)) return(x); i <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x))); i[is.na(i)] <- seq_along(x)[is.na(i)]; x[order(i)] }
  match_brms <- function(ref, available) {
    prefixed <- paste0("b_", ref)
    if (all(prefixed %in% available)) return(prefixed)
    if (all(ref %in% available)) return(ref)
    stop("Could not match brms posterior draws to cross-basis coefficients.")
  }
  check_vcov <- function(V, nms) {
    if (is.null(rownames(V)) || is.null(colnames(V))) stop("The covariance matrix must contain row and column names.")
    miss <- setdiff(nms, intersect(rownames(V), colnames(V)))
    if (length(miss)) stop("Cross-basis coefficients absent from the covariance matrix: ", paste(miss, collapse = ", "), ".")
  }

  s <- spec[[var]]
  if (is.null(s$max_lag) || !is.numeric(s$max_lag) || any(!is.finite(s$max_lag))) stop("Invalid `max_lag` metadata for variable '", var, "'.")
  max_lag <- as.integer(max(s$max_lag))
  if (type == "lag" && (value < 0 || value > max_lag)) stop("For `type = 'lag'`, `value` must lie between 0 and max_lag.")
  if (is.null(s$argvar) || !is.list(s$argvar) || is.null(s$arglag) || !is.list(s$arglag)) stop("Invalid `argvar` or `arglag` metadata for variable '", var, "'.")

  ids <- unique(data[[group]])
  pooled <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    z <- data[data[[group]] == ids[[i]], , drop = FALSE]
    z <- z[order(z[[time]]), , drop = FALSE]
    if (anyDuplicated(z[[time]])) stop("Duplicated time values for group '", ids[[i]], "'.")
    if (nrow(z) < max_lag + 1L) stop("Group '", ids[[i]], "' has insufficient temporal coverage. Required observations: ", max_lag + 1L, ".")
    pooled[[i]] <- c(as.numeric(z[[var]]), rep(NA_real_, max_lag))
  }
  x_pool <- unlist(pooled, use.names = FALSE)
  cb <- dlnm::crossbasis(x_pool, lag = max_lag, argvar = s$argvar, arglag = s$arglag)
  p <- ncol(cb)

  extract_coef_vcov <- function(model) {
    if (inherits(model, "glmmTMB")) return(list(beta = glmmTMB::fixef(model)$cond, vcov = as.matrix(stats::vcov(model)$cond)))
    if (inherits(model, "merMod")) return(list(beta = lme4::fixef(model), vcov = as.matrix(stats::vcov(model))))
    if (is.list(model) && !is.null(model$gam) && inherits(model$gam, "gam")) return(list(beta = stats::coef(model$gam), vcov = as.matrix(stats::vcov(model$gam))))
    if (inherits(model, "lme")) return(list(beta = nlme::fixef(model), vcov = as.matrix(stats::vcov(model))))
    if (inherits(model, "gls") || inherits(model, "gam")) return(list(beta = stats::coef(model), vcov = as.matrix(stats::vcov(model))))
    if (inherits(model, "HLfit")) { V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL); if (is.null(V)) stop("Could not extract vcov from spaMM model."); return(list(beta = spaMM::fixef(model), vcov = V)) }
    if (inherits(model, "brmsfit")) { fe <- brms::fixef(model); b <- fe[, "Estimate"]; names(b) <- rownames(fe); return(list(beta = b, vcov = as.matrix(stats::vcov(model)))) }
    if (inherits(model, "inla")) { b <- model$summary.fixed$mean; if (is.null(names(b))) names(b) <- rownames(model$summary.fixed); V <- diag(model$summary.fixed$sd^2); dimnames(V) <- list(names(b), names(b)); return(list(beta = b, vcov = V)) }
    if (inherits(model, "bdlnm")) { if (is.null(model$coefficients.summary) || is.null(model$coefficients)) stop("The bdlnm model lacks coefficient summaries or draws."); b <- model$coefficients.summary[, "mean"]; if (is.null(names(b))) names(b) <- rownames(model$coefficients.summary); V <- stats::cov(t(model$coefficients)); return(list(beta = b, vcov = V)) }
    b <- stats::coef(model); V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL); if (!is.numeric(b) || is.null(V)) stop("Could not extract coefficients and covariance matrix."); list(beta = b, vcov = V)
  }
  cv <- extract_coef_vcov(fit)
  beta_full <- cv$beta
  V_full <- cv$vcov
  if (is.null(names(beta_full))) stop("The model coefficient vector has no names.")

  candidates <- character(0)
  if (!is.null(cb_cols_fit)) candidates <- sort_cb(grep(paste0("^cb_", var, "_"), cb_cols_fit, value = TRUE))
  if (!length(candidates) && is.data.frame(data_template)) candidates <- sort_cb(grep(paste0("^cb_", var, "_"), names(data_template), value = TRUE))
  if (!length(candidates)) candidates <- sort_cb(grep(paste0("^cb_", var, "_"), names(beta_full), value = TRUE))
  if (!length(candidates) && inherits(fit, "bdlnm")) stop("Could not map the stored bdlnm basis to named cross-basis coefficients for variable '", var, "'.")
  cb_names <- candidates[candidates %in% names(beta_full)]
  if (length(cb_names) != p) stop("Could not align fitted coefficients with the reconstructed cross-basis. Expected ", p, " coefficients but found ", length(cb_names), ".")
  check_vcov(V_full, cb_names)
  beta_sub <- beta_full[cb_names]
  V_sub <- V_full[cb_names, cb_names, drop = FALSE]

  extract_draws <- function(model) {
    if (inherits(model, "brmsfit")) {
      if (!requireNamespace("posterior", quietly = TRUE)) stop("Package 'posterior' is required for brms uncertainty.")
      d <- posterior::as_draws_matrix(model); dn <- match_brms(cb_names, colnames(d)); d <- as.matrix(d[, dn, drop = FALSE]); if (nrow(d) > n_samples) d <- d[sample(seq_len(nrow(d)), n_samples, replace = FALSE), , drop = FALSE]; colnames(d) <- cb_names; return(d)
    }
    if (inherits(model, "inla")) {
      if (!requireNamespace("INLA", quietly = TRUE)) stop("Package 'INLA' is required for INLA uncertainty.")
      post <- tryCatch(INLA::inla.posterior.sample(n = n_samples, result = model), error = function(e) NULL)
      if (is.null(post)) stop("INLA posterior samples could not be generated. Fit with `control.compute = list(config = TRUE)`. ")
      d <- do.call(rbind, lapply(post, function(one) { latent <- one$latent; names(latent) <- gsub(":1$", "", names(latent)); vals <- latent[cb_names]; if (anyNA(vals)) stop("Could not match INLA draws to cross-basis coefficients."); as.numeric(vals) })); colnames(d) <- cb_names; return(d)
    }
    if (inherits(model, "bdlnm")) {
      d <- model$coefficients; if (is.null(dim(d))) d <- matrix(d, ncol = 1L); if (is.null(rownames(d)) || !all(cb_names %in% rownames(d))) stop("Could not match bdlnm draws to cross-basis coefficients."); d <- d[cb_names, , drop = FALSE]; if (ncol(d) > n_samples) d <- d[, sample(seq_len(ncol(d)), n_samples, replace = FALSE), drop = FALSE]; out <- t(d); colnames(out) <- cb_names; return(out)
    }
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required for frequentist uncertainty.")
    d <- MASS::mvrnorm(n_samples, mu = beta_sub, Sigma = V_sub); if (is.null(dim(d))) d <- matrix(d, nrow = 1L); colnames(d) <- cb_names; d
  }

  cr_args <- function(coef, V = diag(0, length(coef))) {
    args <- list(basis = cb, coef = coef, vcov = V, type = type)
    if (type %in% c("lag", "var")) args$value <- value
    do.call(dlnm::crossreduce, args)
  }

  if (!uncertainty) {
    cr <- cr_args(beta_sub, V_sub)
    df <- data.frame(x = cr$predvar, eta = cr$fit, stringsAsFactors = FALSE)
  } else {
    draws <- extract_draws(fit)
    draw_sd <- apply(draws, 2, stats::sd)
    if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) warning("Near-zero coefficient-draw variability for variable '", var, "'. Intervals may collapse.", call. = FALSE)
    ref <- cr_args(beta_sub)
    x_ref <- ref$predvar
    eta_mat <- vapply(seq_len(nrow(draws)), function(i) cr_args(draws[i, ])$fit, numeric(length(x_ref)))
    if (output == "samples") {
      df <- data.frame(x = rep(x_ref, times = ncol(eta_mat)), eta = as.vector(eta_mat), sample = rep(seq_len(ncol(eta_mat)), each = length(x_ref)), stringsAsFactors = FALSE)
    } else {
      qs <- t(apply(eta_mat, 1, safe_q))
      df <- data.frame(x = x_ref, eta = apply(eta_mat, 1, stats::median, na.rm = TRUE), eta_sd = apply(eta_mat, 1, safe_sd), low = qs[, 1], high = qs[, 2], stringsAsFactors = FALSE)
    }
  }

  transform_effect <- function(x) {
    if (scale == "link") return(x)
    if (scale == "response") return(linkinv(x))
    (exp(x) - 1) * 100
  }
  df$effect <- transform_effect(df$eta)
  if ("low" %in% names(df)) {
    df$low_eff <- transform_effect(df$low)
    df$high_eff <- transform_effect(df$high)
  } else {
    df$low_eff <- NA_real_
    df$high_eff <- NA_real_
  }
  df$type <- type
  df$value <- value
  df$scale <- scale
  df$var <- var
  df
}
