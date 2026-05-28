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
#'
#' @return data.frame with predictions
#' @export
predict_outcome <- function(
    fit,
    profiles,
    re = c("population", "conditional"),
    id = NULL,
    allow_new_levels = FALSE,
    type = c("response", "link", "conditional")
) {

  re <- match.arg(re)
  type <- match.arg(type)

  # -----------------------
  # ✅ VALIDATIONS (robust but minimal)
  # -----------------------
  if (is.null(fit)) stop("`fit` cannot be NULL.")

  spec <- attr(fit, "epiexposure_spec")
  nd0  <- attr(fit, "epiexposure_dat_template")
  vars_fit <- attr(fit, "epiexposure_vars")
  id_col <- attr(fit, "epiexposure_id_col")

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
  # Helper
  # -----------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # -----------------------
  # Normalize profiles
  # -----------------------
  if (is.numeric(profiles)) {
    if (length(vars_fit) != 1) {
      stop("Single numeric profile provided but model has multiple exposure variables.")
    }
    profiles <- setNames(list(profiles), vars_fit)

  } else if (!is.list(profiles) || is.null(names(profiles))) {
    stop("profiles must be a numeric vector or named list.")
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
  # Functions
  # -----------------------
  get_cb_cols <- function(dat, v) {
    grep(paste0("^cb_", v, "_"), names(dat), value = TRUE)
  }

  build_cb_row <- function(profile_vec, spec_v) {

    lag_max <- as.integer(spec_v$lag_max)

    if (length(profile_vec) != lag_max + 1L) {
      stop("Profile length mismatch (expected ", lag_max + 1L, ").")
    }

    cb <- dlnm::crossbasis(
      profile_vec,
      lag = lag_max,
      argvar = spec_v$argvar,
      arglag = spec_v$arglag
    )

    as.numeric(cb[lag_max + 1L, ])
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

    nd[cols] <- as.list(cb_row)
  }

  # -----------------------
  # Prediction
  # -----------------------
  want_population <- identical(re, "population")
  pred <- NULL
  model <- fit

  if (inherits(model, "glmmTMB")) {

    re_form <- if (want_population) NA else NULL

    pred <- predict(
      model,
      newdata = nd,
      type = type,
      re.form = re_form,
      allow.new.levels = allow_new_levels
    )

  } else if (inherits(model, "merMod")) {

    re_form <- if (want_population) NA else NULL

    pred <- predict(
      model,
      newdata = nd,
      type = if (type == "link") "link" else "response",
      re.form = re_form,
      allow.new.levels = allow_new_levels
    )

  } else if (inherits(model, "brmsfit")) {

    re_formula <- if (want_population) NA else NULL

    pp <- brms::fitted(
      model,
      newdata = nd,
      re_formula = re_formula,
      summary = TRUE
    )

    pred <- pp[, "Estimate"]

  } else if (inherits(model, "HLfit")) {

    re_form <- if (want_population) NA else NULL

    pred <- predict(
      model,
      newdata = nd,
      type = if (type == "link") "link" else "response",
      re.form = re_form
    )

  } else if (inherits(model, "lme")) {

    level <- if (want_population) 0 else 1

    pred <- nlme::predict.lme(model, newdata = nd, level = level)

  } else if (inherits(model, "gls")) {

    pred <- nlme::predict.gls(model, newdata = nd)

  } else if (inherits(model, "gam")) {

    pred <- mgcv::predict.gam(
      model,
      newdata = nd,
      type = if (type == "link") "link" else "response"
    )

  } else if (inherits(model, "inla")) {

    beta <- model$summary.fixed$mean
    common <- intersect(names(beta), names(nd))

    if (length(common) == 0) {
      stop("No matching covariates for INLA prediction.")
    }

    eta <- as.numeric(as.matrix(nd[, common, drop = FALSE]) %*% beta[common])

    pred <- eta

    if (!want_population) {
      warning("INLA conditional prediction not implemented (fixed-effects only).")
    }

  } else {

    pred <- predict(
      model,
      newdata = nd,
      type = if (type == "link") "link" else "response"
    )
  }

  # -----------------------
  # Output
  # -----------------------
  out <- data.frame(prediction = as.numeric(pred))

  if (!is.null(id)) {
    out[[id_col]] <- id
    out <- out[, c(id_col, "prediction")]
  }

  return(out)
}
