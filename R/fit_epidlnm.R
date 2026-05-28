#' Fit DLNM inferential model
#'
#' @param dat Design matrix containing y_model and cb_* terms
#' @param model_engine Modeling engine
#'   ("glm","glmmTMB","gam","gamm","gls","spamm","brms","inla")
#' @param family Distribution (engine-specific; can be character or family object)
#' @param random_effect Random effect variable (optional)
#' @param epiexposure_spec Optional list describing how crossbasis was built.
#'   Recommended structure (by variable name):
#'   list(
#'     tmean = list(lag_max=85, argvar=list(...), arglag=list(...)),
#'     vpd   = list(lag_max=85, argvar=list(...), arglag=list(...))
#'   )
#'   If NULL, the model is fit normally but predict_outcome() will not be able
#'   to rebuild cb terms from a profile without this spec.
#' @param ... Additional arguments passed to the modeling engine
#'
#' @return Fitted model object (same class as before), with extra attributes
#' @export
fit_epidlnm <- function(dat,
                        model_engine,
                        family,
                        random_effect = NULL,
                        epiexposure_spec = NULL,
                        ...) {

  cb_cols <- grep("^cb_", names(dat), value = TRUE)
  rhs <- paste(cb_cols, collapse = " + ")

  # -------------------------------
  # Formula construction
  # -------------------------------
  fml_fixed <- paste0("y_model ~ 1 + ", rhs)

  fml <- if (!is.null(random_effect) &&
             model_engine %in% c("glmmTMB", "gamm", "spamm", "brms")) {
    as.formula(paste0(fml_fixed, " + (1|", random_effect, ")"))
  } else {
    as.formula(fml_fixed)
  }

  # -------------------------------
  # Helper: infer vars from cb_* names
  # cb_<var>_<k>
  # -------------------------------
  infer_vars <- function(cb_cols) {
    # extract between "cb_" and last "_<number>"
    # works for cb_tmean_1, cb_vpd_12, etc.
    v <- sub("^cb_", "", cb_cols)
    v <- sub("_[0-9]+$", "", v)
    unique(v)
  }
  vars_inferred <- infer_vars(cb_cols)

  # -------------------------------
  # Helper: attach metadata safely (attributes)
  # -------------------------------
  attach_epiexposure_meta <- function(model_obj) {
    # 1) template (one row is enough)
    dat_template <- dat[1, , drop = FALSE]

    attr(model_obj, "epiexposure_engine") <- model_engine
    attr(model_obj, "epiexposure_family") <- family
    attr(model_obj, "epiexposure_cb_cols") <- cb_cols
    attr(model_obj, "epiexposure_vars") <- vars_inferred
    attr(model_obj, "epiexposure_dat_template") <- dat_template

    # random effect grouping column name (if any)
    attr(model_obj, "epiexposure_id_col") <- if (!is.null(random_effect)) random_effect else NULL

    # DLNM spec for rebuilding cb terms from profiles (optional)
    attr(model_obj, "epiexposure_spec") <- epiexposure_spec

    model_obj
  }

  # -------------------------------
  # Frequentist models
  # -------------------------------

  if (model_engine == "glm") {

    fam <- switch(
      family,
      gaussian = gaussian(),
      poisson  = poisson(link = "log"),
      gamma    = Gamma(link = "log"),
      binomial = binomial(link = "logit")
    )

    mod <- glm(fml, data = dat, family = fam, ...)
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "glmmTMB") {

    fam <- switch(
      family,
      beta      = glmmTMB::beta_family(link = "logit"),
      gaussian  = gaussian(),
      poisson   = poisson(link = "log"),
      gamma     = Gamma(link = "log"),
      binomial  = binomial(link = "logit"),
      negbin    = glmmTMB::nbinom2()
    )

    mod <- glmmTMB::glmmTMB(fml, data = dat, family = fam, ...)
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "gam") {

    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for GAM.")
    }

    mod <- mgcv::gam(fml, data = dat, family = family, ...)
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "gamm") {

    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for GAMM.")
    }

    # mgcv::gamm returns a list with $gam and $lme etc.
    mod <- mgcv::gamm(fml, data = dat, family = family, ...)
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "gls") {

    mod <- nlme::gls(fml, data = dat, ...)
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "spamm") {

    if (!requireNamespace("spaMM", quietly = TRUE)) {
      stop("Package 'spaMM' is required for spatial mixed models.")
    }

    mod <- spaMM::fitme(fml, data = dat, family = family, ...)
    return(attach_epiexposure_meta(mod))
  }

  # -------------------------------
  # Bayesian models
  # -------------------------------

  if (model_engine == "brms") {

    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required for Bayesian models.")
    }

    mod <- brms::brm(
      formula = fml,
      data    = dat,
      family  = family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "inla") {

    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for INLA models.")
    }

    mod <- INLA::inla(
      formula = fml,
      data    = dat,
      family  = family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  stop("Unsupported model_engine.")
}
