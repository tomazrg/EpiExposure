#' Fit DLNM inferential model
#'
#' @param dat Design matrix containing y_model and cb_* terms
#' @param model_engine Modeling engine
#'   ("glm","glmmTMB","gam","gamm","gls","spamm","brms","inla","bdlnm")
#' @param family Distribution (engine-specific; can be character or family object)
#' @param random_effect Random effect variable (optional)
#' @param epiexposure_spec Optional list describing how crossbasis was built.
#'   Recommended structure (by variable name):
#'   list(
#'     tmean = list(lag_max=85, argvar=list(...), arglag=list(...)),
#'     vpd   = list(lag_max=85, argvar=list(...), arglag=list(...))
#'   )
#'   If NULL, the model is fit normally but downstream prediction from profiles
#'   will not be available.
#' @param basis_objects Optional named list of original crossbasis/onebasis objects.
#'   Required when model_engine = "bdlnm", because bdlnm::bdlnm() needs the basis
#'   objects explicitly in the formula environment.
#' @param ... Additional arguments passed to the modeling engine
#'
#' @return Fitted model object (same class as before), with extra attributes
#' @export
fit_epidlnm <- function(dat,
                        model_engine,
                        family,
                        random_effect = NULL,
                        epiexposure_spec = NULL,
                        basis_objects = NULL,
                        ...) {

  # -------------------------------
  # Basic validation
  # -------------------------------
  if (!is.data.frame(dat)) {
    stop("`dat` must be a data.frame.")
  }

  if (!"y_model" %in% names(dat)) {
    stop("`dat` must contain a column named 'y_model'.")
  }

  model_engine <- match.arg(
    model_engine,
    choices = c("glm", "glmmTMB", "gam", "gamm", "gls", "spamm", "brms", "inla", "bdlnm")
  )

  cb_cols <- grep("^cb_", names(dat), value = TRUE)

  if (length(cb_cols) == 0 && model_engine != "bdlnm") {
    stop("No cb_* columns found in `dat`.")
  }

  # -------------------------------
  # Helper: infer vars from cb_* names
  # cb_<var>_<k>
  # -------------------------------
  infer_vars <- function(cb_cols) {
    v <- sub("^cb_", "", cb_cols)
    v <- sub("_[0-9]+$", "", v)
    unique(v)
  }

  vars_inferred <- if (length(cb_cols) > 0) infer_vars(cb_cols) else character(0)

  # -------------------------------
  # Helper: attach metadata safely (attributes)
  # -------------------------------
  attach_epiexposure_meta <- function(model_obj) {

    dat_template <- dat[1, , drop = FALSE]

    attr(model_obj, "epiexposure_engine") <- model_engine
    attr(model_obj, "epiexposure_family") <- family
    attr(model_obj, "epiexposure_cb_cols") <- cb_cols
    attr(model_obj, "epiexposure_vars") <- vars_inferred
    attr(model_obj, "epiexposure_dat_template") <- dat_template
    attr(model_obj, "epiexposure_id_col") <- if (!is.null(random_effect)) random_effect else NULL
    attr(model_obj, "epiexposure_spec") <- epiexposure_spec
    attr(model_obj, "epiexposure_basis_objects") <- basis_objects

    model_obj
  }

  # -------------------------------
  # Formula construction for matrix-based engines
  # -------------------------------
  rhs <- if (length(cb_cols) > 0) paste(cb_cols, collapse = " + ") else ""

  fml_fixed <- if (nzchar(rhs)) {
    paste0("y_model ~ 1 + ", rhs)
  } else {
    "y_model ~ 1"
  }

  # Engines with lme4-style random effects
  fml <- if (!is.null(random_effect) &&
             model_engine %in% c("glmmTMB", "gamm", "spamm", "brms")) {
    as.formula(paste0(fml_fixed, " + (1|", random_effect, ")"))
  } else {
    as.formula(fml_fixed)
  }

  # -------------------------------
  # Frequentist models
  # -------------------------------

  if (model_engine == "glm") {

    fam <- if (inherits(family, "family")) {
      family
    } else {
      switch(
        family,
        gaussian = gaussian(),
        poisson  = poisson(link = "log"),
        gamma    = Gamma(link = "log"),
        binomial = binomial(link = "logit"),
        stop("Unsupported family for glm.")
      )
    }

    mod <- glm(fml, data = dat, family = fam, ...)
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "glmmTMB") {

    fam <- if (inherits(family, "family")) {
      family
    } else {
      switch(
        family,
        beta      = glmmTMB::beta_family(link = "logit"),
        gaussian  = gaussian(),
        poisson   = poisson(link = "log"),
        gamma     = Gamma(link = "log"),
        binomial  = binomial(link = "logit"),
        negbin    = glmmTMB::nbinom2(),
        stop("Unsupported family for glmmTMB.")
      )
    }

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
  # Bayesian models - brms
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

  # -------------------------------
  # Bayesian models - INLA
  # random effect as f(random_effect, model = "iid")
  # -------------------------------
  if (model_engine == "inla") {

    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for INLA models.")
    }

    fml_inla <- if (!is.null(random_effect)) {
      as.formula(paste0(fml_fixed, " + f(", random_effect, ", model = 'iid')"))
    } else {
      as.formula(fml_fixed)
    }

    mod <- INLA::inla(
      formula = fml_inla,
      data    = dat,
      family  = family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  # -------------------------------
  # Bayesian DLNM backend - bdlnm
  # IMPORTANT:
  # bdlnm needs the original basis objects in the formula, not only cb_* columns.
  # -------------------------------
  if (model_engine == "bdlnm") {

    if (!requireNamespace("bdlnm", quietly = TRUE)) {
      stop("Package 'bdlnm' is required for model_engine = 'bdlnm'.")
    }

    if (is.null(basis_objects) || !is.list(basis_objects) || is.null(names(basis_objects))) {
      stop(
        "For model_engine = 'bdlnm', you must provide `basis_objects` as a named list ",
        "of original crossbasis/onebasis objects (e.g., output of define_exposure())."
      )
    }

    # validate basis names
    basis_names <- names(basis_objects)
    if (length(basis_names) == 0) {
      stop("`basis_objects` is empty.")
    }

    # build formula using basis object names, not cb_* columns
    rhs_bdlnm <- paste(basis_names, collapse = " + ")

    fml_bdlnm <- if (!is.null(random_effect)) {
      # INLA-style random effect, since bdlnm uses INLA internally
      as.formula(paste0("y_model ~ 1 + ", rhs_bdlnm, " + f(", random_effect, ", model = 'iid')"))
    } else {
      as.formula(paste0("y_model ~ 1 + ", rhs_bdlnm))
    }

    # create environment containing basis objects
    eval_env <- new.env(parent = parent.frame())
    for (nm in basis_names) {
      assign(nm, basis_objects[[nm]], envir = eval_env)
    }
    environment(fml_bdlnm) <- eval_env

    mod <- bdlnm::bdlnm(
      formula = fml_bdlnm,
      data    = dat,
      family  = family,
      ...
    )

    return(attach_epiexposure_meta(mod))
  }

  stop("Unsupported model_engine.")
}
