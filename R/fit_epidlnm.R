#' Fit DLNM inferential model
#'
#' @param dat Design matrix containing y_model and cb_* terms
#' @param model_engine Modeling engine
#'   ("glm","glmmTMB","gam","gamm","gls","spamm","brms","inla")
#' @param family Distribution
#' @param random_effect Random effect variable (optional)
#' @param ... Additional arguments passed to the modeling engine
#'
#' @return Fitted model object
#' @export
fit_epidlnm <- function(dat,
                        model_engine,
                        family,
                        random_effect = NULL,
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

    return(glm(fml, data = dat, family = fam, ...))
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

    return(glmmTMB::glmmTMB(fml, data = dat, family = fam, ...))
  }

  if (model_engine == "gam") {

    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for GAM.")
    }

    return(mgcv::gam(fml, data = dat, family = family, ...))
  }

  if (model_engine == "gamm") {

    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for GAMM.")
    }

    return(mgcv::gamm(fml, data = dat, family = family, ...))
  }

  if (model_engine == "gls") {

    return(nlme::gls(fml, data = dat, ...))
  }

  if (model_engine == "spamm") {

    if (!requireNamespace("spaMM", quietly = TRUE)) {
      stop("Package 'spaMM' is required for spatial mixed models.")
    }

    return(spaMM::fitme(fml, data = dat, family = family, ...))
  }

  # -------------------------------
  # Bayesian models
  # -------------------------------

  if (model_engine == "brms") {

    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required for Bayesian models.")
    }

    return(brms::brm(
      formula = fml,
      data    = dat,
      family  = family,
      ...
    ))
  }

  if (model_engine == "inla") {

    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for INLA models.")
    }

    return(INLA::inla(
      formula = fml,
      data    = dat,
      family  = family,
      ...
    ))
  }

  stop("Unsupported model_engine.")
}
