#' Fit DLNM inferential model
#'
#' @param dat Design matrix
#' @param y_var Response variable
#' @param model_engine "glmmTMB","glm","gls"
#' @param family Distribution
#' @param random_effect Random effect variable
#'
#' @return fitted model
#' @export
fit_epidlnm <- function(dat,
                        y_var,
                        model_engine,
                        family,
                        random_effect = NULL) {
  
  cb_cols <- grep("^cb_", names(dat), value = TRUE)
  rhs <- paste(cb_cols, collapse = " + ")
  
  fml <- if (!is.null(random_effect)) {
    as.formula(paste0("y_model ~ 1 + ", rhs, " + (1|", random_effect, ")"))
  } else {
    as.formula(paste0("y_model ~ 1 + ", rhs))
  }
  
  if (model_engine == "glmmTMB") {
    , family = fam)    fam <- switch(
  } else if (model_engine == "glm") {
    fam <- switch(
      family,
      gamma    = Gamma(link = "log"),
      poisson  = poisson(link = "log"),
      gaussian = gaussian()
    )
    glm(fml, data = dat, family = fam)
  } else {
    nlme::gls(fml, data = dat)
  }
}
family,
beta     = glmmTMB::beta_family(link = "logit"),
gamma    = Gamma(link = "log"),
poisson  = poisson(link = "log"),
gaussian = gaussian()
    )
    