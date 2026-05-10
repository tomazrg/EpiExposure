#' Recenter DLNM effects
#'
#' @param fit Fitted model
#' @param cb_template Crossbasis template
#' @param new_cen New centering value
#'
#' @return crosspred object
#' @export
recenter_effects <- function(fit, cb_template, new_cen) {
  
  if (inherits(fit, "glmmTMB")) {
    beta <- fixef(fit)$cond
    V    <- vcov(fit)$cond
  } else {
    beta <- coef(fit)
    V    <- vcov(fit)
  }
  
  dlnm::crosspred(
    cb_template,
    coef = beta,
    vcov = V,
    cen  = new_cen,
    bylag = 1
  )
}