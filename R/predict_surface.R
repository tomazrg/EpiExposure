#' Predict full DLNM lag-exposure surface
#'
#' @param fit Fitted model
#' @param cb_template Cross#' @param cen Centering value#' @param cb_template Crossbasis template
#' @param probs Quantiles
#'
#' @return list
#' @export
predict_surface <- function(fit,
                            cb_template,
                            wx_values,
                            cen = NULL,
                            probs = seq(0.05, 0.95, by = 0.01)) {
  
  if (inherits(fit, "glmmTMB")) {
    beta <- fixef(fit)$cond
    V    <- vcov(fit)$cond
  } else {
    beta <- coef(fit)
    V    <- vcov(fit)
  }
  
  at_vals <- quantile(wx_values, probs, na.rm = TRUE)
  if (is.null(cen)) cen <- median(wx_values, na.rm = TRUE)
  
  dlnm::crosspred(
    cb_template,
    coef = beta,
    vcov = V,
    at   = at_vals,
    cen  = cen,
    bylag = 1
  )
}
#' @param wx_values Exposure values
