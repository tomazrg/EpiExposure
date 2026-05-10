#' Reduce DLNM effects to one dimension
#'
#' @param fit Fitted model
#' @param cb_template Crossbasis template
#' @param type Reduction type ("overall","lag","var")
#' @param value Fixed lag or exposure value
#' @param family_choice Distribution
#'
#' @return data.frame
#' @export
reduce_effects <- function(fit,
                           cb_template,
                           type = c("overall", "lag", "var"),
                           value = NULL,
                           family_choice) {
  
  type <- match.arg(type)
  
  if (inherits(fit, "glmmTMB")) {
    beta <- fixef(fit)$cond
    V    <- vcov(fit)$cond
  } else {
    beta <- coef(fit)
    V    <- vcov(fit)
  }
  
  cr <- dlnm::crossreduce(
    cb_template,
    coef = beta,
    vcov = V,
    type = type,
    value = value
  )
  
  df <- data.frame(
    x = cr$predvar,
    eta = cr$fit,
    low = cr$low,
    high = cr$high
  )
  
  df$effect <- if (family_choice %in% c("beta","gamma","poisson")) {
    (exp(df$eta) - 1) * 100
  } else {
    df$eta
  }
  
  df
}