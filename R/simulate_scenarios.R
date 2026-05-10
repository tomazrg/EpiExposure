#' Simulate epidemiological scenarios
#'
#' @param fit Fitted model
#' @param cb_templates DLNM templates
#' @param wx_ref Weather reference data
#' @param var Variable to perturb
#' @param values Simulated values
#' @param lag_max Maximum lag
#' @param lag_band Active lags
#'
#' @return data.frame
#' @export
simulate_scenarios <- function(fit,
                               cb_templates,
                               wx_ref,
                               var,
                               values,
                               lag_max,
                               lag_band) {
  
  beta <- fixef(fit)$cond
  vars <- names(cb_templates)
  
  cb_names <- lapply(vars, function(v)
    grep(paste0("^cb_", v, "_"), names(beta), value = TRUE))
  names(cb_names) <- vars
  
  n_ref <- max(table(wx_ref$epi_id))
  
  centers <- sapply(vars, function(v)
    median(wx_ref[[v]], na.rm = TRUE))
  
  make_series <- function(val, cen) {
    x <- rep(cen, n_ref)
    idx <- n_ref - lag_band
    x[idx] <- val
    x
  }
  
  base_cb <- lapply(vars, function(v) {
    cb <- dlnm::crossbasis(
      rep(centers[v], n_ref),
      lag    = lag_max,
      argvar = attr(cb_templates[[v]], "argvar"),
      arglag = attr(cb_templates[[v]], "arglag")
    )
    as.numeric(cb[n_ref, ])
  })
  names(base_cb) <- vars
  
  purrr::map_dfr(values, function(val) {
    
    cb_now <- base_cb
    s <- make_series(val, centers[var])
    cb <- dlnm::crossbasis(
      s,
      lag    = lag_max,
      argvar = attr(cb_templates[[var]], "argvar"),
      arglag = attr(cb_templates[[var]], "arglag")
    )
    cb_now[[var]] <- as.numeric(cb[n_ref, ])
    
    eta <- beta["(Intercept)"]
    for (v in vars) {
      eta <- eta + sum(beta[cb_names[[v]]] * cb_now[[v]])
    }
    
    data.frame(
      value      = val,
      prediction = plogis(eta)
    )
  })
}