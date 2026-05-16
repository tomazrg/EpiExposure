#' Predict full DLNM exposure-lag-response surface
#'
#' Generates the full DLNM surface (exposure × lag) based on
#' fitted model coefficients and cross-basis structure.
#'
#' @param fit Fitted model object from fit_epidlnm()
#' @param cb_template Crossbasis template for one exposure
#' @param wx_values Observed exposure values
#' @param ref Reference exposure definition:
#'   list(method = "median" | "percentile" | "fixed", value = NULL)
#' @param probs Quantiles defining exposure grid
#'
#' @return A crosspred object containing the exposure-lag-response surface
#'
#' @export
predict_surface <- function(fit,
                            cb_template,
                            wx_values,
                            ref = list(method = "median", value = NULL),
                            probs = seq(0.05, 0.95, by = 0.01)) {

  # ------------------------------------------------------------
  # ✅ Extract beta and vcov for ALL model types
  # ------------------------------------------------------------
  extract_beta_vcov <- function(fit) {

    if (inherits(fit, "glmmTMB")) {
      return(list(
        beta = fixef(fit)$cond,
        V    = vcov(fit)$cond
      ))
    }

    if (inherits(fit, "brmsfit")) {
      fe <- brms::fixef(fit)
      beta <- fe[, "Estimate"]
      V <- as.matrix(brms::vcov(fit))
      names(beta) <- rownames(fe)
      return(list(beta = beta, V = V))
    }

    if (inherits(fit, "inla")) {
      beta <- fit$summary.fixed$mean
      V <- diag(fit$summary.fixed$sd^2)
      rownames(V) <- names(beta)
      colnames(V) <- names(beta)
      return(list(beta = beta, V = V))
    }

    if (inherits(fit, "HLfit")) {
      return(list(
        beta = spaMM::fixef(fit),
        V    = spaMM::vcov(fit)
      ))
    }

    return(list(
      beta = coef(fit),
      V    = vcov(fit)
    ))
  }

  bv <- extract_beta_vcov(fit)
  beta <- bv$beta
  V    <- bv$V

  # ------------------------------------------------------------
  # ✅ Keep only DLNM terms
  # ------------------------------------------------------------
  idx <- grepl("^cb_", names(beta))
  beta_sub <- beta[idx]
  V_sub    <- V[idx, idx, drop = FALSE]

  # ------------------------------------------------------------
  # Exposure grid
  # ------------------------------------------------------------
  at_vals <- sort(unique(
    as.numeric(quantile(wx_values, probs = probs, na.rm = TRUE))
  ))

  # ------------------------------------------------------------
  # ✅ NEW: Reference (cen) choice
  # ------------------------------------------------------------
  cen <- switch(
    ref$method,
    median     = median(wx_values, na.rm = TRUE),
    percentile = quantile(wx_values, ref$value, na.rm = TRUE),
    fixed      = ref$value,
    stop("Invalid ref$method")
  )

  cen <- as.numeric(cen)

  # ------------------------------------------------------------
  # ✅ DLNM full surface
  # ------------------------------------------------------------
  cp <- dlnm::crosspred(
    cb_template,
    coef  = beta_sub,
    vcov  = V_sub,
    at    = at_vals,
    cen   = cen,
    bylag = 1
  )

  return(cp)
}
