#' Reduce DLNM effects to one dimension
#'
#' Reduces a DLNM exposure-lag surface into a one-dimensional
#' exposure-response or lag-response relationship.
#'
#' @param fit Fitted model object (glm, glmmTMB, gam, brms, INLA, spaMM, etc.)
#' @param cb_template Crossbasis template for one exposure
#' @param type Reduction type: "overall", "lag", or "var"
#' @param value Optional numeric value (required for "lag" or "var")
#' @param scale Output scale: "link", "response", "percent"
#'
#' @return data.frame with reduced DLNM effects
#'
#' @export
reduce_effects <- function(fit,
                           cb_template,
                           type = c("overall", "lag", "var"),
                           value = NULL,
                           scale = c("percent", "response", "link")) {

  type  <- match.arg(type)
  scale <- match.arg(scale)

  # ------------------------------------------------------------
  # ✅ Extract beta and vcov (ALL supported models)
  # ------------------------------------------------------------
  extract_beta_vcov <- function(fit) {

    # glmmTMB
    if (inherits(fit, "glmmTMB")) {
      return(list(
        beta = fixef(fit)$cond,
        V    = vcov(fit)$cond
      ))
    }

    # brms
    if (inherits(fit, "brmsfit")) {
      fe <- brms::fixef(fit)
      beta <- fe[, "Estimate"]
      V <- as.matrix(brms::vcov(fit))
      names(beta) <- rownames(fe)
      return(list(beta = beta, V = V))
    }

    # INLA
    if (inherits(fit, "inla")) {
      beta <- fit$summary.fixed$mean
      V <- diag(fit$summary.fixed$sd^2)
      rownames(V) <- names(beta)
      colnames(V) <- names(beta)
      return(list(beta = beta, V = V))
    }

    # spaMM
    if (inherits(fit, "HLfit")) {
      return(list(
        beta = spaMM::fixef(fit),
        V    = spaMM::vcov(fit)
      ))
    }

    # Default: glm, gam, gamm, gls
    return(list(
      beta = coef(fit),
      V    = vcov(fit)
    ))
  }

  bv   <- extract_beta_vcov(fit)
  beta <- bv$beta
  V    <- bv$V

  # ------------------------------------------------------------
  # ✅ KEEP ONLY DLNM TERMS (ESSENTIAL)
  # ------------------------------------------------------------
  idx <- grepl("^cb_", names(beta))

  if (!any(idx)) {
    stop("No DLNM (cb_) terms found in model coefficients.")
  }

  beta_sub <- beta[idx]
  V_sub    <- V[idx, idx, drop = FALSE]

  # ------------------------------------------------------------
  # ✅ DLNM reduction
  # ------------------------------------------------------------
  cr <- dlnm::crossreduce(
    cb_template,
    coef  = beta_sub,
    vcov  = V_sub,
    type  = type,
    value = value
  )

  # ------------------------------------------------------------
  # ✅ Build dataframe
  # ------------------------------------------------------------
  df <- data.frame(
    x   = cr$predvar,
    eta = cr$fit,
    low = cr$low,
    high = cr$high
  )

  # ------------------------------------------------------------
  # ✅ Scale transformation
  # ------------------------------------------------------------
  if (scale == "link") {

    df$effect <- df$eta
    df$low_eff <- df$low
    df$high_eff <- df$high

  } else if (scale == "response") {

    df$effect <- exp(df$eta)
    df$low_eff <- exp(df$low)
    df$high_eff <- exp(df$high)

  } else if (scale == "percent") {

    df$effect <- (exp(df$eta) - 1) * 100
    df$low_eff <- (exp(df$low) - 1) * 100
    df$high_eff <- (exp(df$high) - 1) * 100
  }

  # ------------------------------------------------------------
  # ✅ Add metadata
  # ------------------------------------------------------------
  df$type  <- type
  df$value <- value

  df
}
