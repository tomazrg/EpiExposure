#' Summarise DLNM effects
#'
#' Summarises estimated DLNM effects as daily (lag-specific) or accumulated
#' effects across epidemiological lag windows. The function is model-agnostic
#' and can be applied to both frequentist and Bayesian regression engines.
#'
#' @param fit A fitted model object returned by \code{fit_epidlnm()}.
#' @param cb_template Cross-basis template for a single exposure variable.
#' @param wx_values Observed exposure values used to define the prediction grid.
#' @param family_choice Distribution of the response variable. One of
#'   \code{"beta"}, \code{"gamma"}, \code{"poisson"}, \code{"binomial"},
#'   \code{"negbin"}, or \code{"gaussian"}.
#' @param scale Effect scale to summarise. Either \code{"daily"} or
#'   \code{"accumulated"}.
#' @param lag_windows Epidemiological lag windows defined by
#'   \code{define_lag_windows()}. Required if \code{scale = "accumulated"}.
#' @param probs Quantiles used to define the exposure prediction grid.
#' @param ref Reference exposure definition. A list with elements
#'   \code{method} (\code{"median"}, \code{"percentile"}, or \code{"fixed"})
#'   and \code{value}.
#' @param effect_measure Scale used to report effects. One of
#'   \code{"percent"}, \code{"ratio"}, or \code{"linear"}.
#'
#' @return A data.frame containing summarised DLNM effects.
#'
#' @export
summarise_effects <- function(fit,
                              cb_template,
                              wx_values,
                              family_choice,
                              scale = c("daily", "accumulated"),
                              lag_windows = NULL,
                              probs = seq(0.05, 0.95, by = 0.01),
                              ref = list(method = "median", value = NULL),
                              effect_measure = c("percent", "ratio", "linear")) {

  scale <- match.arg(scale)
  effect_measure <- match.arg(effect_measure)

  stopifnot(
    family_choice %in%
      c("beta", "gamma", "poisson", "binomial", "negbin", "gaussian")
  )

  if (scale == "accumulated") {
    if (is.null(lag_windows)) {
      stop("lag_windows must be provided when scale = 'accumulated'")
    }
    stopifnot(all(c("lag_start", "lag_end", "window_id") %in% names(lag_windows)))
  }

  # ------------------------------------------------------------
  # 1) Extract fixed-effect coefficients and VCOV (ALL engines)
  # ------------------------------------------------------------
  extract_beta_vcov <- function(fit) {

    # glmmTMB
    if (inherits(fit, "glmmTMB")) {
      return(list(
        beta = fixef(fit)$cond,
        V    = vcov(fit)$cond
      ))
    }

    # brms (Bayesian)
    if (inherits(fit, "brmsfit")) {
      fe <- brms::fixef(fit)
      beta <- fe[, "Estimate"]
      V <- as.matrix(brms::vcov(fit))
      names(beta) <- rownames(fe)
      return(list(beta = beta, V = V))
    }

    # INLA (Bayesian)
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

    # glm, gam, gamm, gls
    return(list(
      beta = coef(fit),
      V    = vcov(fit)
    ))
  }

  bv <- extract_beta_vcov(fit)
  beta <- bv$beta
  V    <- bv$V

  idx <- grepl("^cb_", names(beta))
  beta_sub <- beta[idx]
  V_sub    <- V[idx, idx, drop = FALSE]

  # ------------------------------------------------------------
  # 2) Exposure prediction grid
  # ------------------------------------------------------------
  at_vals <- sort(unique(
    as.numeric(quantile(wx_values, probs = probs, na.rm = TRUE))
  ))

  # ------------------------------------------------------------
  # 3) Reference (cen)
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
  # 4) DLNM prediction
  # ------------------------------------------------------------
  cp <- dlnm::crosspred(
    cb_template,
    coef  = beta_sub,
    vcov  = V_sub,
    at    = at_vals,
    cen   = cen,
    bylag = 1
  )

  mf <- cp$matfit
  lag_idx <- suppressWarnings(as.integer(gsub("lag", "", colnames(mf))))
  if (anyNA(lag_idx)) lag_idx <- 0:(ncol(mf) - 1)

  # ------------------------------------------------------------
  # Effect transformation
  # ------------------------------------------------------------
  transform_effect <- function(eta) {
    if (effect_measure == "linear") return(eta)
    if (effect_measure == "ratio")  return(exp(eta))
    return((exp(eta) - 1) * 100)
  }

  # ------------------------------------------------------------
  # 5) DAILY EFFECTS
  # ------------------------------------------------------------
  if (scale == "daily") {

    df <- as.data.frame(mf)
    colnames(df) <- paste0("lag_", lag_idx)

    return(
      df |>
        dplyr::mutate(value = at_vals) |>
        tidyr::pivot_longer(
          cols = starts_with("lag_"),
          names_to = "lag",
          values_to = "eta"
        ) |>
        dplyr::mutate(
          lag    = as.integer(gsub("lag_", "", lag)),
          effect = transform_effect(eta),
          scale  = "daily"
        )
    )
  }

  # ------------------------------------------------------------
  # 6) ACCUMULATED EFFECTS
  # ------------------------------------------------------------
  if (scale == "accumulated") {

    return(
      purrr::map_dfr(seq_len(nrow(lag_windows)), function(i) {

        sel <- which(
          lag_idx >= lag_windows$lag_start[i] &
            lag_idx <= lag_windows$lag_end[i]
        )

        if (length(sel) == 0) {
          warning(
            "Lag window ", lag_windows$window_id[i],
            " contains no lags; skipping."
          )
          return(NULL)
        }

        eta_sum <- rowSums(mf[, sel, drop = FALSE])

        data.frame(
          value  = at_vals,
          window = lag_windows$window_id[i],
          effect = transform_effect(eta_sum),
          scale  = "accumulated"
        )
      })
    )
  }
}
