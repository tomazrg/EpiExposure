#' Summarise DLNM effects (daily or accumulated)#' model
#' @param cb_template Crossbasis template for one variable
#' @param wx_values Observed exposure values
#' @param family_choice Distribution ("beta","gamma","poisson","gaussian")
#' @param scale Effect scale: "daily" or "accumulated"
#' @param lag_windows Required if scale = "accumulated"
#' @param probs Quantiles used to define exposure grid
#' @param ref Reference definition:
#'   list(method = "median" | "percentile" | "fixed", value = NULL)
#'
#' @return data.frame
#'
#' @export
summarise_effects <- function(fit,
                              cb_template,
                              wx_values,
                              family_choice,
                              scale = c("daily", "accumulated"),
                              lag_windows = NULL,
                              probs = seq(0.05, 0.95, by = 0.01),
                              ref = list(method = "median", value = NULL)) {

  scale <- match.arg(scale)

  # ------------------------------------------------------------
  # 1) Extract coefficients and vcov
  # ------------------------------------------------------------
  if (inherits(fit, "glmmTMB")) {
    beta <- fixef(fit)$cond
    V    <- vcov(fit)$cond
  } else {
    beta <- coef(fit)
    V    <- vcov(fit)
  }

  idx <- grepl("^cb_", names(beta))
  beta_sub <- beta[idx]
  V_sub    <- V[idx, idx, drop = FALSE]

  # ------------------------------------------------------------
  # 2) Exposure grid
  # ------------------------------------------------------------
  at_vals <- sort(unique(
    as.numeric(quantile(wx_values, probs = probs, na.rm = TRUE))
  ))

  # ------------------------------------------------------------
  # 3) Reference (cen) selection
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
  # 4) DLNM prediction (lag-specific surface)
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
  # 5) DAILY EFFECT
  # ------------------------------------------------------------
  if (scale == "daily") {

    df <- as.data.frame(mf)
    colnames(df) <- paste0("lag_", lag_idx)

    return(
      df |>
        dplyr::mutate(value = at_vals) |>
        tidyr::pivot_longer(
          cols      = starts_with("lag_"),
          names_to  = "lag",
          values_to = "eta"
        ) |>
        dplyr::mutate(
          lag = as.integer(gsub("lag_", "", lag)),
          effect = dplyr::case_when(
            family_choice %in% c("beta", "gamma", "poisson") ~
              (exp(eta) - 1) * 100,
            TRUE ~ eta
          )
        )
    )
  }

  # ------------------------------------------------------------
  # 6) ACCUMULATED EFFECT
  # ------------------------------------------------------------
  if (scale == "accumulated") {

    if (is.null(lag_windows)) {
      stop("lag_windows must be provided when scale = 'accumulated'")
    }

    return(
      purrr::map_dfr(seq_len(nrow(lag_windows)), function(i) {

        sel <- which(
          lag_idx >= lag_windows$lag_start[i] &
            lag_idx <= lag_windows$lag_end[i]
        )

        eta_sum <- rowSums(mf[, sel, drop = FALSE])

        data.frame(
          value  = at_vals,
          window = lag_windows$window_id[i],
          effect = dplyr::case_when(
            family_choice %in% c("beta", "gamma", "poisson") ~
              (exp(eta_sum) - 1) * 100,
            TRUE ~ eta_sum
          )
        )
      })
    )
  }
}

#'
