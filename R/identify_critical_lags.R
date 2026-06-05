#' Identify critical lags based on daily DLNM effects
#'
#' Supports three input types:
#' 1) Deterministic (no uncertainty)
#' 2) Summary (with sd / CI)
#' 3) Samples (posterior / simulated draws)
#'
#' @param daily_df Output from summarise_effects(scale = "daily").
#'   Must contain at least: lag, effect.
#'   Optional columns: sample, var, effect_sd, effect_lower, effect_upper
#' @param metric Summary metric used to quantify lag importance.
#'   One of "max", "mean", or "absmean".
#'
#' @return A data.frame with lag importance ranking
#' @export
identify_critical_lags <- function(
    daily_df,
    metric = c("max", "mean", "absmean")
) {

  metric <- match.arg(metric)

  # ------------------------------------------------------------
  # Basic input check
  # ------------------------------------------------------------
  if (!all(c("lag", "effect") %in% names(daily_df))) {
    stop("daily_df must contain columns 'lag' and 'effect'.")
  }

  has_samples <- "sample" %in% names(daily_df)
  has_var     <- "var" %in% names(daily_df)

  # ------------------------------------------------------------
  # Helper: compute lag summary
  # ------------------------------------------------------------
  summarise_lag <- function(df) {

    out <- df |>
      dplyr::group_by(lag) |>
      dplyr::summarise(
        mean_effect    = mean(effect, na.rm = TRUE),
        max_effect     = max(effect, na.rm = TRUE),
        absmean_effect = mean(abs(effect), na.rm = TRUE),
        .groups = "drop"
      )

    out$score <- dplyr::case_when(
      metric == "max"     ~ abs(out$max_effect),
      metric == "mean"    ~ abs(out$mean_effect),
      metric == "absmean" ~ out$absmean_effect
    )

    out
  }

  # ============================================================
  # ✅ CASE 1 & 2: deterministic OR summary
  # ============================================================
  if (!has_samples) {

    if (!has_var) {

      # single variable
      lag_summary <- summarise_lag(daily_df)

      return(
        lag_summary |>
          dplyr::arrange(dplyr::desc(score))
      )

    } else {

      # multi-variable
      lag_summary <- daily_df |>
        dplyr::group_by(var, lag) |>
        dplyr::group_modify(~ summarise_lag(.x)) |>
        dplyr::ungroup()

      return(
        lag_summary |>
          dplyr::arrange(var, dplyr::desc(score))
      )
    }
  }

  # ============================================================
  # ✅ CASE 3: samples (uncertainty)
  # ============================================================

  if (!has_var) {

    # -----------------------
    # sample-level summary
    # -----------------------
    lag_sample <- daily_df |>
      dplyr::group_by(sample, lag) |>
      dplyr::group_modify(~ summarise_lag(.x)) |>
      dplyr::ungroup()

    # -----------------------
    # aggregate across samples
    # -----------------------
    lag_summary <- lag_sample |>
      dplyr::group_by(lag) |>
      dplyr::summarise(
        score_mean  = mean(score, na.rm = TRUE),
        score_sd    = stats::sd(score, na.rm = TRUE),
        score_lower = stats::quantile(score, 0.025, na.rm = TRUE),
        score_upper = stats::quantile(score, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    return(
      lag_summary |>
        dplyr::arrange(dplyr::desc(score_mean))
    )

  } else {

    # -----------------------
    # multivariable samples
    # -----------------------
    lag_sample <- daily_df |>
      dplyr::group_by(sample, var, lag) |>
      dplyr::group_modify(~ summarise_lag(.x)) |>
      dplyr::ungroup()

    lag_summary <- lag_sample |>
      dplyr::group_by(var, lag) |>
      dplyr::summarise(
        score_mean  = mean(score, na.rm = TRUE),
        score_sd    = stats::sd(score, na.rm = TRUE),
        score_lower = stats::quantile(score, 0.025, na.rm = TRUE),
        score_upper = stats::quantile(score, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    return(
      lag_summary |>
        dplyr::arrange(var, dplyr::desc(score_mean))
    )
  }
}
