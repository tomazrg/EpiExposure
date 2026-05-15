#' Identify critical lags based on daily DLNM effects
#'
#' Identifies lags with the strongest disease association by summarising
#' daily DLNM effects across exposure levels.
#'
#' @param daily_df Output from summarise_effects(scale = "daily").
#'   Must contain columns: lag, effect.
#' @param metric Summary metric used to quantify lag importance.
#'   One of "max", "mean", or "absmean".
#'
#' @return A data.frame with one row per lag and summary statistics
#'   used to identify critical lags.
#'
#' @export
identify_critical_lags <- function(daily_df,
                                   metric = c("max", "mean", "absmean")) {

  metric <- match.arg(metric)

  # Basic input check
  stopifnot(all(c("lag", "effect") %in% names(daily_df)))

  # ------------------------------------------------------------
  # Summarise daily effects by lag
  # ------------------------------------------------------------
  lag_summary <- daily_df |>
    dplyr::group_by(lag) |>
    dplyr::summarise(
      mean_effect = mean(effect, na.rm = TRUE),
      max_effect  = max(effect, na.rm = TRUE),
      absmean_effect = mean(abs(effect), na.rm = TRUE),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # Select metric used to identify critical lags
  # ------------------------------------------------------------
  lag_summary <- lag_summary |>
    dplyr::mutate(
      score = dplyr::case_when(
        metric == "max"     ~ abs(max_effect),
        metric == "mean"    ~ abs(mean_effect),
        metric == "absmean" ~ absmean_effect
      )
    )

  # ------------------------------------------------------------
  # Rank lags by importance
  # ------------------------------------------------------------
  lag_summary |>
    dplyr::arrange(dplyr::desc(score))
}
