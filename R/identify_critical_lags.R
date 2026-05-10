#' Identify critical lags based on daily effects#' Identifyr::group_by(lag) |>
dplyr::summarise(
  mean_effect = mean(effect, na.rm = TRUE),
  max_effect  = max(effect, na.rm = TRUE),
  .groups = "drop"
) |>
  dplyr::arrange(desc(abs(max_effect)))
}
#'
#' @param daily_df Output from summarise_effects(scale="daily")
#'
#' @return data.frame
#' @export
identify_critical_lags <- function(daily_df) {
  
  daily_df |>
    