#' Compute lag contribution to cumulative effect
#'
#' @param daily_df Output from summarise_effects(scale="daily")
#'
#' @return data.frame
#' @export
lag_contribution <- function(daily_df) {
  
  total <- sum(daily_df$effect, na.rm = TRUE)
  
  daily_df |>
    dplyr::group_by(lag) |>
    dplyr::summarise(
      contribution = sum(effect, na.rm = TRUE) / total,
      .groups = "drop"
    )
}