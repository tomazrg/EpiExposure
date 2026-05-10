#' Compare accumulated effects between lag windows
#'
#' @param accumulated_df Output from summarise_effects(scale="accumulated")
#' @param window_ref Reference window id
#'
#' @return data.frame
#' @export
compare_lag_windows <- function(accumulated_df, window_ref = 1) {
  
  ref <- accumulated_df |>
    dplyr::filter(window == window_ref) |>
    dplyr::select(value, ref_effect = effect)
  
  accumulated_df |>
    dplyr::left_join(ref, by = "value") |>
    dplyr::mutate(diff = effect - ref_effect,
                  ratio = effect / ref_effect)
}