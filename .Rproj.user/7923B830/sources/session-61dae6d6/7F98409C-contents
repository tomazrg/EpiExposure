#' Compare accumulated DLNM effects between lag windows
#'
#' Compares accumulated effects across epidemiological lag windows
#' using a reference window. Differences and ratios are computed
#' relative to the reference.
#'
#' @param accumulated_df Output from summarise_effects(scale = "accumulated").
#'   Must contain columns: value, window, effect.
#' @param window_ref Reference window (numeric index or character label, e.g. 1 or "W1")
#' @param eps Small constant to avoid division by zero in ratio calculation.
#'
#' @return A data.frame with pairwise comparisons relative to the reference window.
#'
#' @export
compare_lag_windows <- function(accumulated_df,
                                window_ref = 1,
                                eps = 1e-10) {

  # ------------------------------------------------------------
  # Basic input checks
  # ------------------------------------------------------------
  stopifnot(all(c("value", "window", "effect") %in% names(accumulated_df)))

  # ------------------------------------------------------------
  # ✅ NEW: Handle numeric → character conversion
  # ------------------------------------------------------------
  if (is.numeric(window_ref)) {

    # extract prefix automatically from data (ex: "W")
    prefix <- gsub("[0-9]+$", "", accumulated_df$window[1])

    window_ref <- paste0(prefix, window_ref)
  }

  # ------------------------------------------------------------
  # Check if reference exists
  # ------------------------------------------------------------
  if (!window_ref %in% accumulated_df$window) {
    stop("window_ref not found in accumulated_df$window")
  }

  # ------------------------------------------------------------
  # Extract reference window
  # ------------------------------------------------------------
  ref_df <- accumulated_df |>
    dplyr::filter(window == window_ref) |>
    dplyr::select(value, ref_effect = effect)

  # ------------------------------------------------------------
  # Join with all windows
  # ------------------------------------------------------------
  out <- accumulated_df |>
    dplyr::left_join(ref_df, by = "value")

  # ------------------------------------------------------------
  # Compute comparisons
  # ------------------------------------------------------------
  out |>
    dplyr::mutate(
      diff = effect - ref_effect,
      ratio = effect / (ref_effect + eps),
      ratio_percent = (ratio - 1) * 100
    )
}
