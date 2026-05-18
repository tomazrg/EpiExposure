#' Quantify lag contributions to the cumulative DLNM effect
#'
#' Decomposes the cumulative DLNM effect into relative contributions
#' of individual lags based on daily lag-specific effects.
#'
#' @param daily_df Output from summarise_effects(scale = "daily").
#'   Must contain columns: lag, effect.
#' @param lag_window Optional vector of length 2 specifying the lag
#'   interval c(lag_start, lag_end). If NULL, the full lag range is used.
#' @param absolute Logical. If TRUE (default), contributions are calculated
#'   using absolute effects to avoid sign cancellation.
#'
#' @return A data.frame with lag-specific contributions to the cumulative effect.
#'
#' @export
lag_contribution <- function(daily_df,
                             lag_window = NULL,
                             absolute = TRUE) {

  # ------------------------------------------------------------
  # Basic input checks
  # ------------------------------------------------------------
  stopifnot(all(c("lag", "effect") %in% names(daily_df)))

  # ------------------------------------------------------------
  # Select lag interval
  # ------------------------------------------------------------
  if (!is.null(lag_window)) {
    stopifnot(length(lag_window) == 2)

    daily_df <- daily_df |>
      dplyr::filter(
        lag >= min(lag_window),
        lag <= max(lag_window)
      )
  }

  # ------------------------------------------------------------
  # Collapse exposure dimension (sum over x)
  # ------------------------------------------------------------
  lag_effect <- daily_df |>
    dplyr::group_by(lag) |>
    dplyr::summarise(
      eta_lag = sum(effect, na.rm = TRUE),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # Compute total cumulative effect
  # ------------------------------------------------------------
  if (absolute) {
    total_effect <- sum(abs(lag_effect$eta_lag), na.rm = TRUE)
  } else {
    total_effect <- sum(lag_effect$eta_lag, na.rm = TRUE)
  }

  if (total_effect == 0) {
    warning("Total cumulative effect is zero; contributions cannot be computed.")
    return(NULL)
  }

  # ------------------------------------------------------------
  # Calculate lag contributions
  # ------------------------------------------------------------
  lag_effect |>
    dplyr::mutate(
      contribution = if (absolute) {
        abs(eta_lag) / total_effect
      } else {
        eta_lag / total_effect
      }
    ) |>
    dplyr::mutate(
      contribution_percent = 100 * contribution
    )
}
