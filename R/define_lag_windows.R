#' Define epidemiological lag windows
#'
#' Creates lag windows to summarise cumulative DLNM effects.
#' If no cut points are provided, a single window from
#' lag 0 to lag_max is returned.
#'
#' @param lag_max Maximum lag (integer)
#' @param cuts Optional numeric vector of cut points defining windows
#'
#' @return data.frame with columns:
#'   - window_id
#'   - lag_start
#'   - lag_end
#'
#' @export
define_lag_windows <- function(lag_max, cuts = NULL) {

  # -------------------------------
  # Basic input checks (optional)
  # -------------------------------
  stopifnot(is.numeric(lag_max), lag_max > 0)

  if (!is.null(cuts)) {
    stopifnot(is.numeric(cuts))
  }

  # -------------------------------
  # Single cumulative window
  # -------------------------------
  if (is.null(cuts) || length(cuts) == 0) {

    return(
      data.frame(
        window_id = 1,
        lag_start = 0,
        lag_end   = lag_max
      )
    )
  }

  # -------------------------------
  # Multiple windows
  # -------------------------------
  cuts <- sort(unique(as.integer(cuts)))
  cuts <- cuts[cuts > 0 & cuts < lag_max]

  # If cuts collapse after filtering
  if (length(cuts) == 0) {

    message(
      "No valid cut points found; using a single cumulative lag window (0–lag_max)."
    )

    return(
      data.frame(
        window_id = 1,
        lag_start = 0,
        lag_end   = lag_max
      )
    )
  }

  starts <- c(0, cuts + 1)
  ends   <- c(cuts, lag_max)

  data.frame(
    window_id = seq_along(starts),
    lag_start = starts,
    lag_end   = ends
  )
}
