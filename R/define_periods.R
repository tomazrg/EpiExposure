#' Define epidemiological lag periods
#'
#' Creates lag periods used to summarise cumulative DLNM effects.
#'
#' Lag periods are defined on the retrospective lag scale, where
#' lag 0 represents the most recent observation and `lag_max`
#' represents the oldest observation included in the exposure history.
#'
#' @param lag_max Maximum lag. Must be a positive integer.
#' @param cuts Optional numeric vector of cut points defining lag periods.
#' @param prefix Character string used to label periods.
#'   Default: `"W"`.
#'
#' @return A data.frame with columns:
#' \itemize{
#'   \item `period`
#'   \item `lag_start`
#'   \item `lag_end`
#' }
#'
#' @details
#' If `cuts = NULL`, a single cumulative lag period is returned:
#'
#' \preformatted{
#' lag_start = 0
#' lag_end   = lag_max
#' }
#'
#' If cut points are supplied, periods are created as contiguous,
#' non-overlapping lag intervals:
#'
#' \preformatted{
#' cuts = c(7, 14)
#'
#' W1 = 0 - 7
#' W2 = 8 - 14
#' W3 = 15 - lag_max
#' }
#'
#' Invalid cut points outside the interval
#' `[1, lag_max - 1]` are ignored.
#'
#' @export
define_periods <- function(
    lag_max,
    cuts = NULL,
    prefix = "W"
) {

  # =========================================================
  # Normalize lag dimension
  # =========================================================

  lag_max <- as.integer(
    max(lag_max)
  )

  # =========================================================
  # Validation
  # =========================================================

  if (
    is.na(lag_max)
  ) {
    stop(
      "`lag_max` must be a finite positive integer."
    )
  }

  if (
    !is.numeric(lag_max) ||
    length(lag_max) != 1L ||
    !is.finite(lag_max) ||
    lag_max <= 0
  ) {
    stop(
      "`lag_max` must be a positive integer."
    )
  }

  if (
    !is.character(prefix) ||
    length(prefix) != 1L ||
    is.na(prefix) ||
    prefix == ""
  ) {
    stop(
      "`prefix` must be a non-empty character string."
    )
  }

  if (!is.null(cuts)) {

    if (!is.numeric(cuts)) {
      stop(
        "`cuts` must be numeric."
      )
    }
  }

  # =========================================================
  # Helper
  # =========================================================

  build_df <- function(
    starts,
    ends
  ) {

    if (length(starts) != length(ends)) {

      stop(
        "Internal error: mismatched period boundaries."
      )
    }

    ids <- paste0(
      prefix,
      seq_along(starts)
    )

    data.frame(
      period = ids,
      lag_start = starts,
      lag_end = ends,
      stringsAsFactors = FALSE
    )
  }

  # =========================================================
  # Single cumulative period
  # =========================================================

  if (
    is.null(cuts) ||
    length(cuts) == 0L
  ) {

    return(
      build_df(
        starts = 0,
        ends = lag_max
      )
    )
  }

  # =========================================================
  # Multiple lag periods
  # =========================================================

  cuts <- sort(
    unique(
      as.integer(cuts)
    )
  )

  cuts <- cuts[
    cuts > 0 &
      cuts < lag_max
  ]

  if (length(cuts) == 0L) {

    message(
      paste0(
        "No valid cut points found; using a single ",
        "cumulative period (0 to lag_max)."
      )
    )

    return(
      build_df(
        starts = 0,
        ends = lag_max
      )
    )
  }

  starts <- c(
    0,
    cuts + 1L
  )

  ends <- c(
    cuts,
    lag_max
  )

  if (any(starts > ends)) {

    stop(
      "Invalid period definition (start > end)."
    )
  }

  build_df(
    starts = starts,
    ends = ends
  )
}
