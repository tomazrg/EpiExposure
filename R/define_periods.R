#' Define epidemiological lag periods
#'
#' Creates lag periods used to summarise cumulative DLNM effects.
#'
#' Lag periods are defined on the retrospective lag scale, where
#' lag 0 represents the most recent observation and `max_lag`
#' represents the oldest observation included in the exposure history.
#'
#' @param max_lag Maximum lag. Must be a positive integer.
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
#' lag_end   = max_lag
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
#' W3 = 15 - max_lag
#' }
#'
#' Invalid cut points outside the interval
#' `[1, max_lag - 1]` are ignored.
#'
#' @export
define_periods <- function(
    max_lag,
    cuts = NULL,
    prefix = "W"
) {

  # =========================================================
  # Normalize lag dimension
  # =========================================================

  max_lag <- as.integer(
    max(max_lag)
  )

  # =========================================================
  # Validation
  # =========================================================

  if (
    is.na(max_lag)
  ) {
    stop(
      "`max_lag` must be a finite positive integer."
    )
  }

  if (
    !is.numeric(max_lag) ||
    length(max_lag) != 1L ||
    !is.finite(max_lag) ||
    max_lag <= 0
  ) {
    stop(
      "`max_lag` must be a positive integer."
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
        ends = max_lag
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
      cuts < max_lag
  ]

  if (length(cuts) == 0L) {

    message(
      paste0(
        "No valid cut points found; using a single ",
        "cumulative period (0 to max_lag)."
      )
    )

    return(
      build_df(
        starts = 0,
        ends = max_lag
      )
    )
  }

  starts <- c(
    0,
    cuts + 1L
  )

  ends <- c(
    cuts,
    max_lag
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
