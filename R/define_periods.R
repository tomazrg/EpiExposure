#' Define epidemiological lag periods
#'
#' Creates lag periods to summarise cumulative DLNM effects.
#'
#' @param lag_max Maximum lag (integer)
#' @param cuts Optional numeric vector of cut points defining periods
#' @param prefix Character string used to label periods (default = "W")
#'
#' @return data.frame with columns:
#'   - period
#'   - lag_start
#'   - lag_end
#'
#' @export
define_periods <- function(lag_max, cuts = NULL, prefix = "W") {

  # -------------------------------
  # Input checks
  # -------------------------------
  if (!is.character(prefix) || length(prefix) != 1) {
    stop("`prefix` must be a single character string.")
  }

  if (!is.numeric(lag_max) || length(lag_max) != 1 || lag_max <= 0) {
    stop("`lag_max` must be a positive integer.")
  }

  lag_max <- as.integer(lag_max)

  if (!is.null(cuts)) {
    if (!is.numeric(cuts)) {
      stop("`cuts` must be numeric.")
    }
  }

  # -------------------------------
  # Helper
  # -------------------------------
  build_df <- function(starts, ends) {

    ids <- paste0(prefix, seq_along(starts))

    data.frame(
      period    = ids,
      lag_start = starts,
      lag_end   = ends,
      stringsAsFactors = FALSE
    )
  }

  # -------------------------------
  # Single period
  # -------------------------------
  if (is.null(cuts) || length(cuts) == 0) {
    return(build_df(0, lag_max))
  }

  # -------------------------------
  # Multiple periods
  # -------------------------------
  cuts <- sort(unique(as.integer(cuts)))
  cuts <- cuts[cuts > 0 & cuts < lag_max]

  if (length(cuts) == 0) {

    message(
      "No valid cut points found; using a single cumulative period (0 to lag_max)."
    )

    return(build_df(0, lag_max))
  }

  starts <- c(0, cuts + 1)
  ends   <- c(cuts, lag_max)

  build_df(starts, ends)
}
