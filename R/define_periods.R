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

  # =========================================================
  # ✅ AJUSTE 1 — DIMENSÃO TEMPORAL (CRÍTICO)
  # =========================================================
  # 🔵 garante consistência com todo o pacote (evita c(0, 85))
  lag_max <- as.integer(max(lag_max))

  # =========================================================
  # ✅ CHECKS
  # =========================================================
  if (!is.character(prefix) || length(prefix) != 1) {
    stop("`prefix` must be a single character string.")
  }

  # 🔵 reforço: inteiro positivo
  if (!is.numeric(lag_max) || length(lag_max) != 1 || !is.finite(lag_max) || lag_max <= 0) {
    stop("`lag_max` must be a positive integer.")
  }

  if (!is.null(cuts)) {
    if (!is.numeric(cuts)) {
      stop("`cuts` must be numeric.")
    }
  }

  # =========================================================
  # ✅ HELPER
  # =========================================================
  build_df <- function(starts, ends) {

    # 🔵 checagem de consistência interna
    if (length(starts) != length(ends)) {
      stop("Internal error: mismatched period boundaries.")
    }

    ids <- paste0(prefix, seq_along(starts))

    data.frame(
      period    = ids,
      lag_start = starts,
      lag_end   = ends,
      stringsAsFactors = FALSE
    )
  }

  # =========================================================
  # ✅ SINGLE PERIOD
  # =========================================================
  if (is.null(cuts) || length(cuts) == 0) {
    return(build_df(0, lag_max))
  }

  # =========================================================
  # ✅ MULTIPLE PERIODS
  # =========================================================
  cuts <- sort(unique(as.integer(cuts)))

  # 🔵 garantir cortes válidos dentro do intervalo
  cuts <- cuts[cuts > 0 & cuts < lag_max]

  if (length(cuts) == 0) {

    message(
      "No valid cut points found; using a single cumulative period (0 to lag_max)."
    )

    return(build_df(0, lag_max))
  }

  # 🔵 definição consistente (inclusive)
  starts <- c(0, cuts + 1)
  ends   <- c(cuts, lag_max)

  # 🔵 checagem final de consistência
  if (any(starts > ends)) {
    stop("Invalid period definition (start > end).")
  }

  build_df(starts, ends)
}
