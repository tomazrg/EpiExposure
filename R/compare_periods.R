#' Compare accumulated DLNM effects between periods
#'
#' Compares accumulated effects across epidemiological periods
#' using a reference period. Differences and ratios are computed
#' relative to the reference.
#'
#' @param accumulated_df Output from summarise_effects(scale = "accumulated").
#'   Must contain columns: value, period, effect.
#' @param period_ref Reference period (numeric index or character label, e.g. 1 or "W1")
#' @param eps Small constant to avoid division by zero in ratio calculation.
#'
#' @return A data.frame with comparisons relative to the reference period.
#'
#' @export
compare_periods <- function(
    accumulated_df,
    period_ref = 1,
    eps = 1e-10
) {

  # ------------------------------------------------------------
  # Basic input checks
  # ------------------------------------------------------------
  req_cols <- c("value", "period", "effect")

  if (!all(req_cols %in% names(accumulated_df))) {
    stop("accumulated_df must contain columns: value, period, effect.")
  }

  if (!is.numeric(accumulated_df$effect)) {
    stop("'effect' column must be numeric.")
  }

  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive numeric scalar.")
  }

  # ------------------------------------------------------------
  # ✅ Check duplicates (value, period)
  # ------------------------------------------------------------
  dup_check <- accumulated_df |>
    dplyr::count(value, period, name = "n") |>
    dplyr::filter(n > 1)

  if (nrow(dup_check) > 0) {
    stop("Found duplicated rows for (value, period). ",
         "Each combination must be unique.")
  }

  # ------------------------------------------------------------
  # ✅ Normalize period_ref if numeric
  # ------------------------------------------------------------
  if (is.numeric(period_ref)) {

    example_period <- accumulated_df$period[1]

    if (!grepl("[0-9]+$", example_period)) {
      stop("period labels must end in numeric index (e.g. W1, W2). ",
           "Provide period_ref as character instead.")
    }

    prefix <- sub("[0-9]+$", "", example_period)

    period_ref <- paste0(prefix, period_ref)
  }

  # ------------------------------------------------------------
  # Check reference period exists
  # ------------------------------------------------------------
  if (!period_ref %in% accumulated_df$period) {
    stop("period_ref not found in accumulated_df$period.")
  }

  # ------------------------------------------------------------
  # Extract reference values
  # ------------------------------------------------------------
  ref_df <- accumulated_df |>
    dplyr::filter(period == period_ref) |>
    dplyr::select(value, ref_effect = effect)

  # ------------------------------------------------------------
  # Join reference to all periods
  # ------------------------------------------------------------
  out <- accumulated_df |>
    dplyr::left_join(ref_df, by = "value")

  # ------------------------------------------------------------
  # Safety check after join
  # ------------------------------------------------------------
  if (any(!is.finite(out$ref_effect))) {
    stop("Join with reference failed for some values. ",
         "Ensure consistent `value` across periods.")
  }

  # ------------------------------------------------------------
  # Compute comparison metrics
  # ------------------------------------------------------------
  out <- out |>
    dplyr::mutate(
      reference = period_ref,
      diff = effect - ref_effect,
      ratio = effect / (ref_effect + eps),
      ratio_percent = (ratio - 1) * 100
    )

  # ------------------------------------------------------------
  # ✅ Order output
  # ------------------------------------------------------------
  out <- out |>
    dplyr::arrange(value, period)

  return(out)
}
