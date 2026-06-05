#' Quantify lag contributions to the cumulative DLNM effect
#'
#' Supports:
#' 1) Deterministic input
#' 2) Summary input (effect with sd/CI)
#' 3) Samples input (with 'sample' column)
#'
#' @param daily_df Output from summarise_effects(scale = "daily")
#' @param lag_window Optional lag interval c(start, end)
#' @param absolute Logical (default TRUE)
#'
#' @return data.frame
#' @export
lag_contribution <- function(
    daily_df,
    lag_window = NULL,
    absolute = TRUE
) {

  if (!all(c("lag", "effect") %in% names(daily_df))) {
    stop("daily_df must contain columns 'lag' and 'effect'.")
  }

  has_samples <- "sample" %in% names(daily_df)
  has_var     <- "var" %in% names(daily_df)

  # ------------------------------------------------------------
  # Subset lag window
  # ------------------------------------------------------------
  if (!is.null(lag_window)) {
    if (length(lag_window) != 2) stop("lag_window must have length 2.")

    daily_df <- daily_df |>
      dplyr::filter(lag >= min(lag_window), lag <= max(lag_window))
  }

  # ------------------------------------------------------------
  # Helper: core computation (single dataset)
  # ------------------------------------------------------------
  compute_contribution <- function(df) {

    lag_effect <- df |>
      dplyr::group_by(lag) |>
      dplyr::summarise(
        eta_lag = sum(effect, na.rm = TRUE),
        .groups = "drop"
      )

    total_effect <- if (absolute) {
      sum(abs(lag_effect$eta_lag), na.rm = TRUE)
    } else {
      sum(lag_effect$eta_lag, na.rm = TRUE)
    }

    if (is.na(total_effect) || total_effect == 0) {
      return(NULL)
    }

    lag_effect |>
      dplyr::mutate(
        contribution = if (absolute) {
          abs(eta_lag) / total_effect
        } else {
          eta_lag / total_effect
        },
        contribution_percent = 100 * contribution
      )
  }

  # ============================================================
  # ✅ CASE 1 & 2: deterministic or summary
  # ============================================================
  if (!has_samples) {

    if (!has_var) {

      out <- compute_contribution(daily_df)
      if (is.null(out)) {
        warning("Total cumulative effect is zero.")
        return(NULL)
      }

      return(out)

    } else {

      out <- daily_df |>
        dplyr::group_by(var) |>
        dplyr::group_modify(~ compute_contribution(.x)) |>
        dplyr::ungroup()

      if (is.null(out) || nrow(out) == 0) {
        warning("Total cumulative effect is zero.")
        return(NULL)
      }

      return(out)
    }
  }

  # ============================================================
  # ✅ CASE 3: samples (uncertainty)
  # ============================================================

  if (!has_var) {

    # sample-level contributions
    lag_sample <- daily_df |>
      dplyr::group_by(sample) |>
      dplyr::group_modify(~ compute_contribution(.x)) |>
      dplyr::ungroup()

    lag_summary <- lag_sample |>
      dplyr::group_by(lag) |>
      dplyr::summarise(
        contribution_mean  = mean(contribution, na.rm = TRUE),
        contribution_sd    = stats::sd(contribution, na.rm = TRUE),
        contribution_lower = stats::quantile(contribution, 0.025, na.rm = TRUE),
        contribution_upper = stats::quantile(contribution, 0.975, na.rm = TRUE),

        contribution_percent_mean  = mean(contribution_percent, na.rm = TRUE),
        contribution_percent_sd    = stats::sd(contribution_percent, na.rm = TRUE),
        contribution_percent_lower = stats::quantile(contribution_percent, 0.025, na.rm = TRUE),
        contribution_percent_upper = stats::quantile(contribution_percent, 0.975, na.rm = TRUE),

        .groups = "drop"
      )

    return(lag_summary)

  } else {

    lag_sample <- daily_df |>
      dplyr::group_by(sample, var) |>
      dplyr::group_modify(~ compute_contribution(.x)) |>
      dplyr::ungroup()

    lag_summary <- lag_sample |>
      dplyr::group_by(var, lag) |>
      dplyr::summarise(
        contribution_mean  = mean(contribution, na.rm = TRUE),
        contribution_sd    = stats::sd(contribution, na.rm = TRUE),
        contribution_lower = stats::quantile(contribution, 0.025, na.rm = TRUE),
        contribution_upper = stats::quantile(contribution, 0.975, na.rm = TRUE),

        contribution_percent_mean  = mean(contribution_percent, na.rm = TRUE),
        contribution_percent_sd    = stats::sd(contribution_percent, na.rm = TRUE),
        contribution_percent_lower = stats::quantile(contribution_percent, 0.025, na.rm = TRUE),
        contribution_percent_upper = stats::quantile(contribution_percent, 0.975, na.rm = TRUE),

        .groups = "drop"
      )

    return(lag_summary)
  }
}
