#' Compare predicted outcomes between multiple exposure scenarios
#'
#' Compares predicted outcomes across two or more exposure scenarios,
#' using `predict_outcome()` as the computational backend.
#'
#' This function supports both deterministic comparisons and uncertainty
#' propagation. When `uncertainty = TRUE`, predictions are computed at the
#' sample level and comparisons are derived from these simulated values.
#'
#' If `output = "summary"`, the central estimate is computed as the median
#' of the sample-based distributions, and interval limits are derived from
#' empirical quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** although traditional terminology might suggest "mean",
#' all central estimates in summary outputs correspond to the *median*
#' when uncertainty is propagated, ensuring robustness under asymmetric
#' distributions.
#'
#' @param fit Fitted model (output of `fit_epidlnm()`).
#'
#' @param profiles Can be:
#'   - a named list of profiles (recommended), e.g.:
#'     `list(scenario1 = ..., scenario2 = ...)`, or
#'   - `NULL` when using `profiles1` and `profiles2` for backward compatibility.
#'
#' @param profiles1 (legacy) First profile (used only if `profiles` is NULL).
#'
#' @param profiles2 (legacy) Second profile (used only if `profiles` is NULL).
#'
#' @param re Character. Prediction level:
#'   - `"population"`: excludes random effects (default).
#'   - `"conditional"`: includes random effects where supported.
#'
#' @param id Optional character string indicating the column used as identifier.
#'
#' @param allow_new_levels Logical. Passed to `predict_outcome()`.
#'
#' @param type Character. Scale of prediction:
#'   `"response"` (default), `"link"`, or `"conditional"`.
#'
#' @param uncertainty Logical. If `TRUE`, uncertainty is propagated using
#'   sample-based predictions.
#'
#' @param output Character. Output type when `uncertainty = TRUE`:
#'   - `"summary"`: returns median-based summaries (default)
#'   - `"samples"`: returns all simulated samples
#'
#' @param n_samples Integer. Number of samples used for uncertainty propagation.
#'
#' @param eps Small positive constant used for numerical stability.
#'
#' @return A data.frame with pairwise comparisons including:
#'   - `scenario1`, `scenario2`
#'   - `pred1`, `pred2`
#'   - `diff`
#'   - `percent_change`
#'   - `ratio`
#'
#' When `uncertainty = TRUE`:
#' - `"samples"`: returns sample-level comparisons
#' - `"summary"`: returns median-based estimates, standard deviation,
#'   and empirical interval limits
#'
#' @details
#' When `uncertainty = TRUE`, this function always operates on simulated
#' prediction samples obtained from `predict_outcome(output = "samples")`.
#'
#' Summaries are then computed as:
#' - central estimate: median
#' - uncertainty intervals: empirical quantiles
#'
#' This approach ensures coherent uncertainty propagation for both
#' frequentist (simulation-based) and Bayesian (posterior-based) models,
#' avoiding incorrect analytic variance approximations.
#'
#' @export
compare_predictions <- function(
    fit,
    profiles = NULL,
    profiles1 = NULL,
    profiles2 = NULL,
    re = c("population", "conditional"),
    id = NULL,
    allow_new_levels = FALSE,
    type = c("response", "link", "conditional"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    eps = 1e-12
) {

  re <- match.arg(re)
  type <- match.arg(type)
  output <- match.arg(output)

  if (!is.numeric(eps) || eps <= 0) stop("eps must be positive.")

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # -----------------------
  # BACKWARD COMPAT
  # -----------------------
  if (is.null(profiles)) {

    if (is.null(profiles1) || is.null(profiles2)) {
      stop("Provide either `profiles` OR both `profiles1` and `profiles2`.")
    }

    profiles <- list(
      scenario1 = profiles1,
      scenario2 = profiles2
    )
  }

  if (!is.list(profiles) || is.null(names(profiles))) {
    stop("`profiles` must be a named list.")
  }

  if (length(profiles) < 2) {
    stop("At least two profiles are required.")
  }

  # ------------------------------------------------------------
  # ALWAYS USE SAMPLES IF UNCERTAINTY
  # ------------------------------------------------------------
  use_output <- if (uncertainty) "samples" else "summary"

  preds_list <- lapply(names(profiles), function(nm) {

    p <- predict_outcome(
      fit = fit,
      profiles = profiles[[nm]],
      re = re,
      id = id,
      allow_new_levels = allow_new_levels,
      type = type,
      uncertainty = uncertainty,
      output = use_output,
      n_samples = n_samples
    )

    p$scenario <- nm
    p
  })

  preds <- do.call(rbind, preds_list)

  # ============================================================
  # NO UNCERTAINTY
  # ============================================================
  if (!uncertainty) {

    wide <- reshape(
      preds,
      idvar = id %||% NULL,
      timevar = "scenario",
      direction = "wide"
    )

    combs <- combn(names(profiles), 2, simplify = FALSE)

    out_list <- lapply(combs, function(cb) {

      s1 <- cb[1]
      s2 <- cb[2]

      p1 <- wide[[paste0("prediction.", s1)]]
      p2 <- wide[[paste0("prediction.", s2)]]

      df <- data.frame(
        scenario1 = s1,
        scenario2 = s2,
        pred1 = p1,
        pred2 = p2
      )

      df$diff <- df$pred2 - df$pred1
      df$percent_change <- df$diff / (abs(df$pred1) + eps) * 100
      df$ratio <- df$pred2 / (df$pred1 + eps)

      if (!is.null(id)) {
        df[[id]] <- wide[[id]]
        df <- df[, c(id, "scenario1", "scenario2", "pred1", "pred2",
                     "diff", "percent_change", "ratio")]
      }

      df
    })

    return(do.call(rbind, out_list))
  }

  # ============================================================
  # UNCERTAINTY (SAMPLE-BASED)
  # ============================================================
  combs <- combn(names(profiles), 2, simplify = FALSE)

  out_samples <- lapply(combs, function(cb) {

    s1 <- cb[1]
    s2 <- cb[2]

    p1 <- preds[preds$scenario == s1, ]
    p2 <- preds[preds$scenario == s2, ]

    by_cols <- "sample"
    if (!is.null(id)) by_cols <- c(id, "sample")

    m <- merge(p1, p2, by = by_cols, suffixes = c("_1", "_2"))

    names(m)[names(m) == "prediction_1"] <- "pred1"
    names(m)[names(m) == "prediction_2"] <- "pred2"

    m$diff <- m$pred2 - m$pred1
    m$percent_change <- m$diff / (abs(m$pred1) + eps) * 100
    m$ratio <- m$pred2 / (m$pred1 + eps)

    m$scenario1 <- s1
    m$scenario2 <- s2

    m
  })

  samples_df <- do.call(rbind, out_samples)

  if (output == "samples") {
    return(samples_df)
  }

  # ============================================================
  # SUMMARY (MEDIAN-BASED)
  # ============================================================
  group_vars <- intersect(c(id, "scenario1", "scenario2"), names(samples_df))

  summary_df <- samples_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      pred1 = stats::median(.data$pred1, na.rm = TRUE),
      pred2 = stats::median(.data$pred2, na.rm = TRUE),

      diff = stats::median(.data$diff, na.rm = TRUE),
      diff_sd = stats::sd(.data$diff, na.rm = TRUE),
      diff_lower = stats::quantile(.data$diff, 0.025, na.rm = TRUE),
      diff_upper = stats::quantile(.data$diff, 0.975, na.rm = TRUE),

      percent_change = stats::median(.data$percent_change, na.rm = TRUE),
      ratio = stats::median(.data$ratio, na.rm = TRUE),

      .groups = "drop"
    )

  as.data.frame(summary_df)
}
