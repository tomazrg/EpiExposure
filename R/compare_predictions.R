#' Compare predicted outcomes between exposure scenarios
#'
#' Compares predicted outcomes across two or more exposure scenarios using
#' `predict_outcomes()` as the computational backend. Comparisons are performed
#' pairwise for all scenarios supplied in `profiles`.
#'
#' Exposure profiles must be supplied in chronological order, from the earliest
#' observation to the most recent observation. Validation of profile lengths,
#' reconstruction of the DLNM cross-basis, coefficient alignment, conversion
#' to the internal lag representation, model-engine handling, and uncertainty
#' propagation are delegated to `predict_outcomes()`.
#'
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param profiles Named list of exposure scenarios. Each top-level element
#'   represents one scenario and must contain the profile structure expected by
#'   `predict_outcomes()`. Profiles must be supplied in chronological order.
#' @param profiles1 Legacy profile for the first scenario. Used together with
#'   `profiles2` only when `profiles = NULL`.
#' @param profiles2 Legacy profile for the second scenario. Used together with
#'   `profiles1` only when `profiles = NULL`.
#' @param re Character. Random-effect prediction level: `"population"` or
#'   `"conditional"`.
#' @param id Optional character scalar identifying the grouping column returned
#'   by `predict_outcomes()`. When supplied, predictions and coefficient draws
#'   are matched by this column.
#' @param allow_new_levels Logical. Passed to `predict_outcomes()`.
#' @param type Character. Prediction scale passed to `predict_outcomes()`:
#'   `"response"` or `"link"`. Conditional predictions are selected with
#'   `re = "conditional"`, not through `type`.
#' @param uncertainty Logical. If `TRUE`, comparisons are calculated for every
#'   coefficient draw returned by `predict_outcomes()`.
#' @param output Character. When `uncertainty = TRUE`, `"summary"` returns
#'   median-based summaries and uncertainty intervals, while `"samples"`
#'   returns draw-level comparisons. When `uncertainty = FALSE`, a deterministic
#'   comparison is returned regardless of `output`.
#' @param n_samples Positive integer number of coefficient draws requested from
#'   `predict_outcomes()` when `uncertainty = TRUE`.
#' @param eps Positive finite numeric constant used to stabilize percentage and
#'   ratio calculations when the first prediction is close to zero.
#'
#' @return A data.frame.
#'
#' Without uncertainty, the output contains `scenario1`, `scenario2`, `pred1`,
#' `pred2`, `diff`, `percent_change`, and `ratio`, plus the column named by `id`
#' when applicable.
#'
#' With `uncertainty = TRUE` and `output = "samples"`, the output contains one
#' row per matched draw and scenario pair.
#'
#' With `uncertainty = TRUE` and `output = "summary"`, the output contains
#' medians, standard deviations, and empirical 95 percent intervals for `pred1`,
#' `pred2`, `diff`, `percent_change`, and `ratio`.
#'
#' @details
#' For a comparison between scenario 1 and scenario 2, the function calculates:
#'
#' \deqn{diff = pred_2 - pred_1}
#'
#' \deqn{percent\_change = 100(pred_2-pred_1)/(|pred_1|+\epsilon)}
#'
#' \deqn{ratio = pred_2/(pred_1+\epsilon)}
#'
#' Scenario order follows the order of names in `profiles`. Therefore, the sign
#' of `diff` and `percent_change`, and the direction of `ratio`, depend on that
#' ordering.
#'
#' Uncertainty summaries are calculated from draw-level comparisons. Temporary
#' names are used for medians inside `dplyr::summarise()` so that the original
#' draw columns are not overwritten before their standard deviations and
#' quantiles are calculated.
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
    type = c("response", "link"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    eps = 1e-12
) {
  re <- match.arg(re)
  type <- match.arg(type)
  output <- match.arg(output)

  if (is.null(fit)) stop("`fit` cannot be NULL.")

  if (!is.logical(allow_new_levels) ||
      length(allow_new_levels) != 1L ||
      is.na(allow_new_levels)) {
    stop("`allow_new_levels` must be TRUE or FALSE.")
  }

  if (!is.logical(uncertainty) ||
      length(uncertainty) != 1L ||
      is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }

  if (!is.numeric(n_samples) ||
      length(n_samples) != 1L ||
      !is.finite(n_samples) ||
      n_samples <= 0 ||
      n_samples != as.integer(n_samples)) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  if (uncertainty && n_samples < 2L) {
    stop("`n_samples` must be at least 2 when `uncertainty = TRUE`.")
  }

  if (!is.numeric(eps) ||
      length(eps) != 1L ||
      !is.finite(eps) ||
      eps <= 0) {
    stop("`eps` must be a single positive finite numeric value.")
  }

  if (!is.null(id) &&
      (!is.character(id) || length(id) != 1L || is.na(id) || id == "")) {
    stop("`id` must be NULL or a single non-empty column name.")
  }

  # Backward compatibility -------------------------------------------------
  if (is.null(profiles)) {
    if (is.null(profiles1) || is.null(profiles2)) {
      stop("Provide either `profiles` or both `profiles1` and `profiles2`.")
    }
    profiles <- list(
      scenario1 = profiles1,
      scenario2 = profiles2
    )
  }

  # Validate scenarios -----------------------------------------------------
  if (!is.list(profiles) || length(profiles) < 2L) {
    stop("`profiles` must be a named list containing at least two scenarios.")
  }

  if (is.null(names(profiles)) ||
      anyNA(names(profiles)) ||
      any(names(profiles) == "")) {
    stop("Every scenario in `profiles` must have a non-empty name.")
  }

  if (anyDuplicated(names(profiles))) {
    stop("Scenario names in `profiles` must be unique.")
  }

  scenario_names <- names(profiles)
  combinations <- utils::combn(scenario_names, 2, simplify = FALSE)

  # Prediction backend -----------------------------------------------------
  backend_output <- if (uncertainty) "samples" else "summary"

  predictions_list <- lapply(
    scenario_names,
    function(scenario_name) {
      prediction <- predict_outcomes(
        fit = fit,
        profiles = profiles[[scenario_name]],
        re = re,
        id = id,
        allow_new_levels = allow_new_levels,
        type = type,
        uncertainty = uncertainty,
        output = backend_output,
        n_samples = n_samples
      )

      if (!is.data.frame(prediction)) {
        prediction <- as.data.frame(prediction)
      }

      if (!"prediction" %in% names(prediction)) {
        stop(
          "`predict_outcomes()` did not return a `prediction` column for scenario '",
          scenario_name, "'."
        )
      }

      if (!is.numeric(prediction$prediction) ||
          any(!is.finite(prediction$prediction))) {
        stop(
          "Predictions for scenario '", scenario_name,
          "' must be finite numeric values."
        )
      }

      if (!is.null(id) && !id %in% names(prediction)) {
        stop(
          "`predict_outcomes()` did not return the requested id column '",
          id, "' for scenario '", scenario_name, "'."
        )
      }

      if (uncertainty && !"sample" %in% names(prediction)) {
        stop(
          "`predict_outcomes()` did not return a `sample` column for scenario '",
          scenario_name, "'."
        )
      }

      prediction$scenario <- scenario_name
      prediction
    }
  )
  names(predictions_list) <- scenario_names

  # Helpers ----------------------------------------------------------------
  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(0)
    stats::sd(x)
  }

  safe_quantile <- function(x, probability) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    stats::quantile(
      x,
      probs = probability,
      na.rm = TRUE,
      names = FALSE
    )
  }

  add_comparison_metrics <- function(data_frame) {
    data_frame$diff <- data_frame$pred2 - data_frame$pred1
    data_frame$percent_change <-
      data_frame$diff / (abs(data_frame$pred1) + eps) * 100
    data_frame$ratio <- data_frame$pred2 / (data_frame$pred1 + eps)
    data_frame
  }

  # Deterministic comparisons ----------------------------------------------
  if (!uncertainty) {
    output_list <- lapply(
      combinations,
      function(comparison) {
        scenario1 <- comparison[1]
        scenario2 <- comparison[2]
        prediction1 <- predictions_list[[scenario1]]
        prediction2 <- predictions_list[[scenario2]]

        if (is.null(id)) {
          if (nrow(prediction1) != nrow(prediction2)) {
            stop(
              "Scenarios '", scenario1, "' and '", scenario2,
              "' returned different numbers of prediction rows. Supply `id` ",
              "to match predictions explicitly."
            )
          }

          comparison_data <- data.frame(
            scenario1 = scenario1,
            scenario2 = scenario2,
            pred1 = prediction1$prediction,
            pred2 = prediction2$prediction,
            stringsAsFactors = FALSE
          )
        } else {
          if (anyDuplicated(prediction1[[id]]) ||
              anyDuplicated(prediction2[[id]])) {
            stop(
              "The id column must uniquely identify deterministic predictions ",
              "within each scenario."
            )
          }

          first <- prediction1[, c(id, "prediction"), drop = FALSE]
          second <- prediction2[, c(id, "prediction"), drop = FALSE]
          names(first)[names(first) == "prediction"] <- "pred1"
          names(second)[names(second) == "prediction"] <- "pred2"

          comparison_data <- merge(
            first,
            second,
            by = id,
            all = FALSE,
            sort = FALSE
          )

          if (nrow(comparison_data) != nrow(first) ||
              nrow(comparison_data) != nrow(second)) {
            stop(
              "Predictions could not be matched one-to-one by id for scenarios '",
              scenario1, "' and '", scenario2, "'."
            )
          }

          comparison_data$scenario1 <- scenario1
          comparison_data$scenario2 <- scenario2
          comparison_data <- comparison_data[
            ,
            c(id, "scenario1", "scenario2", "pred1", "pred2"),
            drop = FALSE
          ]
        }

        add_comparison_metrics(comparison_data)
      }
    )

    output_data <- do.call(rbind, output_list)
    rownames(output_data) <- NULL
    return(output_data)
  }

  # Draw-level uncertainty comparisons ------------------------------------
  sample_outputs <- lapply(
    combinations,
    function(comparison) {
      scenario1 <- comparison[1]
      scenario2 <- comparison[2]
      prediction1 <- predictions_list[[scenario1]]
      prediction2 <- predictions_list[[scenario2]]

      matching_columns <- "sample"
      if (!is.null(id)) matching_columns <- c(id, "sample")

      if (anyDuplicated(prediction1[matching_columns]) ||
          anyDuplicated(prediction2[matching_columns])) {
        stop(
          "The matching columns do not uniquely identify prediction draws ",
          "within each scenario."
        )
      }

      first <- prediction1[
        ,
        c(matching_columns, "prediction"),
        drop = FALSE
      ]
      second <- prediction2[
        ,
        c(matching_columns, "prediction"),
        drop = FALSE
      ]

      names(first)[names(first) == "prediction"] <- "pred1"
      names(second)[names(second) == "prediction"] <- "pred2"

      comparison_data <- merge(
        first,
        second,
        by = matching_columns,
        all = FALSE,
        sort = FALSE
      )

      if (nrow(comparison_data) != nrow(first) ||
          nrow(comparison_data) != nrow(second)) {
        stop(
          "Prediction draws could not be matched one-to-one for scenarios '",
          scenario1, "' and '", scenario2, "'."
        )
      }

      comparison_data$scenario1 <- scenario1
      comparison_data$scenario2 <- scenario2
      comparison_data <- add_comparison_metrics(comparison_data)

      output_columns <- c(
        if (!is.null(id)) id,
        "sample",
        "scenario1",
        "scenario2",
        "pred1",
        "pred2",
        "diff",
        "percent_change",
        "ratio"
      )

      comparison_data[, output_columns, drop = FALSE]
    }
  )

  samples_data <- do.call(rbind, sample_outputs)
  rownames(samples_data) <- NULL

  if (output == "samples") return(samples_data)

  # Median-based uncertainty summaries ------------------------------------
  grouping_variables <- c(
    if (!is.null(id)) id,
    "scenario1",
    "scenario2"
  )

  # IMPORTANT: central estimates use temporary names. Reusing `pred1`,
  # `pred2`, `diff`, `percent_change`, or `ratio` as an output name before
  # calculating SDs and quantiles would mask the original draw-level column
  # inside dplyr::summarise().
  summary_data <- samples_data |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(grouping_variables))
    ) |>
    dplyr::summarise(
      pred1_estimate = stats::median(.data$pred1, na.rm = TRUE),
      pred1_sd = safe_sd(.data$pred1),
      pred1_lower = safe_quantile(.data$pred1, 0.025),
      pred1_upper = safe_quantile(.data$pred1, 0.975),

      pred2_estimate = stats::median(.data$pred2, na.rm = TRUE),
      pred2_sd = safe_sd(.data$pred2),
      pred2_lower = safe_quantile(.data$pred2, 0.025),
      pred2_upper = safe_quantile(.data$pred2, 0.975),

      diff_estimate = stats::median(.data$diff, na.rm = TRUE),
      diff_sd = safe_sd(.data$diff),
      diff_lower = safe_quantile(.data$diff, 0.025),
      diff_upper = safe_quantile(.data$diff, 0.975),

      percent_change_estimate = stats::median(
        .data$percent_change,
        na.rm = TRUE
      ),
      percent_change_sd = safe_sd(.data$percent_change),
      percent_change_lower = safe_quantile(.data$percent_change, 0.025),
      percent_change_upper = safe_quantile(.data$percent_change, 0.975),

      ratio_estimate = stats::median(.data$ratio, na.rm = TRUE),
      ratio_sd = safe_sd(.data$ratio),
      ratio_lower = safe_quantile(.data$ratio, 0.025),
      ratio_upper = safe_quantile(.data$ratio, 0.975),

      .groups = "drop"
    )

  # Restore the established public output names after all summaries have been
  # calculated from the original draw-level columns.
  names(summary_data)[names(summary_data) == "pred1_estimate"] <- "pred1"
  names(summary_data)[names(summary_data) == "pred2_estimate"] <- "pred2"
  names(summary_data)[names(summary_data) == "diff_estimate"] <- "diff"
  names(summary_data)[names(summary_data) == "percent_change_estimate"] <-
    "percent_change"
  names(summary_data)[names(summary_data) == "ratio_estimate"] <- "ratio"

  summary_data <- summary_data[
    ,
    c(
      grouping_variables,
      "pred1", "pred1_sd", "pred1_lower", "pred1_upper",
      "pred2", "pred2_sd", "pred2_lower", "pred2_upper",
      "diff", "diff_sd", "diff_lower", "diff_upper",
      "percent_change", "percent_change_sd",
      "percent_change_lower", "percent_change_upper",
      "ratio", "ratio_sd", "ratio_lower", "ratio_upper"
    ),
    drop = FALSE
  ]

  as.data.frame(summary_data)
}
