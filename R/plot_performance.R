#' Plot out-of-fold model performance
#'
#' Creates performance plots for models evaluated by `find_bestfit()`.
#' For non-binary outcomes, the default is an observed-versus-predicted
#' plot. For binomial outcomes, the default is a receiver operating
#' characteristic (ROC) curve.
#'
#' The function retrieves out-of-fold predictions from the
#' `"predictions_by_model"` attribute, or from the `"predictions"`
#' attribute when the model-specific list is unavailable. Model-level
#' validation metrics are retrieved from columns of the main
#' `find_bestfit()` result.
#'
#' @param object A non-empty data frame returned by `find_bestfit()`. It must
#'   contain a unique `model_id` column, the requested metric columns, and
#'   retained out-of-fold predictions.
#' @param model_id Optional vector of model identifiers to plot. If `NULL`,
#'   every retained model with available predictions is plotted.
#' @param metrics Unique character vector naming validation metrics stored as
#'   numeric columns of `object`. Metrics are displayed in the supplied order.
#'   Common non-binary metrics include `"CCC"`, `"Cb"`, `"rho"`, `"RMSE"`,
#'   and `"MAE"`. Binomial metrics include `"ROC_AUC"`, `"Brier"`,
#'   `"LogLoss"`, `"Accuracy"`, `"Balanced_Accuracy"`, `"Sensitivity"`,
#'   `"Specificity"`, `"Precision"`, `"F1"`, and `"MCC"`.
#' @param type Character. Plot type selected from `"auto"`,
#'   `"observed_predicted"`, `"roc"`, `"threshold"`, or `"confusion"`.
#'   `"auto"` uses an ROC curve for binomial outcomes and an
#'   observed-versus-predicted plot otherwise.
#' @param threshold Optional probability threshold strictly between 0 and 1.
#'   If `NULL`, `attr(object, "threshold")` is used when available; otherwise
#'   `0.5` is used. It is used by `"threshold"` and `"confusion"` plots and
#'   identifies the selected operating point on ROC curves.
#' @param x Character scalar naming the observed-value column in prediction
#'   data. Default is `"observed"`.
#' @param y Character scalar naming the predicted-value column in prediction
#'   data. For binomial outcomes, this column must contain probabilities in
#'   `[0, 1]`. Default is `"predicted"`.
#' @param scale_factor Finite non-zero numeric scalar multiplying both axes
#'   only for `type = "observed_predicted"`. Use `100` to display proportions
#'   as percentages.
#' @param ncol Optional positive integer passed to `ggplot2::facet_wrap()`.
#' @param metric_digits Non-negative integer controlling metric-label
#'   precision.
#' @param x_lab,y_lab Optional axis labels. Plot-type-specific labels are used
#'   when either argument is `NULL`.
#' @param show_fit Logical. For observed-versus-predicted plots, add a linear
#'   fitted line.
#' @param fit_se Logical. Display an uncertainty band around the fitted line.
#' @param fit_color Color of the fitted regression line.
#' @param show_identity Logical. For observed-versus-predicted plots, add the
#'   1:1 line. For ROC curves, add the no-discrimination diagonal.
#' @param identity_color Color of the identity or no-discrimination line.
#' @param identity_linetype Line type of the identity or no-discrimination
#'   line.
#' @param label_x,label_y Numeric positions for metric labels. Infinite
#'   defaults place labels close to a panel corner where supported.
#' @param label_hjust,label_vjust Horizontal and vertical justification for
#'   metric labels.
#' @param label_size Numeric text size for metric labels and confusion-matrix
#'   counts.
#' @param base_size Base font size passed to `ggplot2::theme_bw()`.
#' @param point_args Optional named list of additional arguments passed to
#'   `ggplot2::geom_point()` for observed-versus-predicted plots.
#' @param ... Additional named arguments passed to `ggplot2::geom_point()` for
#'   observed-versus-predicted plots.
#'
#' @return A `ggplot` object that can be extended with additional ggplot2
#'   layers, scales, labels, or themes.
#'
#' @details
#' `plot_performance()` does not refit models or recalculate their stored
#' validation metrics. Predictions are obtained from:
#'
#' `attr(object, "predictions_by_model")`
#'
#' or, when that list is unavailable:
#'
#' `attr(object, "predictions")`
#'
#' For binomial outcomes, ROC curves are reconstructed directly from the
#' observed classes and out-of-fold predicted probabilities. Sensitivity and
#' specificity are calculated over the distinct probability thresholds
#' represented in the selected predictions.
#'
#' The ROC curve uses false-positive rate (`1 - specificity`) on the x-axis
#' and sensitivity on the y-axis. A highlighted point identifies performance
#' at `threshold`.
#'
#' Threshold plots display sensitivity and specificity across probability
#' thresholds and add a vertical line at the selected threshold.
#'
#' Confusion matrices classify probabilities using
#' `predicted >= threshold` and report classification counts.
#'
#' The ROC and threshold plots require both observed classes. The observed
#' values must be coded as `0` and `1`, and predicted values must be finite
#' probabilities in `[0, 1]`.
#'
#' @examples
#' performance_results <- data.frame(
#'   model_id = c("model_1", "model_2"),
#'   CCC = c(0.91, 0.86),
#'   Cb = c(0.96, 0.93),
#'   rho = c(0.95, 0.92),
#'   RMSE = c(0.07, 0.10),
#'   MAE = c(0.05, 0.08)
#' )
#'
#' continuous_predictions <- list(
#'   model_1 = data.frame(
#'     observed = c(0.05, 0.12, 0.20, 0.31, 0.43, 0.55),
#'     predicted = c(0.06, 0.14, 0.18, 0.29, 0.46, 0.53)
#'   ),
#'   model_2 = data.frame(
#'     observed = c(0.05, 0.12, 0.20, 0.31, 0.43, 0.55),
#'     predicted = c(0.09, 0.16, 0.24, 0.27, 0.39, 0.49)
#'   )
#' )
#'
#' attr(
#'   performance_results,
#'   "predictions_by_model"
#' ) <- continuous_predictions
#'
#' continuous_plot <- plot_performance(
#'   performance_results,
#'   metrics = c("CCC", "Cb", "rho"),
#'   scale_factor = 100,
#'   x_lab = "Observed (%)",
#'   y_lab = "Predicted (%)",
#'   size = 2.5,
#'   alpha = 0.8
#' )
#'
#' binary_results <- data.frame(
#'   model_id = c("model_1", "model_2"),
#'   ROC_AUC = c(0.88, 0.82),
#'   Brier = c(0.14, 0.18)
#' )
#'
#' binary_predictions <- list(
#'   model_1 = data.frame(
#'     observed = c(0, 0, 0, 1, 1, 1),
#'     predicted = c(0.08, 0.21, 0.39, 0.58, 0.77, 0.91)
#'   ),
#'   model_2 = data.frame(
#'     observed = c(0, 0, 0, 1, 1, 1),
#'     predicted = c(0.18, 0.31, 0.48, 0.52, 0.68, 0.79)
#'   )
#' )
#'
#' attr(binary_results, "outcome_type") <- "binary"
#' attr(binary_results, "threshold") <- 0.5
#' attr(binary_results, "predictions_by_model") <- binary_predictions
#'
#' roc_plot <- plot_performance(
#'   binary_results,
#'   metrics = c("ROC_AUC", "Brier")
#' )
#'
#' threshold_plot <- plot_performance(
#'   binary_results,
#'   model_id = "model_1",
#'   metrics = "ROC_AUC",
#'   type = "threshold"
#' )
#'
#' confusion_plot <- plot_performance(
#'   binary_results,
#'   model_id = "model_1",
#'   metrics = c("ROC_AUC", "Brier"),
#'   type = "confusion"
#' )
#'
#' stopifnot(
#'   inherits(continuous_plot, "ggplot"),
#'   inherits(roc_plot, "ggplot"),
#'   inherits(threshold_plot, "ggplot"),
#'   inherits(confusion_plot, "ggplot")
#' )
#'
#' @export
plot_performance <- function(
    object,
    model_id = NULL,
    metrics = "CCC",
    type = c(
      "auto",
      "observed_predicted",
      "roc",
      "threshold",
      "confusion"
    ),
    threshold = NULL,
    x = "observed",
    y = "predicted",
    scale_factor = 1,
    ncol = NULL,
    metric_digits = 2L,
    x_lab = NULL,
    y_lab = NULL,
    show_fit = TRUE,
    fit_se = FALSE,
    fit_color = "orange",
    show_identity = TRUE,
    identity_color = "red",
    identity_linetype = "dashed",
    label_x = -Inf,
    label_y = Inf,
    label_hjust = -0.1,
    label_vjust = 1.2,
    label_size = 3.5,
    base_size = 11,
    point_args = list(),
    ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required by `plot_performance()`.",
      call. = FALSE
    )
  }

  type <- match.arg(type)

  scalar_name <- function(value) {
    is.character(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      nzchar(value)
  }

  scalar_flag <- function(value) {
    is.logical(value) &&
      length(value) == 1L &&
      !is.na(value)
  }

  whole_scalar <- function(value) {
    is.numeric(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      is.finite(value) &&
      value == as.integer(value)
  }

  scalar_color <- function(value, argument) {
    if (!is.character(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !nzchar(value)) {
      stop(
        "`",
        argument,
        "` must be one non-empty color value.",
        call. = FALSE
      )
    }
  }

  if (!is.data.frame(object) || !nrow(object)) {
    stop(
      "`object` must be a non-empty data frame returned by ",
      "`find_bestfit()`.",
      call. = FALSE
    )
  }

  if (!"model_id" %in% names(object)) {
    stop(
      "`object` must contain a `model_id` column.",
      call. = FALSE
    )
  }

  if (anyNA(object$model_id) ||
      anyDuplicated(as.character(object$model_id))) {
    stop(
      "`object$model_id` must contain unique non-missing identifiers.",
      call. = FALSE
    )
  }

  if (!is.character(metrics) ||
      !length(metrics) ||
      anyNA(metrics) ||
      any(!nzchar(metrics)) ||
      anyDuplicated(metrics)) {
    stop(
      "`metrics` must contain unique, non-empty metric names.",
      call. = FALSE
    )
  }

  missing_metrics <- setdiff(
    metrics,
    names(object)
  )

  if (length(missing_metrics)) {
    stop(
      "Metric column(s) not found in `object`: ",
      paste(missing_metrics, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  nonnumeric_metrics <- metrics[
    !vapply(
      object[metrics],
      is.numeric,
      logical(1)
    )
  ]

  if (length(nonnumeric_metrics)) {
    stop(
      "Metric column(s) must be numeric: ",
      paste(nonnumeric_metrics, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!scalar_name(x)) {
    stop(
      "`x` must be one non-empty column name.",
      call. = FALSE
    )
  }

  if (!scalar_name(y)) {
    stop(
      "`y` must be one non-empty column name.",
      call. = FALSE
    )
  }

  if (!is.numeric(scale_factor) ||
      length(scale_factor) != 1L ||
      is.na(scale_factor) ||
      !is.finite(scale_factor) ||
      scale_factor == 0) {
    stop(
      "`scale_factor` must be one finite non-zero numeric value.",
      call. = FALSE
    )
  }

  if (!is.null(ncol) &&
      (
        !whole_scalar(ncol) ||
        ncol < 1L
      )) {
    stop(
      "`ncol` must be NULL or a positive integer.",
      call. = FALSE
    )
  }

  if (!whole_scalar(metric_digits) ||
      metric_digits < 0L) {
    stop(
      "`metric_digits` must be a non-negative integer.",
      call. = FALSE
    )
  }

  metric_digits <- as.integer(metric_digits)

  if (!scalar_flag(show_fit)) {
    stop(
      "`show_fit` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!scalar_flag(fit_se)) {
    stop(
      "`fit_se` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!scalar_flag(show_identity)) {
    stop(
      "`show_identity` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  scalar_color(
    fit_color,
    "fit_color"
  )

  scalar_color(
    identity_color,
    "identity_color"
  )

  if (!is.list(point_args)) {
    stop(
      "`point_args` must be a list.",
      call. = FALSE
    )
  }

  if (length(point_args) &&
      (
        is.null(names(point_args)) ||
        anyNA(names(point_args)) ||
        any(!nzchar(names(point_args)))
      )) {
    stop(
      "All entries in `point_args` must be named.",
      call. = FALSE
    )
  }

  dots <- list(...)

  if (length(dots) &&
      (
        is.null(names(dots)) ||
        anyNA(names(dots)) ||
        any(!nzchar(names(dots)))
      )) {
    stop(
      "All arguments supplied through `...` must be named.",
      call. = FALSE
    )
  }

  duplicated_point_args <- intersect(
    names(point_args),
    names(dots)
  )

  if (length(duplicated_point_args)) {
    stop(
      "Point argument(s) supplied in both `point_args` and `...`: ",
      paste(duplicated_point_args, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  point_args <- c(
    point_args,
    dots
  )

  predictions_by_model <- attr(
    object,
    "predictions_by_model",
    exact = TRUE
  )

  if (!is.list(predictions_by_model) ||
      !length(predictions_by_model)) {
    prediction_data <- attr(
      object,
      "predictions",
      exact = TRUE
    )

    if (!is.data.frame(prediction_data) ||
        !nrow(prediction_data) ||
        !"model_id" %in% names(prediction_data)) {
      stop(
        "`object` does not contain usable `predictions_by_model` or ",
        "`predictions` attributes. Use the original object returned by ",
        "`find_bestfit()` rather than an exported/imported table.",
        call. = FALSE
      )
    }

    predictions_by_model <- split(
      prediction_data,
      as.character(prediction_data$model_id),
      drop = TRUE
    )
  }

  available_prediction_ids <- names(
    predictions_by_model
  )

  if (is.null(available_prediction_ids) ||
      anyNA(available_prediction_ids) ||
      any(!nzchar(available_prediction_ids)) ||
      anyDuplicated(available_prediction_ids)) {
    stop(
      "`predictions_by_model` must have unique non-empty model-ID names.",
      call. = FALSE
    )
  }

  object_ids <- as.character(
    object$model_id
  )

  if (is.null(model_id)) {
    selected_ids <- object_ids[
      object_ids %in% available_prediction_ids
    ]
  } else {
    if (!length(model_id) ||
        anyNA(model_id) ||
        anyDuplicated(as.character(model_id))) {
      stop(
        "`model_id` must contain unique, non-missing identifiers.",
        call. = FALSE
      )
    }

    selected_ids <- as.character(
      model_id
    )

    missing_in_object <- setdiff(
      selected_ids,
      object_ids
    )

    if (length(missing_in_object)) {
      stop(
        "Selected model_id value(s) not found in `object`: ",
        paste(missing_in_object, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    missing_predictions <- setdiff(
      selected_ids,
      available_prediction_ids
    )

    if (length(missing_predictions)) {
      stop(
        "Prediction data are unavailable for model_id value(s): ",
        paste(missing_predictions, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (!length(selected_ids)) {
    stop(
      "No selected models have prediction data.",
      call. = FALSE
    )
  }

  prediction_rows <- lapply(
    selected_ids,
    function(id) {
      current <- predictions_by_model[[id]]

      if (!is.data.frame(current) ||
          !nrow(current)) {
        stop(
          "Prediction data for model_id ",
          id,
          " are empty or invalid.",
          call. = FALSE
        )
      }

      missing_xy <- setdiff(
        c(x, y),
        names(current)
      )

      if (length(missing_xy)) {
        stop(
          "Prediction data for model_id ",
          id,
          " lack column(s): ",
          paste(missing_xy, collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      if (!is.numeric(current[[x]]) ||
          !is.numeric(current[[y]])) {
        stop(
          "The selected observed and predicted columns must be numeric.",
          call. = FALSE
        )
      }

      if (anyNA(current[[x]]) ||
          anyNA(current[[y]]) ||
          any(!is.finite(current[[x]])) ||
          any(!is.finite(current[[y]]))) {
        stop(
          "Observed and predicted values must be finite and non-missing.",
          call. = FALSE
        )
      }

      current$model_id <- id
      current
    }
  )

  prediction_data <- do.call(
    rbind,
    prediction_rows
  )

  rownames(prediction_data) <- NULL

  prediction_data$model_id <- factor(
    as.character(prediction_data$model_id),
    levels = selected_ids
  )

  outcome_type <- attr(
    object,
    "outcome_type",
    exact = TRUE
  )

  if (!is.null(outcome_type)) {
    outcome_type <- as.character(
      outcome_type
    )
  }

  observed_values <- prediction_data[[x]]
  predicted_values <- prediction_data[[y]]

  binary_from_data <-
    all(observed_values %in% c(0, 1)) &&
    all(predicted_values >= 0) &&
    all(predicted_values <= 1)

  binary_outcome <-
    identical(outcome_type, "binary") ||
    binary_from_data

  if (identical(type, "auto")) {
    type <- if (binary_outcome) {
      "roc"
    } else {
      "observed_predicted"
    }
  }

  binary_types <- c(
    "roc",
    "threshold",
    "confusion"
  )

  if (type %in% binary_types &&
      !binary_outcome) {
    stop(
      "`type = ",
      sQuote(type),
      "` requires a binomial outcome with observed classes coded 0/1 ",
      "and predicted probabilities in [0, 1].",
      call. = FALSE
    )
  }

  if (binary_outcome) {
    if (!all(observed_values %in% c(0, 1))) {
      stop(
        "Binomial performance plots require observed outcomes coded 0/1.",
        call. = FALSE
      )
    }

    if (any(predicted_values < 0 |
            predicted_values > 1)) {
      stop(
        "Binomial performance plots require predicted probabilities in ",
        "[0, 1].",
        call. = FALSE
      )
    }
  }

  if (type %in% c("roc", "threshold") &&
      length(unique(observed_values)) < 2L) {
    stop(
      "ROC and threshold plots require both observed classes 0 and 1.",
      call. = FALSE
    )
  }

  if (type %in% c("roc", "threshold") &&
      all(predicted_values %in% c(0, 1))) {
    warning(
      "Predicted values contain only 0 and 1. ROC and threshold curves ",
      "will have limited probability-resolution information.",
      call. = FALSE
    )
  }

  stored_threshold <- attr(
    object,
    "threshold",
    exact = TRUE
  )

  if (is.null(threshold)) {
    threshold <- if (
      is.numeric(stored_threshold) &&
      length(stored_threshold) == 1L &&
      !is.na(stored_threshold) &&
      is.finite(stored_threshold) &&
      stored_threshold > 0 &&
      stored_threshold < 1
    ) {
      as.numeric(stored_threshold)
    } else {
      0.5
    }
  }

  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      is.na(threshold) ||
      !is.finite(threshold) ||
      threshold <= 0 ||
      threshold >= 1) {
    stop(
      "`threshold` must be one finite probability strictly between 0 and 1.",
      call. = FALSE
    )
  }

  metric_rows <- object[
    match(selected_ids, object_ids),
    ,
    drop = FALSE
  ]

  metric_rows$model_id <- factor(
    as.character(metric_rows$model_id),
    levels = selected_ids
  )

  format_metric <- function(
    metric_name,
    value
  ) {
    formatted_value <- if (is.na(value)) {
      "NA"
    } else {
      formatC(
        value,
        format = "f",
        digits = metric_digits
      )
    }

    display_name <- gsub(
      "_",
      " ",
      metric_name,
      fixed = TRUE
    )

    paste0(
      display_name,
      " = ",
      formatted_value
    )
  }

  metric_rows$.metric_label <- vapply(
    seq_len(nrow(metric_rows)),
    function(i) {
      paste(
        vapply(
          metrics,
          function(metric_name) {
            format_metric(
              metric_name,
              metric_rows[[metric_name]][i]
            )
          },
          character(1)
        ),
        collapse = "\n"
      )
    },
    character(1)
  )

  compute_binary_curve <- function(
    observed,
    predicted,
    model_identifier
  ) {
    thresholds <- sort(
      unique(
        c(
          1,
          predicted,
          0
        )
      ),
      decreasing = TRUE
    )

    rows <- lapply(
      thresholds,
      function(current_threshold) {
        predicted_class <- as.integer(
          predicted >= current_threshold
        )

        tp <- sum(
          observed == 1 &
            predicted_class == 1
        )

        fn <- sum(
          observed == 1 &
            predicted_class == 0
        )

        tn <- sum(
          observed == 0 &
            predicted_class == 0
        )

        fp <- sum(
          observed == 0 &
            predicted_class == 1
        )

        sensitivity <- if (
          tp + fn > 0
        ) {
          tp / (tp + fn)
        } else {
          NA_real_
        }

        specificity <- if (
          tn + fp > 0
        ) {
          tn / (tn + fp)
        } else {
          NA_real_
        }

        data.frame(
          model_id = model_identifier,
          threshold = current_threshold,
          sensitivity = sensitivity,
          specificity = specificity,
          false_positive_rate =
            1 - specificity,
          stringsAsFactors = FALSE
        )
      }
    )

    output <- do.call(
      rbind,
      rows
    )

    rownames(output) <- NULL
    output
  }

  add_metric_labels <- function(
    plot_object,
    x_position,
    y_position,
    hjust_value = label_hjust,
    vjust_value = label_vjust
  ) {
    plot_object +
      ggplot2::geom_text(
        data = metric_rows,
        mapping = ggplot2::aes(
          label = .data$.metric_label
        ),
        x = x_position,
        y = y_position,
        hjust = hjust_value,
        vjust = vjust_value,
        size = label_size,
        inherit.aes = FALSE
      )
  }

  if (identical(
    type,
    "observed_predicted"
  )) {
    prediction_data$.plot_x <-
      prediction_data[[x]] *
      scale_factor

    prediction_data$.plot_y <-
      prediction_data[[y]] *
      scale_factor

    if (is.null(x_lab)) {
      x_lab <- tools::toTitleCase(
        gsub(
          "_",
          " ",
          x
        )
      )
    }

    if (is.null(y_lab)) {
      y_lab <- tools::toTitleCase(
        gsub(
          "_",
          " ",
          y
        )
      )
    }

    plot_object <- ggplot2::ggplot(
      prediction_data,
      ggplot2::aes(
        x = .data$.plot_x,
        y = .data$.plot_y
      )
    )

    point_call <- c(
      list(
        mapping = NULL,
        data = NULL,
        inherit.aes = TRUE
      ),
      point_args
    )

    plot_object <- plot_object +
      do.call(
        ggplot2::geom_point,
        point_call
      )

    if (show_fit) {
      plot_object <- plot_object +
        ggplot2::geom_smooth(
          method = "lm",
          formula = y ~ x,
          se = fit_se,
          colour = fit_color
        )
    }

    if (show_identity) {
      plot_object <- plot_object +
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          colour = identity_color,
          linetype = identity_linetype
        )
    }

    plot_object <- add_metric_labels(
      plot_object,
      label_x,
      label_y
    )

    return(
      plot_object +
        ggplot2::facet_wrap(
          ggplot2::vars(
            .data$model_id
          ),
          ncol = ncol
        ) +
        ggplot2::theme_bw(
          base_size = base_size
        ) +
        ggplot2::labs(
          x = x_lab,
          y = y_lab
        ) +
        ggplot2::theme(
          legend.position = "none"
        )
    )
  }

  curve_rows <- lapply(
    selected_ids,
    function(id) {
      current <- prediction_data[
        as.character(
          prediction_data$model_id
        ) == id,
        ,
        drop = FALSE
      ]

      compute_binary_curve(
        observed = current[[x]],
        predicted = current[[y]],
        model_identifier = id
      )
    }
  )

  curve_data <- do.call(
    rbind,
    curve_rows
  )

  rownames(curve_data) <- NULL

  curve_data$model_id <- factor(
    as.character(curve_data$model_id),
    levels = selected_ids
  )

  operating_rows <- lapply(
    selected_ids,
    function(id) {
      current <- prediction_data[
        as.character(
          prediction_data$model_id
        ) == id,
        ,
        drop = FALSE
      ]

      observed <- as.integer(
        current[[x]]
      )

      predicted_class <- as.integer(
        current[[y]] >= threshold
      )

      tp <- sum(
        observed == 1 &
          predicted_class == 1
      )

      fn <- sum(
        observed == 1 &
          predicted_class == 0
      )

      tn <- sum(
        observed == 0 &
          predicted_class == 0
      )

      fp <- sum(
        observed == 0 &
          predicted_class == 1
      )

      sensitivity <- if (
        tp + fn > 0
      ) {
        tp / (tp + fn)
      } else {
        NA_real_
      }

      specificity <- if (
        tn + fp > 0
      ) {
        tn / (tn + fp)
      } else {
        NA_real_
      }

      data.frame(
        model_id = id,
        threshold = threshold,
        sensitivity = sensitivity,
        specificity = specificity,
        false_positive_rate = 1 - specificity,
        stringsAsFactors = FALSE
      )
    }
  )

  operating_data <- do.call(
    rbind,
    operating_rows
  )

  rownames(operating_data) <- NULL

  operating_data$model_id <- factor(
    as.character(operating_data$model_id),
    levels = selected_ids
  )

  if (identical(type, "roc")) {
    if (is.null(x_lab)) {
      x_lab <- "False positive rate (1 - specificity)"
    }

    if (is.null(y_lab)) {
      y_lab <- "Sensitivity"
    }

    plot_object <- ggplot2::ggplot(
      curve_data,
      ggplot2::aes(
        x = .data$false_positive_rate,
        y = .data$sensitivity
      )
    ) +
      ggplot2::geom_path(
        linewidth = 0.8,
        colour = "#2C7FB8"
      ) +
      ggplot2::geom_point(
        data = operating_data,
        mapping = ggplot2::aes(
          x = .data$false_positive_rate,
          y = .data$sensitivity
        ),
        inherit.aes = FALSE,
        size = 2.5,
        colour = "#D95F0E"
      )

    if (show_identity) {
      plot_object <- plot_object +
        ggplot2::geom_abline(
          slope = 1,
          intercept = 0,
          colour = identity_color,
          linetype = identity_linetype
        )
    }

    plot_object <- add_metric_labels(
      plot_object,
      -Inf,
      Inf,
      hjust_value = -0.1,
      vjust_value = 1.2
    )

    return(
      plot_object +
        ggplot2::facet_wrap(
          ggplot2::vars(
            .data$model_id
          ),
          ncol = ncol
        ) +
        ggplot2::coord_equal(
          xlim = c(0, 1),
          ylim = c(0, 1),
          expand = FALSE
        ) +
        ggplot2::scale_x_continuous(
          breaks = seq(
            0,
            1,
            by = 0.2
          )
        ) +
        ggplot2::scale_y_continuous(
          breaks = seq(
            0,
            1,
            by = 0.2
          )
        ) +
        ggplot2::theme_bw(
          base_size = base_size
        ) +
        ggplot2::labs(
          x = x_lab,
          y = y_lab
        ) +
        ggplot2::theme(
          legend.position = "none",
          panel.grid.minor = ggplot2::element_blank()
        )
    )
  }

  if (identical(type, "threshold")) {
    if (is.null(x_lab)) {
      x_lab <- "Probability threshold"
    }

    if (is.null(y_lab)) {
      y_lab <- "Metric value"
    }

    threshold_long <- rbind(
      data.frame(
        model_id = curve_data$model_id,
        threshold = curve_data$threshold,
        metric = "Sensitivity",
        value = curve_data$sensitivity,
        stringsAsFactors = FALSE
      ),
      data.frame(
        model_id = curve_data$model_id,
        threshold = curve_data$threshold,
        metric = "Specificity",
        value = curve_data$specificity,
        stringsAsFactors = FALSE
      )
    )

    threshold_long$model_id <- factor(
      as.character(
        threshold_long$model_id
      ),
      levels = selected_ids
    )

    threshold_long$metric <- factor(
      threshold_long$metric,
      levels = c(
        "Sensitivity",
        "Specificity"
      )
    )

    return(
      ggplot2::ggplot(
        threshold_long,
        ggplot2::aes(
          x = .data$threshold,
          y = .data$value,
          colour = .data$metric
        )
      ) +
        ggplot2::geom_line(
          linewidth = 0.8
        ) +
        ggplot2::geom_vline(
          xintercept = threshold,
          colour = "#D95F0E",
          linetype = "dashed"
        ) +
        ggplot2::facet_wrap(
          ggplot2::vars(
            .data$model_id
          ),
          ncol = ncol
        ) +
        ggplot2::scale_x_continuous(
          limits = c(0, 1),
          breaks = seq(
            0,
            1,
            by = 0.2
          )
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(
            0,
            1,
            by = 0.2
          )
        ) +
        ggplot2::scale_colour_manual(
          values = c(
            Sensitivity = "#2C7FB8",
            Specificity = "#D95F0E"
          )
        ) +
        ggplot2::theme_bw(
          base_size = base_size
        ) +
        ggplot2::labs(
          x = x_lab,
          y = y_lab,
          colour = NULL
        ) +
        ggplot2::theme(
          legend.position = "bottom",
          panel.grid.minor = ggplot2::element_blank()
        )
    )
  }

  confusion_rows <- lapply(
    selected_ids,
    function(id) {
      current <- prediction_data[
        as.character(
          prediction_data$model_id
        ) == id,
        ,
        drop = FALSE
      ]

      observed_class <- as.integer(
        current[[x]]
      )

      predicted_class <- as.integer(
        current[[y]] >= threshold
      )

      grid <- expand.grid(
        observed_class = c(0L, 1L),
        predicted_class = c(0L, 1L),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      grid$count <- vapply(
        seq_len(nrow(grid)),
        function(i) {
          sum(
            observed_class ==
              grid$observed_class[i] &
              predicted_class ==
              grid$predicted_class[i]
          )
        },
        integer(1)
      )

      grid$model_id <- id
      grid
    }
  )

  confusion_data <- do.call(
    rbind,
    confusion_rows
  )

  rownames(confusion_data) <- NULL

  confusion_data$model_id <- factor(
    as.character(
      confusion_data$model_id
    ),
    levels = selected_ids
  )

  confusion_data$observed_class <- factor(
    confusion_data$observed_class,
    levels = c(0, 1),
    labels = c(
      "0",
      "1"
    )
  )

  confusion_data$predicted_class <- factor(
    confusion_data$predicted_class,
    levels = c(0, 1),
    labels = c(
      "0",
      "1"
    )
  )

  if (is.null(x_lab)) {
    x_lab <- paste0(
      "Predicted class (threshold = ",
      formatC(
        threshold,
        format = "f",
        digits = metric_digits
      ),
      ")"
    )
  }

  if (is.null(y_lab)) {
    y_lab <- "Observed class"
  }

  ggplot2::ggplot(
    confusion_data,
    ggplot2::aes(
      x = .data$predicted_class,
      y = .data$observed_class,
      fill = .data$count
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = .data$count
      ),
      size = label_size
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(
        .data$model_id
      ),
      ncol = ncol
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#2C7FB8"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(
      base_size = base_size
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      fill = "Count"
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}
