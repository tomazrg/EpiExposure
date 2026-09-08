#' Plot observed versus predicted model performance
#'
#' Creates observed-versus-predicted performance plots for models evaluated by
#' `find_bestfit()`. The function retrieves out-of-fold predictions from the
#' `"predictions_by_model"` attribute and model-level validation metrics from
#' the main `find_bestfit()` result.
#'
#' One panel is produced per selected model. Validation metrics are displayed
#' inside each panel in the order supplied through `metrics`. For example,
#' `metrics = c("CCC", "Cb", "rho")` displays CCC first, followed by Cb and
#' rho on separate lines.
#'
#' @param object A data frame returned by `find_bestfit()`. The object must
#'   contain a `model_id` column and the requested validation-metric columns,
#'   and must retain the `"predictions_by_model"` attribute.
#' @param model_id Optional vector of model identifiers to plot. If `NULL`, all
#'   models retained in `object` that have prediction data are plotted. Models
#'   are displayed in the order supplied by `model_id`, or in the current order
#'   of `object` when `model_id = NULL`.
#' @param metrics Character vector naming validation metrics stored as columns
#'   in `object`, such as `"CCC"`, `"Cb"`, `"rho"`, `"RMSE"`, `"MAE"`,
#'   `"ROC_AUC"`, or `"Brier"`. Metrics are printed in the supplied order.
#' @param x Character scalar naming the observed-value column in the prediction
#'   data. Default is `"observed"`.
#' @param y Character scalar naming the predicted-value column in the prediction
#'   data. Default is `"predicted"`.
#' @param scale_factor Finite numeric scalar multiplying both plotted axes.
#'   Use `100` to display proportions as percentages and `1` to retain the
#'   original response scale.
#' @param ncol Optional positive integer passed to `ggplot2::facet_wrap()`.
#' @param metric_digits Non-negative integer number of decimal places used in
#'   metric labels.
#' @param x_lab Optional x-axis label. If `NULL`, the name supplied to `x` is
#'   converted to title case.
#' @param y_lab Optional y-axis label. If `NULL`, the name supplied to `y` is
#'   converted to title case.
#' @param show_fit Logical. If `TRUE`, adds a fitted regression line using
#'   `ggplot2::geom_smooth(method = "lm")`.
#' @param fit_se Logical. Display the uncertainty band around the fitted line.
#' @param fit_color Color of the fitted regression line.
#' @param show_identity Logical. If `TRUE`, adds the 1:1 identity line.
#' @param identity_color Color of the identity line.
#' @param identity_linetype Line type of the identity line.
#' @param label_x,label_y Numeric positions for metric labels. Defaults are
#'   `-Inf` and `Inf`, placing labels at the upper-left corner of each panel.
#' @param label_hjust,label_vjust Horizontal and vertical justification for
#'   metric labels.
#' @param label_size Numeric text size for metric labels.
#' @param base_size Base font size passed to `ggplot2::theme_bw()`.
#' @param point_args Optional named list of additional arguments passed to
#'   `ggplot2::geom_point()`. This is useful when point customization should be
#'   supplied programmatically.
#' @param ... Additional arguments passed directly to
#'   `ggplot2::geom_point()`, such as `size`, `shape`, `colour`, `alpha`,
#'   `fill`, or `stroke`.
#'
#' @return A `ggplot` object. The plot can be extended with additional ggplot2
#'   layers, scales, labels, or themes using `+`.
#'
#' @details
#' `plot_performance()` does not recalculate predictions or performance
#' metrics. Observed and predicted values are obtained from:
#'
#' `attr(object, "predictions_by_model")`
#'
#' The metric annotations are obtained from columns of `object`. Consequently,
#' the function preserves the exact out-of-fold values and validation metrics
#' calculated by `find_bestfit()`.
#'
#' The default axes use the original response scale. For Beta or Binomial
#' outcomes stored as proportions, use `scale_factor = 100`, together with axis
#' labels such as `x_lab = "Observed (%)"` and `y_lab = "Predicted (%)"`.
#'
#' @examples
#' \dontrun{
#' plot_performance(
#'   best_models,
#'   model_id = best_models$model_id[1:6],
#'   metrics = c("CCC", "Cb", "rho"),
#'   scale_factor = 100,
#'   ncol = 3,
#'   x_lab = "Observed (%)",
#'   y_lab = "Predicted (%)",
#'   size = 2.5,
#'   alpha = 0.8
#' )
#'
#' plot_performance(
#'   best_models,
#'   metrics = c("RMSE", "MAE"),
#'   point_args = list(shape = 21, fill = "steelblue", colour = "black")
#' ) +
#'   ggplot2::theme(text = ggplot2::element_text(face = "bold"))
#' }
#'
#' @export
plot_performance <- function(
    object,
    model_id = NULL,
    metrics = "CCC",
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
    stop("Package 'ggplot2' is required by `plot_performance()`.", call. = FALSE)
  }
  
  scalar_name <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
  }
  
  scalar_flag <- function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
  }
  
  whole_scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == as.integer(value)
  }
  
  if (!is.data.frame(object) || !nrow(object)) {
    stop("`object` must be a non-empty data frame returned by `find_bestfit()`.",
         call. = FALSE)
  }
  
  if (!"model_id" %in% names(object)) {
    stop("`object` must contain a `model_id` column.", call. = FALSE)
  }
  
  if (anyNA(object$model_id) || anyDuplicated(object$model_id)) {
    stop("`object$model_id` must contain unique non-missing identifiers.",
         call. = FALSE)
  }
  
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      any(!nzchar(metrics)) || anyDuplicated(metrics)) {
    stop("`metrics` must contain unique, non-empty metric names.", call. = FALSE)
  }
  
  missing_metrics <- setdiff(metrics, names(object))
  if (length(missing_metrics)) {
    stop(
      "Metric column(s) not found in `object`: ",
      paste(missing_metrics, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  nonnumeric_metrics <- metrics[!vapply(object[metrics], is.numeric, logical(1))]
  if (length(nonnumeric_metrics)) {
    stop(
      "Metric column(s) must be numeric: ",
      paste(nonnumeric_metrics, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  if (!scalar_name(x)) stop("`x` must be one non-empty column name.", call. = FALSE)
  if (!scalar_name(y)) stop("`y` must be one non-empty column name.", call. = FALSE)
  
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      is.na(scale_factor) || !is.finite(scale_factor)) {
    stop("`scale_factor` must be one finite numeric value.", call. = FALSE)
  }
  
  if (!is.null(ncol) && (!whole_scalar(ncol) || ncol < 1L)) {
    stop("`ncol` must be NULL or a positive integer.", call. = FALSE)
  }
  
  if (!whole_scalar(metric_digits) || metric_digits < 0L) {
    stop("`metric_digits` must be a non-negative integer.", call. = FALSE)
  }
  metric_digits <- as.integer(metric_digits)
  
  if (!scalar_flag(show_fit)) stop("`show_fit` must be TRUE or FALSE.", call. = FALSE)
  if (!scalar_flag(fit_se)) stop("`fit_se` must be TRUE or FALSE.", call. = FALSE)
  if (!scalar_flag(show_identity)) {
    stop("`show_identity` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.list(point_args)) {
    stop("`point_args` must be a list.", call. = FALSE)
  }
  
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("All arguments supplied through `...` must be named.", call. = FALSE)
  }
  
  duplicated_point_args <- intersect(names(point_args), names(dots))
  if (length(duplicated_point_args)) {
    stop(
      "Point argument(s) supplied in both `point_args` and `...`: ",
      paste(duplicated_point_args, collapse = ", "), ".",
      call. = FALSE
    )
  }
  point_args <- c(point_args, dots)
  
  predictions_by_model <- attr(object, "predictions_by_model", exact = TRUE)
  
  if (!is.list(predictions_by_model) || !length(predictions_by_model)) {
    prediction_data <- attr(object, "predictions", exact = TRUE)
    
    if (!is.data.frame(prediction_data) || !nrow(prediction_data) ||
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
  
  available_prediction_ids <- names(predictions_by_model)
  if (is.null(available_prediction_ids) || any(!nzchar(available_prediction_ids))) {
    stop("`predictions_by_model` must be named by `model_id`.", call. = FALSE)
  }
  
  object_ids <- as.character(object$model_id)
  
  if (is.null(model_id)) {
    selected_ids <- object_ids[object_ids %in% available_prediction_ids]
  } else {
    if (!length(model_id) || anyNA(model_id) || anyDuplicated(model_id)) {
      stop("`model_id` must contain unique, non-missing identifiers.",
           call. = FALSE)
    }
    selected_ids <- as.character(model_id)
    
    missing_in_object <- setdiff(selected_ids, object_ids)
    if (length(missing_in_object)) {
      stop(
        "Selected model_id value(s) not found in `object`: ",
        paste(missing_in_object, collapse = ", "), ".",
        call. = FALSE
      )
    }
    
    missing_predictions <- setdiff(selected_ids, available_prediction_ids)
    if (length(missing_predictions)) {
      stop(
        "Prediction data are unavailable for model_id value(s): ",
        paste(missing_predictions, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  
  if (!length(selected_ids)) {
    stop("No selected models have prediction data.", call. = FALSE)
  }
  
  prediction_rows <- lapply(selected_ids, function(id) {
    current <- predictions_by_model[[id]]
    if (!is.data.frame(current) || !nrow(current)) {
      stop("Prediction data for model_id ", id, " are empty or invalid.",
           call. = FALSE)
    }
    
    missing_xy <- setdiff(c(x, y), names(current))
    if (length(missing_xy)) {
      stop(
        "Prediction data for model_id ", id,
        " lack column(s): ", paste(missing_xy, collapse = ", "), ".",
        call. = FALSE
      )
    }
    
    if (!is.numeric(current[[x]]) || !is.numeric(current[[y]])) {
      stop("The selected `x` and `y` prediction columns must be numeric.",
           call. = FALSE)
    }
    
    current$model_id <- id
    current$.plot_x <- current[[x]] * scale_factor
    current$.plot_y <- current[[y]] * scale_factor
    current
  })
  
  prediction_data <- do.call(rbind, prediction_rows)
  rownames(prediction_data) <- NULL
  
  if (any(!is.finite(prediction_data$.plot_x)) ||
      any(!is.finite(prediction_data$.plot_y))) {
    stop("The selected observed and predicted values must be finite.",
         call. = FALSE)
  }
  
  prediction_data$model_id <- factor(
    as.character(prediction_data$model_id),
    levels = selected_ids
  )
  
  metric_rows <- object[match(selected_ids, object_ids), , drop = FALSE]
  metric_rows$model_id <- factor(
    as.character(metric_rows$model_id),
    levels = selected_ids
  )
  
  format_metric <- function(metric_name, value) {
    formatted_value <- if (is.na(value)) {
      "NA"
    } else {
      formatC(value, format = "f", digits = metric_digits)
    }
    paste0(metric_name, " = ", formatted_value)
  }
  
  metric_rows$.metric_label <- vapply(
    seq_len(nrow(metric_rows)),
    function(i) {
      paste(
        vapply(
          metrics,
          function(metric_name) {
            format_metric(metric_name, metric_rows[[metric_name]][i])
          },
          character(1)
        ),
        collapse = "\n"
      )
    },
    character(1)
  )
  
  if (is.null(x_lab)) {
    x_lab <- tools::toTitleCase(gsub("_", " ", x))
  }
  if (is.null(y_lab)) {
    y_lab <- tools::toTitleCase(gsub("_", " ", y))
  }
  
  plot_object <- ggplot2::ggplot(
    prediction_data,
    ggplot2::aes(x = .data$.plot_x, y = .data$.plot_y)
  )
  
  point_call <- c(
    list(mapping = NULL, data = NULL, inherit.aes = TRUE),
    point_args
  )
  plot_object <- plot_object + do.call(ggplot2::geom_point, point_call)
  
  if (show_fit) {
    plot_object <- plot_object + ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = fit_se,
      colour = fit_color
    )
  }
  
  if (show_identity) {
    plot_object <- plot_object + ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      colour = identity_color,
      linetype = identity_linetype
    )
  }
  
  plot_object <- plot_object +
    ggplot2::geom_text(
      data = metric_rows,
      mapping = ggplot2::aes(label = .data$.metric_label),
      x = label_x,
      y = label_y,
      hjust = label_hjust,
      vjust = label_vjust,
      size = label_size,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$model_id),
      ncol = ncol
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab
    ) +
    ggplot2::theme(
      legend.position = "none"
    )
  
  plot_object
}
