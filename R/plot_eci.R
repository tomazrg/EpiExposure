#' Plot ECI and lag-specific contribution results
#'
#' Creates publication-ready plots from objects returned by `compute_eci()` and
#' `compute_ecilag()`. Either object may be supplied alone. When both are
#' supplied, the overall ECI perspective and the lag-specific perspective are
#' combined and labelled "(a)" and "(b)", respectively.
#'
#' The function recognizes draw-level output (`output = "samples"`) through a
#' `sample` column and summary output (`output = "summary"`) through the
#' corresponding estimate and interval columns. For `compute_ecilag()` objects,
#' `by_lag_samples` and `by_lag` are extracted internally.
#'
#' @param eci_overal Optional data frame returned by `compute_eci()`.
#' @param eci_lag Optional object returned by `compute_ecilag()`. A list with
#'   `by_lag_samples` and/or `by_lag` is handled internally. A lag-level data
#'   frame is also accepted.
#' @param facet_scales Character scalar controlling facet scales. One of
#'   `"fixed"`, `"free"`, `"free_x"`, or `"free_y"`. Default is `"free_y"`.
#' @param overall_nrow Positive integer number of facet rows for the overall ECI
#'   panel. Default is `3`.
#' @param lag_ncol Positive integer number of facet columns for lag panels.
#'   Default is `1`.
#' @param sample_smooth Logical scalar. If `TRUE`, draw-level curves are shown
#'   with `geom_smooth(se = FALSE)`, reproducing the historical manual plot.
#'   If `FALSE`, draw-level curves are connected with `geom_line()`.
#' @param smooth_method Optional smoothing method passed to `geom_smooth()`.
#'   The default `NULL` lets ggplot2 select its standard method.
#' @param smooth_formula Optional formula passed to `geom_smooth()`. Default is
#'   `NULL`.
#' @param linewidth Positive finite line width. Default is `0.7`.
#' @param alpha Finite number in `[0, 1]` controlling uncertainty-ribbon alpha.
#'   Default is `0.20`.
#' @param overall_x_label,overall_y_label,lag_x_label,lag_y_label,
#'   exposure_y_label Axis labels.
#' @param exposure_colors Character vector of colors used for exposure-
#'   difference curves. Named vectors are supported. The default reproduces the
#'   historical manual palette.
#' @param base_size Positive finite base font size. Default is `10`.
#' @param panel_labels Character vector of length two used when both
#'   perspectives are combined. Default is `c("(a)", "(b)")`.
#' @param label_size Positive finite panel-label size. Default is `12`.
#' @param rel_widths Optional positive numeric vector controlling relative panel
#'   widths. By default, all displayed panels have equal width. Supply two values
#'   when two panels are displayed or three values when the overall ECI, lag
#'   contribution, and exposure-difference panels are displayed together.
#'
#' @return A ggplot/cowplot object. When both inputs are supplied, the returned
#'   object combines the overall ECI panel and the lag-specific panel. When only
#'   `eci_lag` is supplied, contribution and exposure-difference panels are
#'   combined. Component plots are stored in the `epiexposure_plot_components`
#'   attribute.
#'
#' @importFrom rlang .data
#' @export
plot_eci <- function(
    eci_overal = NULL,
    eci_lag = NULL,
    facet_scales = c("free_y", "fixed", "free", "free_x"),
    overall_nrow = 3L,
    lag_ncol = 1L,
    sample_smooth = TRUE,
    smooth_method = NULL,
    smooth_formula = NULL,
    linewidth = 0.7,
    alpha = 0.20,
    overall_x_label = "ECI",
    overall_y_label = "Predicted (%)",
    lag_x_label = "lag",
    lag_y_label = "Contribution (%)",
    exposure_y_label = "Exposure difference",
    exposure_colors = c("#4C72B0", "#C44E52", "#55A868"),
    base_size = 10,
    panel_labels = c("(a)", "(b)"),
    label_size = 12,
    rel_widths = NULL
) {
  stopf <- function(...) stop(..., call. = FALSE)
  scalar_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  scalar_number <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  }
  positive_integer <- function(x) {
    scalar_number(x) && x >= 1 && x == as.integer(x)
  }
  require_columns <- function(data, columns, object_name) {
    missing <- setdiff(columns, names(data))
    if (length(missing)) {
      stopf("`", object_name, "` is missing required column(s): ",
            paste(missing, collapse = ", "), ".")
    }
  }
  first_existing <- function(data, candidates) {
    hit <- candidates[candidates %in% names(data)]
    if (length(hit)) hit[[1L]] else NULL
  }
  validate_frame <- function(x, name) {
    if (!is.data.frame(x) || !nrow(x)) {
      stopf("`", name, "` must be a non-empty data.frame.")
    }
    x
  }
  add_sample_scale <- function(plot, data) {
    if (is.numeric(data$sample)) {
      plot + ggplot2::scale_color_viridis_c()
    } else {
      plot + ggplot2::scale_color_viridis_d()
    }
  }
  add_sample_curves <- function(plot, mapping, data) {
    if (isTRUE(sample_smooth)) {
      args <- list(
        mapping = mapping,
        data = data,
        se = FALSE,
        linetype = 1,
        linewidth = linewidth
      )
      if (!is.null(smooth_method)) args$method <- smooth_method
      if (!is.null(smooth_formula)) args$formula <- smooth_formula
      plot + do.call(ggplot2::geom_smooth, args)
    } else {
      plot + ggplot2::geom_line(
        mapping = mapping,
        data = data,
        linewidth = linewidth
      )
    }
  }

  if (is.null(eci_overal) && is.null(eci_lag)) {
    stopf("Supply at least one of `eci_overal` or `eci_lag`.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stopf("Package 'ggplot2' is required by `plot_eci()`.")
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stopf("Package 'cowplot' is required by `plot_eci()`.")
  }

  facet_scales <- match.arg(facet_scales)
  if (!positive_integer(overall_nrow)) stopf("`overall_nrow` must be a positive integer.")
  if (!positive_integer(lag_ncol)) stopf("`lag_ncol` must be a positive integer.")
  if (!scalar_flag(sample_smooth)) stopf("`sample_smooth` must be TRUE or FALSE.")
  if (!scalar_number(linewidth) || linewidth <= 0) stopf("`linewidth` must be positive.")
  if (!scalar_number(alpha) || alpha < 0 || alpha > 1) stopf("`alpha` must lie in [0, 1].")
  if (!scalar_number(base_size) || base_size <= 0) stopf("`base_size` must be positive.")
  if (!is.character(panel_labels) ||
      anyNA(panel_labels) ||
      length(panel_labels) < 1L) {
    stopf(
      "`panel_labels` must contain at least one character label."
    )
  }
  if (!scalar_number(label_size) || label_size <= 0) stopf("`label_size` must be positive.")
  if (!is.null(rel_widths) &&
      (!is.numeric(rel_widths) || !length(rel_widths) ||
       anyNA(rel_widths) || any(!is.finite(rel_widths)) ||
       any(rel_widths <= 0))) {
    stopf("`rel_widths` must be NULL or contain positive finite numbers.")
  }

  common_theme <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size, face = "bold"),
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "white", colour = "black")
    )

  overall_plot <- NULL
  lag_contribution_plot <- NULL
  lag_exposure_plot <- NULL
  lag_panel <- NULL

  if (!is.null(eci_overal)) {
    overall_data <- validate_frame(eci_overal, "eci_overal")
    require_columns(overall_data, c("var", "ECI_weighted", "predicted"), "eci_overal")
    overall_is_samples <- "sample" %in% names(overall_data)

    overall_plot <- ggplot2::ggplot(overall_data)

    if (overall_is_samples) {
      overall_plot <- add_sample_curves(
        overall_plot,
        ggplot2::aes(
          x = .data$ECI_weighted,
          y = .data$predicted,
          group = .data$sample,
          color = .data$sample
        ),
        overall_data
      )
      overall_plot <- add_sample_scale(overall_plot, overall_data)
    } else {
      eci_lower <- first_existing(overall_data, c("ECI_weighted_lower", "ECI_lower"))
      eci_upper <- first_existing(overall_data, c("ECI_weighted_upper", "ECI_upper"))
      pred_lower <- first_existing(overall_data, c("predicted_lower", "prediction_lower"))
      pred_upper <- first_existing(overall_data, c("predicted_upper", "prediction_upper"))

      if (!is.null(pred_lower) && !is.null(pred_upper)) {
        overall_plot <- overall_plot + ggplot2::geom_ribbon(
          ggplot2::aes(
            x = .data$ECI_weighted,
            ymin = .data[[pred_lower]],
            ymax = .data[[pred_upper]],
            group = .data$var
          ),
          alpha = alpha,
          fill = "grey60",
          colour = NA
        )
      }

      overall_plot <- overall_plot + ggplot2::geom_line(
        ggplot2::aes(x = .data$ECI_weighted, y = .data$predicted, group = .data$var),
        linewidth = linewidth,
        colour = "black"
      )

      if (!is.null(eci_lower) && !is.null(eci_upper)) {
        overall_plot <- overall_plot + ggplot2::geom_errorbarh(
          ggplot2::aes(
            y = .data$predicted,
            xmin = .data[[eci_lower]],
            xmax = .data[[eci_upper]]
          ),
          height = 0,
          alpha = max(alpha, 0.35),
          linewidth = linewidth * 0.6,
          colour = "grey35"
        )
      }
    }

    overall_plot <- overall_plot +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$var),
        nrow = as.integer(overall_nrow),
        scales = facet_scales
      ) +
      ggplot2::labs(x = overall_x_label, y = overall_y_label) +
      common_theme
  }

  if (!is.null(eci_lag)) {
    if (is.data.frame(eci_lag)) {
      lag_base <- validate_frame(eci_lag, "eci_lag")
      lag_samples <- if ("sample" %in% names(lag_base)) lag_base else NULL
      lag_summary <- if (!"sample" %in% names(lag_base)) lag_base else NULL
    } else if (is.list(eci_lag)) {
      lag_samples <- eci_lag$by_lag_samples
      lag_summary <- eci_lag$by_lag
      if (!is.null(lag_samples)) lag_samples <- validate_frame(lag_samples, "eci_lag$by_lag_samples")
      if (!is.null(lag_summary)) lag_summary <- validate_frame(lag_summary, "eci_lag$by_lag")
    } else {
      stopf("`eci_lag` must be a compute_ecilag result or a non-empty data.frame.")
    }

    contribution_data <- if (!is.null(lag_samples)) lag_samples else lag_summary
    if (is.null(contribution_data)) {
      stopf("`eci_lag` does not contain `by_lag_samples` or `by_lag` data.")
    }
    require_columns(contribution_data, c("var", "lag"), "lag contribution data")

    contribution_is_samples <- "sample" %in% names(contribution_data)
    contribution_col <- if (contribution_is_samples) {
      first_existing(contribution_data, c("percent_contribution", "contribution"))
    } else {
      first_existing(contribution_data, c("contribution", "percent_contribution"))
    }
    if (is.null(contribution_col)) {
      stopf("Lag contribution data must contain `percent_contribution` or `contribution`.")
    }

    lag_contribution_plot <- ggplot2::ggplot(contribution_data)

    if (contribution_is_samples) {
      lag_contribution_plot <- add_sample_curves(
        lag_contribution_plot,
        ggplot2::aes(
          x = .data$lag,
          y = .data[[contribution_col]],
          group = .data$sample,
          color = .data$sample
        ),
        contribution_data
      )
      lag_contribution_plot <- add_sample_scale(lag_contribution_plot, contribution_data)
    } else {
      contribution_lower <- first_existing(
        contribution_data,
        c("contribution_lower", "contributon_lower", "percent_contribution_lower")
      )
      contribution_upper <- first_existing(
        contribution_data,
        c("contribution_upper", "percent_contribution_upper")
      )

      if (!is.null(contribution_lower) && !is.null(contribution_upper)) {
        lag_contribution_plot <- lag_contribution_plot + ggplot2::geom_ribbon(
          ggplot2::aes(
            x = .data$lag,
            ymin = .data[[contribution_lower]],
            ymax = .data[[contribution_upper]],
            group = .data$var
          ),
          alpha = alpha,
          fill = "grey60",
          colour = NA
        )
      }

      lag_contribution_plot <- lag_contribution_plot + ggplot2::geom_line(
        ggplot2::aes(x = .data$lag, y = .data[[contribution_col]], group = .data$var),
        linewidth = linewidth,
        colour = "black"
      )
    }

    lag_contribution_plot <- lag_contribution_plot +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$var),
        ncol = as.integer(lag_ncol),
        scales = facet_scales
      ) +
      ggplot2::labs(x = lag_x_label, y = lag_y_label) +
      common_theme

    exposure_data <- lag_summary
    if (is.null(exposure_data) && !is.null(lag_samples) &&
        "exposure_minus_reference" %in% names(lag_samples)) {
      grouping <- intersect(c("var", "lag", "epi_id", "group_level"), names(lag_samples))
      if (!length(grouping)) grouping <- c("var", "lag")
      exposure_data <- stats::aggregate(
        lag_samples$exposure_minus_reference,
        by = lag_samples[grouping],
        FUN = stats::median,
        na.rm = TRUE
      )
      names(exposure_data)[ncol(exposure_data)] <- "exposure_minus_reference"
    }

    if (!is.null(exposure_data) &&
        all(c("var", "lag", "exposure_minus_reference") %in% names(exposure_data))) {
      lag_exposure_plot <- ggplot2::ggplot(
        exposure_data,
        ggplot2::aes(
          x = .data$lag,
          y = .data$exposure_minus_reference,
          color = .data$var,
          group = .data$var
        )
      ) +
        ggplot2::geom_line(linewidth = linewidth) +
        ggplot2::facet_wrap(
          ggplot2::vars(.data$var),
          ncol = as.integer(lag_ncol),
          scales = facet_scales
        ) +
        ggplot2::scale_color_manual(values = exposure_colors) +
        ggplot2::labs(x = lag_x_label, y = exposure_y_label) +
        common_theme
    }

    lag_panel <- if (is.null(lag_exposure_plot)) {
      lag_contribution_plot
    } else {
      cowplot::plot_grid(
        lag_contribution_plot,
        lag_exposure_plot,
        nrow = 1,
        rel_widths = if (is.null(rel_widths)) c(1, 1) else rel_widths[seq_len(2L)]
      )
    }
  }

  components <- list(
    overall = overall_plot,
    lag_contribution = lag_contribution_plot,
    lag_exposure = lag_exposure_plot,
    lag_panel = lag_panel
  )

  # Combine the original ggplot objects directly instead of nesting the two
  # lag plots inside a second grid. This gives every displayed panel the same
  # default width and reproduces `plot_grid(p1, p2, p3, nrow = 1)`.
  if (!is.null(overall_plot) && !is.null(lag_contribution_plot) &&
      !is.null(lag_exposure_plot)) {
    widths <- if (is.null(rel_widths)) c(1, 1, 1) else rel_widths
    if (length(widths) != 3L) {
      stopf("`rel_widths` must contain three values when all three panels are displayed.")
    }
    labels_use <- if (length(panel_labels) >= 3L) {
      panel_labels[1:3]
    } else {
      c(panel_labels, rep("", 3L - length(panel_labels)))
    }

    out <- cowplot::plot_grid(
      overall_plot,
      lag_contribution_plot,
      lag_exposure_plot,
      nrow = 1,
      labels = labels_use,
      label_size = label_size,
      rel_widths = widths
    )
  } else if (!is.null(overall_plot) && !is.null(lag_contribution_plot)) {
    widths <- if (is.null(rel_widths)) c(1, 1) else rel_widths
    if (length(widths) != 2L) {
      stopf("`rel_widths` must contain two values when two panels are displayed.")
    }
    out <- cowplot::plot_grid(
      overall_plot,
      lag_contribution_plot,
      nrow = 1,
      labels = panel_labels[1:2],
      label_size = label_size,
      rel_widths = widths
    )
  } else if (!is.null(overall_plot)) {
    out <- overall_plot
  } else if (!is.null(lag_contribution_plot) && !is.null(lag_exposure_plot)) {
    widths <- if (is.null(rel_widths)) c(1, 1) else rel_widths
    if (length(widths) != 2L) {
      stopf("`rel_widths` must contain two values when two lag panels are displayed.")
    }
    out <- cowplot::plot_grid(
      lag_contribution_plot,
      lag_exposure_plot,
      nrow = 1,
      rel_widths = widths
    )
  } else {
    out <- lag_contribution_plot
  }

  attr(out, "epiexposure_plot_components") <- components
  attr(out, "epiexposure_plot_type") <- if (!is.null(overall_plot) && !is.null(lag_panel)) {
    "overall_and_lag"
  } else if (!is.null(overall_plot)) {
    "overall"
  } else {
    "lag"
  }
  out
}
