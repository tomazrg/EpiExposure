# ============================================================================
# plot_effects()
# Default design reproduces the manually created EpiExposure figures:
#   * lag-specific: independent plot for each variable and metric
#   * effect: blue-white-red diverging palette, symmetric clipping at 2-98%
#   * delta: viridis palette, symmetric clipping at 2-98%
#   * effect row above delta row; x-axis labels hidden on the upper row
#   * independent legends for every variable/metric panel
#   * period-specific: smoothed sample curves, free facets, viridis colors,
#     white facet strips, no legend
# ============================================================================

.resolve_plot_theme <- function(theme, base_size) {
  if (inherits(theme, "theme")) {
    return(theme)
  }

  theme <- match.arg(theme, c("bw", "minimal", "classic"))

  switch(
    theme,
    bw = ggplot2::theme_bw(base_size = base_size),
    minimal = ggplot2::theme_minimal(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size)
  )
}

.resolve_named_labels <- function(x, keys, defaults, arg = "labels") {
  out <- defaults[keys]

  missing_default <- is.na(out)
  out[missing_default] <- keys[missing_default]
  names(out) <- keys

  if (is.null(x)) {
    return(out)
  }

  if (!is.character(x)) {
    stop("`", arg, "` must be a character vector.", call. = FALSE)
  }

  if (!is.null(names(x)) && any(nzchar(names(x)))) {
    supplied_names <- names(x)[nzchar(names(x))]
    idx <- intersect(supplied_names, keys)

    # Extra named labels are intentionally ignored so a complete named vector
    # can still be supplied when only a subset of variables/metrics is plotted.
    if (length(idx) == 0L) {
      stop(
        "None of the names supplied in `", arg,
        "` match the variables/metrics being plotted.",
        call. = FALSE
      )
    }

    out[idx] <- x[idx]
    return(out)
  }

  if (length(x) == length(keys)) {
    out[] <- x
    return(out)
  }

  if (length(x) == 1L && length(keys) == 1L) {
    out[] <- x
    return(out)
  }

  stop(
    "When more than one variable/metric is plotted, `", arg,
    "` must be a named vector or an unnamed vector with the same length.",
    call. = FALSE
  )
}

.safe_symmetric_limit <- function(x, probs = c(0.02, 0.98)) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(1)
  }

  lims <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  max_abs <- max(abs(lims), na.rm = TRUE)

  if (!is.finite(max_abs) || max_abs <= 0) {
    max_abs <- max(abs(x), na.rm = TRUE)
  }

  if (!is.finite(max_abs) || max_abs <= 0) {
    max_abs <- 1
  }

  max_abs
}

.default_period_order <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]

  if (length(x) == 0L) {
    return(character())
  }

  # Reproduces Period4, Period3, Period2, Period1 when names follow PeriodN.
  if (all(grepl("^Period[0-9]+$", x))) {
    n <- as.integer(sub("^Period", "", x))
    return(x[order(n, decreasing = TRUE)])
  }

  rev(x)
}

#' Plot lag-specific or period-specific effects
#'
#' @param data Data frame returned by summarise_effects().
#' @param scale Plot scale. Use only \"lag\" or \"period\".
#' @param vars Character vector with variables to plot. NULL uses all variables.
#' @param metric For lag plots, one or both of c("effect", "delta").
#' @param delta_multiplier Multiplier applied to delta before plotting.
#' @param metric_labels Named labels for the lag color legends.
#' @param ylab Lag y-axis labels. Prefer a named vector keyed by variable.
#' @param metric_ncol Number of metric blocks per row. Default 1 places effect
#'   above delta, matching the manual figure.
#' @param vars_ncol Number of variable panels per metric block. NULL uses all
#'   selected variables in one row.
#' @param effect_palette Three colors: low, mid, high.
#' @param delta_palette NULL uses viridis; otherwise a vector of >= 2 colors.
#' @param delta_viridis_option Viridis option used when delta_palette is NULL.
#' @param clip_quantiles Quantiles used to clip lag fill values.
#' @param effect_breaks Contour breaks for effect.
#' @param delta_bins Number of contour bins for delta.
#' @param lag_reverse Reverse the lag axis.
#' @param lag_xlab X-axis label for lag plots.
#' @param legend_position Named positions for effect and delta legends.
#' @param lag_theme One of "bw", "minimal", "classic", or a ggplot2 theme.
#' @param lag_base_size Base size for lag theme.
#' @param axis_title_size Axis-title size for lag plots.
#' @param axis_text_size Axis-text size for lag plots.
#' @param contour_colour Contour color.
#' @param contour_linewidth Contour line width.
#' @param zero_linewidth Delta zero-contour line width.
#' @param effect_zero_text_size Effect zero-contour label size.
#' @param delta_zero_text_size Delta zero-contour label size.
#' @param interpolate Passed to geom_raster().
#' @param period_order Period facet order. NULL reproduces reversed PeriodN order.
#' @param vars_order Period variable order. NULL follows `vars`.
#' @param period_response Response column for period plots; default "eta".
#' @param period_exclude_samples Samples removed from period plot. Default 10
#'   reproduces the supplied manual plot; use NULL to keep every sample.
#' @param period_palette NULL uses viridis; otherwise vector of >= 2 colors.
#' @param period_viridis_option Viridis option for period curves.
#' @param period_xlab X-axis label for period plots.
#' @param period_ylab Y-axis label for period plots.
#' @param period_theme One of "bw", "minimal", "classic", or a ggplot2 theme.
#' @param period_base_size Base size for period theme. Default 11 matches theme_bw().
#' @param period_linewidth Line width passed to geom_smooth().
#' @param period_show_legend Whether to show the period color legend.
#' @param panel_labels Character vector used to label the lag-specific metric
#'   blocks in the combined plot. Default is `c("(a)", "(b)")`, so the
#'   first selected metric block is labelled "(a)" and the second "(b)".
#'   Labels are applied only when `scale = "lag"` and `output = "plot"`;
#'   `output = "list"` returns the original unlabelled component plots.
#' @param label_size Positive finite size used for `panel_labels`. Default is `12`.
#' @param output Either "plot" or "list". For lag, "list" returns every panel.
#'
#' @return A ggplot/cowplot object, or a list of ggplots when output = "list".
#' @export
plot_effects <- function(
    data,
    scale = c("lag", "period"),
    vars = NULL,
    metric = c("effect", "delta"),
    delta_multiplier = 100,
    metric_labels = c(effect = "Effect", delta = "Delta"),
    ylab = NULL,
    metric_ncol = 1,
    vars_ncol = NULL,
    effect_palette = c(low = "#3B4CC0", mid = "white", high = "#B40426"),
    delta_palette = NULL,
    delta_viridis_option = "D",
    clip_quantiles = c(0.02, 0.98),
    effect_breaks = c(-15, -10, -5, 0, 5, 10, 15),
    delta_bins = 10,
    lag_reverse = TRUE,
    lag_xlab = "Lag",
    legend_position = c(effect = "top", delta = "bottom"),
    lag_theme = "bw",
    lag_base_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    contour_colour = "black",
    contour_linewidth = 0.25,
    zero_linewidth = 0.8,
    effect_zero_text_size = 3.2,
    delta_zero_text_size = 3,
    interpolate = TRUE,
    period_order = NULL,
    vars_order = NULL,
    period_response = "eta",
    period_exclude_samples = 10,
    period_palette = NULL,
    period_viridis_option = "D",
    period_xlab = "Exposure",
    period_ylab = "Response",
    period_theme = "bw",
    period_base_size = 11,
    period_linewidth = 1,
    period_show_legend = FALSE,
    panel_labels = c("(a)", "(b)"),
    label_size = 12,
    output = c("plot", "list")
) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  scale <- match.arg(scale)
  output <- match.arg(output)

  if (!is.character(panel_labels) ||
      !length(panel_labels) ||
      anyNA(panel_labels)) {
    stop(
      "`panel_labels` must contain at least one non-missing character label.",
      call. = FALSE
    )
  }

  if (!is.numeric(label_size) ||
      length(label_size) != 1L ||
      is.na(label_size) ||
      !is.finite(label_size) ||
      label_size <= 0) {
    stop("`label_size` must be one positive finite number.", call. = FALSE)
  }

  if (!"var" %in% names(data)) {
    stop("`data` must contain a `var` column.", call. = FALSE)
  }

  available_vars <- unique(as.character(data$var))
  available_vars <- available_vars[!is.na(available_vars)]

  if (is.null(vars)) {
    preferred <- c("tmean", "rain", "wetness")
    vars <- c(
      intersect(preferred, available_vars),
      setdiff(available_vars, preferred)
    )
  } else {
    vars <- unique(as.character(vars))
    missing_vars <- setdiff(vars, available_vars)
    if (length(missing_vars) > 0L) {
      stop(
        "Variable(s) not found in `data`: ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }
  }

  dat <- data[as.character(data$var) %in% vars, , drop = FALSE]

  # ==========================================================================
  # LAG-SPECIFIC EFFECTS
  # ==========================================================================
  if (scale == "lag") {

    required <- c("lag", "value")
    missing_cols <- setdiff(required, names(dat))
    if (length(missing_cols) > 0L) {
      stop(
        "Lag plots require column(s): ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    metric <- unique(as.character(metric))
    bad_metric <- setdiff(metric, c("effect", "delta"))
    if (length(bad_metric) > 0L) {
      stop(
        "Unknown metric(s): ", paste(bad_metric, collapse = ", "),
        ". Use `effect`, `delta`, or both.",
        call. = FALSE
      )
    }

    # Fixed default order: effect first, delta second.
    metric <- intersect(c("effect", "delta"), metric)
    if (length(metric) == 0L) {
      stop("Select at least one of `effect` or `delta`.", call. = FALSE)
    }

    missing_metrics <- setdiff(metric, names(dat))
    if (length(missing_metrics) > 0L) {
      stop(
        "Metric column(s) not found in `data`: ",
        paste(missing_metrics, collapse = ", "),
        call. = FALSE
      )
    }

    if (length(clip_quantiles) != 2L ||
        any(!is.finite(clip_quantiles)) ||
        clip_quantiles[1] < 0 ||
        clip_quantiles[2] > 1 ||
        clip_quantiles[1] >= clip_quantiles[2]) {
      stop(
        "`clip_quantiles` must contain two increasing probabilities between 0 and 1.",
        call. = FALSE
      )
    }

    if (length(effect_palette) != 3L) {
      stop("`effect_palette` must contain exactly three colors.", call. = FALSE)
    }

    if (!is.null(delta_palette) && length(delta_palette) < 2L) {
      stop("`delta_palette` must contain at least two colors.", call. = FALSE)
    }

    if (!is.numeric(metric_ncol) || length(metric_ncol) != 1L ||
        !is.finite(metric_ncol) || metric_ncol < 1) {
      stop("`metric_ncol` must be a positive integer.", call. = FALSE)
    }
    metric_ncol <- as.integer(metric_ncol)

    if (is.null(vars_ncol)) {
      vars_ncol <- length(vars)
    }
    if (!is.numeric(vars_ncol) || length(vars_ncol) != 1L ||
        !is.finite(vars_ncol) || vars_ncol < 1) {
      stop("`vars_ncol` must be a positive integer.", call. = FALSE)
    }
    vars_ncol <- as.integer(vars_ncol)

    default_ylab <- c(
      tmean = "Mean temperature (\u00B0C)",
      rain = "Daily precipitation (mm)",
      wetness = "Leaf wetness (%)"
    )

    ylab_use <- .resolve_named_labels(
      x = ylab,
      keys = vars,
      defaults = default_ylab,
      arg = "ylab"
    )

    default_metric_labels <- c(effect = "Effect", delta = "Delta")
    metric_labels_use <- .resolve_named_labels(
      x = metric_labels,
      keys = metric,
      defaults = default_metric_labels,
      arg = "metric_labels"
    )

    default_legend_position <- c(effect = "top", delta = "bottom")
    legend_position_use <- .resolve_named_labels(
      x = legend_position,
      keys = metric,
      defaults = default_legend_position,
      arg = "legend_position"
    )

    theme_lag <- .resolve_plot_theme(lag_theme, lag_base_size)

    make_lag_panel <- function(v, m, hide_x = FALSE) {

      df <- dat[as.character(dat$var) == v, , drop = FALSE]

      z <- df[[m]]
      if (m == "delta") {
        z <- z * delta_multiplier
      }

      max_abs <- .safe_symmetric_limit(z, probs = clip_quantiles)

      df$z_raw <- z
      df$z_fill <- pmax(pmin(z, max_abs), -max_abs)

      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = .data[["lag"]],
          y = .data[["value"]]
        )
      ) +
        ggplot2::geom_raster(
          ggplot2::aes(
            fill = .data[["z_fill"]]
          ),
          interpolate = interpolate
        )

      if (m == "effect") {
        p <- p +
          ggplot2::geom_contour(
            ggplot2::aes(z = .data[["z_raw"]]),
            color = contour_colour,
            linewidth = contour_linewidth,
            breaks = effect_breaks
          ) +
          metR::geom_text_contour(
            ggplot2::aes(z = .data[["z_raw"]]),
            breaks = 0,
            stroke = 0.15,
            size = effect_zero_text_size
          ) +
          ggplot2::scale_fill_gradient2(
            low = unname(effect_palette[1]),
            mid = unname(effect_palette[2]),
            high = unname(effect_palette[3]),
            midpoint = 0,
            limits = c(-max_abs, max_abs),
            oob = scales::squish,
            name = unname(metric_labels_use[[m]]),
            guide = ggplot2::guide_colorbar(
              barwidth = grid::unit(4, "cm"),
              barheight = grid::unit(0.5, "cm")
            )
          )
      } else {
        p <- p +
          ggplot2::geom_contour(
            ggplot2::aes(z = .data[["z_raw"]]),
            color = contour_colour,
            linewidth = contour_linewidth,
            bins = delta_bins
          ) +
          ggplot2::geom_contour(
            ggplot2::aes(z = .data[["z_raw"]]),
            breaks = 0,
            color = contour_colour,
            linewidth = zero_linewidth
          ) +
          metR::geom_text_contour(
            ggplot2::aes(z = .data[["z_raw"]]),
            breaks = 0,
            stroke = 0.15,
            size = delta_zero_text_size
          )

        if (is.null(delta_palette)) {
          p <- p + ggplot2::scale_fill_viridis_c(
            option = delta_viridis_option,
            name = unname(metric_labels_use[[m]])
          )
        } else {
          p <- p + ggplot2::scale_fill_gradientn(
            colours = delta_palette,
            name = unname(metric_labels_use[[m]])
          )
        }
      }

      if (isTRUE(lag_reverse)) {
        p <- p + ggplot2::scale_x_reverse(
          breaks = pretty(df$lag, n = 8),
          expand = c(0, 0)
        )
      } else {
        p <- p + ggplot2::scale_x_continuous(
          breaks = pretty(df$lag, n = 8),
          expand = c(0, 0)
        )
      }

      p <- p +
        ggplot2::scale_y_continuous(
          breaks = pretty(df$value, n = 6),
          expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(
          x = if (hide_x) "" else lag_xlab,
          y = unname(ylab_use[[v]])
        ) +
        theme_lag +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(
            face = "bold",
            size = axis_title_size
          ),
          axis.text = ggplot2::element_text(size = axis_text_size),
          legend.title = ggplot2::element_text(face = "bold"),
          legend.position = unname(legend_position_use[[m]])
        )

      if (hide_x) {
        p <- p + ggplot2::theme(
          axis.text.x = ggplot2::element_blank()
        )
      }

      p
    }

    panels <- list()
    panel_names <- character()

    # metric-major ordering reproduces:
    # tmean effect, rain effect, wetness effect,
    # tmean delta,  rain delta,  wetness delta.
    for (m in metric) {
      for (v in vars) {
        hide_x <- length(metric) > 1L &&
          metric_ncol == 1L &&
          m != metric[length(metric)]

        panels[[length(panels) + 1L]] <- make_lag_panel(v, m, hide_x = hide_x)
        panel_names <- c(panel_names, paste(v, m, sep = "__"))
      }
    }

    names(panels) <- panel_names

    if (output == "list") {
      return(panels)
    }

    # With metric_ncol = 1 and vars_ncol = length(vars), this is exactly the
    # manual 2-row layout. metric_ncol = 2 can place two metric blocks side by side.
    final_ncol <- min(length(panels), vars_ncol * metric_ncol)

    # Label metric blocks without changing panel order, panel content, scales,
    # smoothing, or layout. Because panels are stored in metric-major order, the
    # first panel of each metric block receives the corresponding label.
    labels_use <- if (length(panel_labels) >= length(metric)) {
      panel_labels[seq_along(metric)]
    } else {
      c(
        panel_labels,
        rep("", length(metric) - length(panel_labels))
      )
    }

    panel_labels_use <- rep("", length(panels))
    metric_start <- seq.int(
      from = 1L,
      by = length(vars),
      length.out = length(metric)
    )
    panel_labels_use[metric_start] <- labels_use

    return(
      cowplot::plot_grid(
        plotlist = panels,
        ncol = final_ncol,
        labels = panel_labels_use,
        label_size = label_size,
        label_x = -0.02
      )
    )
  }

  # ==========================================================================
  # PERIOD-SPECIFIC EFFECTS
  # ==========================================================================

  required <- c("period", "value", period_response)
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0L) {
    stop(
      "Period plots require column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(period_order)) {
    period_order <- .default_period_order(dat$period)
  } else {
    period_order <- as.character(period_order)
  }

  available_periods <- unique(as.character(dat$period))
  missing_periods <- setdiff(period_order, available_periods)
  if (length(missing_periods) > 0L) {
    stop(
      "Period(s) not found in `data`: ",
      paste(missing_periods, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(vars_order)) {
    vars_order <- vars
  } else {
    vars_order <- as.character(vars_order)
  }

  missing_period_vars_sel <- setdiff(vars_order, vars)
  if (length(missing_period_vars_sel) > 0L) {
    stop(
      "Variable(s) in `vars_order` were not selected in `vars`: ",
      paste(missing_period_vars_sel, collapse = ", "),
      call. = FALSE
    )
  }

  df <- dat

  if ("sample" %in% names(df) && !is.null(period_exclude_samples)) {
    df <- df[!(df$sample %in% period_exclude_samples), , drop = FALSE]
  }

  df$period <- factor(as.character(df$period), levels = period_order)
  df$var <- factor(as.character(df$var), levels = vars_order)
  df$response_value <- df[[period_response]]

  theme_period <- .resolve_plot_theme(period_theme, period_base_size)

  if ("sample" %in% names(df)) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[["value"]],
        y = .data[["response_value"]],
        group = .data[["sample"]],
        color = .data[["sample"]]
      )
    ) +
      ggplot2::geom_smooth(
        se = FALSE,
        linewidth = period_linewidth
      )

    if (is.null(period_palette)) {
      p <- p + ggplot2::scale_color_viridis_c(
        option = period_viridis_option
      )
    } else {
      if (length(period_palette) < 2L) {
        stop("`period_palette` must contain at least two colors.", call. = FALSE)
      }
      p <- p + ggplot2::scale_color_gradientn(
        colours = period_palette
      )
    }
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[["value"]],
        y = .data[["response_value"]],
        group = 1
      )
    ) +
      ggplot2::geom_smooth(
        se = FALSE,
        linewidth = period_linewidth
      )
  }

  p <- p +
    ggplot2::facet_wrap(
      stats::as.formula("var ~ period"),
      scales = "free"
    ) +
    theme_period +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black"
      ),
      legend.position = if (isTRUE(period_show_legend)) "right" else "none"
    ) +
    ggplot2::labs(
      x = period_xlab,
      y = period_ylab
    )

  if (output == "list") {
    return(list(period = p))
  }

  p
}
