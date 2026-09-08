#' Plot predictions across epidemiological scenarios
#'
#' Creates scenario-response plots from the output of `simulate_scenarios()`.
#' The function supports sample-level predictions returned by
#' `output = "samples"` and summarized predictions returned by
#' `output = "summary"`.
#'
#' For sample-level output, one smoothed line is drawn for each coefficient or
#' posterior sample and colors are mapped to `sample`. The default color scale
#' is `ggplot2::scale_color_viridis_c()`. For summary output, the central
#' prediction is displayed as a line and, when available, `lower` and `upper`
#' are displayed as an uncertainty ribbon.
#'
#' @param data A non-empty data frame returned by `simulate_scenarios()`.
#' @param x Character scalar naming the scenario or exposure variable displayed
#'   on the horizontal axis, such as `"rain"`.
#' @param y Character scalar naming the prediction column. Default is
#'   `"prediction"`.
#' @param facet Optional character vector naming variables used to define facet
#'   panels, such as `c("wetness", "tmean")`. If `NULL`, no faceting is used.
#' @param sample Character scalar naming the sample identifier used when
#'   `data` contains draw-specific predictions. Default is `"sample"`.
#' @param lower,upper Character scalars naming the lower and upper uncertainty
#'   interval columns used for summary output. Defaults are `"lower"` and
#'   `"upper"`.
#' @param mode Character. `"auto"`, `"samples"`, or `"summary"`. With
#'   `"auto"`, sample mode is selected when the `sample` column is present;
#'   otherwise summary mode is used.
#' @param scale_factor Finite numeric scalar multiplying `y`, `lower`, and
#'   `upper`. Use `100` to display proportional predictions as percentages.
#' @param facet_ncol,facet_nrow Optional positive integers controlling the facet
#'   arrangement.
#' @param color_scale Optional continuous ggplot2 color scale used in sample
#'   mode. When `NULL`, `ggplot2::scale_color_viridis_c()` is used. Any
#'   compatible continuous color scale may be supplied, including a customized
#'   `scale_color_viridis_c()`, `scale_color_gradient()`,
#'   `scale_color_gradientn()`, or `scale_color_distiller()`.
#' @param color_label Legend title used by the default color scale. Default is
#'   `"Sample"`.
#' @param line_color Color used for the central prediction line in summary mode.
#' @param ribbon_fill Fill color used for the uncertainty ribbon in summary
#'   mode.
#' @param ribbon_alpha Alpha transparency used for the uncertainty ribbon.
#' @param show_ribbon Logical. Display an uncertainty ribbon in summary mode
#'   when the `lower` and `upper` columns exist.
#' @param smooth Logical. If `TRUE`, sample trajectories and summary predictions
#'   are drawn with `ggplot2::geom_smooth()`. If `FALSE`, they are drawn with
#'   `ggplot2::geom_line()`.
#' @param smooth_method Smoothing method passed to `ggplot2::geom_smooth()`.
#'   Default is `"loess"`, matching ggplot2 behavior for smaller datasets.
#' @param smooth_formula Optional smoothing formula. If `NULL`, ggplot2 uses the
#'   default formula for the selected method.
#' @param smooth_se Logical. Display the uncertainty band calculated by
#'   `geom_smooth()`. The default is `FALSE` because scenario uncertainty is
#'   represented by draws or by the supplied interval columns.
#' @param line_width Numeric line width used for trajectories.
#' @param line_alpha Numeric alpha transparency used for trajectories.
#' @param x_breaks,y_breaks Optional numeric axis breaks.
#' @param x_limits,y_limits Optional numeric vectors of length two defining axis
#'   limits. Scale limits remove observations outside the specified range.
#' @param x_expand,y_expand Expansion passed to the continuous axes. The default
#'   x expansion reproduces the narrow margins used in the original plot.
#' @param x_label,y_label Axis labels. If `x_label = NULL`, a title-case label is
#'   generated from `x`. Default `y_label` is `"Predicted"`.
#' @param title,subtitle,caption Optional plot annotations.
#' @param theme A ggplot2 theme. When `NULL`,
#'   `ggplot2::theme_bw(base_size = base_size)` is used.
#' @param base_size Base font size used by the default theme.
#' @param text_face Font face used for common plot text.
#' @param axis_title_size Axis-title size.
#' @param strip_background Facet-strip background element.
#' @param legend_position Legend position. Default is `"none"`, matching the
#'   original sample-trajectory plot.
#' @param line_args Optional named list of additional arguments passed to the
#'   line or smooth layer.
#' @param ribbon_args Optional named list of additional arguments passed to
#'   `ggplot2::geom_ribbon()` in summary mode.
#' @param add Optional ggplot2 component or list of components added after the
#'   default plot is assembled.
#' @param ... Additional named arguments passed to the trajectory layer. Entries
#'   in `line_args` take precedence over entries supplied through `...`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' `plot_scenarios()` does not calculate predictions. Predictions must first be
#' produced by `simulate_scenarios()`.
#'
#' In sample mode, `sample` is converted to numeric for the default continuous
#' viridis scale. Every sample is retained as an independent trajectory through
#' `group = sample`. If a custom discrete color scale is desired, set
#' `color_scale` to that scale and use `add` to replace or extend the default
#' mapping as needed.
#'
#' In summary mode, the function looks for the columns specified by `lower` and
#' `upper`. When both are present and `show_ribbon = TRUE`, the interval is
#' plotted as a ribbon. If either column is absent, only the central prediction
#' is plotted.
#'
#' @examples
#' \dontrun{
#' plot_scenarios(
#'   data = pred_scenarios,
#'   x = "rain",
#'   y = "prediction",
#'   facet = c("wetness", "tmean"),
#'   scale_factor = 100,
#'   x_limits = c(0, 15),
#'   x_label = "Exposure",
#'   y_label = "Predicted",
#'   legend_position = "none"
#' )
#'
#' plot_scenarios(
#'   pred_scenarios,
#'   x = "rain",
#'   facet = c("wetness", "tmean"),
#'   color_scale = ggplot2::scale_color_viridis_c(
#'     option = "magma",
#'     name = "Posterior sample"
#'   ),
#'   facet_ncol = 2,
#'   line_width = 0.7,
#'   line_alpha = 0.8
#' )
#'
#' plot_scenarios(
#'   data = pred_summary,
#'   x = "rain",
#'   mode = "summary",
#'   facet = c("wetness", "tmean"),
#'   scale_factor = 100,
#'   line_color = "navy",
#'   ribbon_fill = "skyblue",
#'   ribbon_alpha = 0.25
#' )
#' }
#'
#' @export
plot_scenarios <- function(
    data,
    x,
    y = "prediction",
    facet = NULL,
    sample = "sample",
    lower = "lower",
    upper = "upper",
    mode = c("auto", "samples", "summary"),
    scale_factor = 1,
    facet_ncol = NULL,
    facet_nrow = NULL,
    color_scale = NULL,
    color_label = "Sample",
    line_color = "black",
    ribbon_fill = "grey70",
    ribbon_alpha = 0.25,
    show_ribbon = TRUE,
    smooth = TRUE,
    smooth_method = "loess",
    smooth_formula = NULL,
    smooth_se = FALSE,
    line_width = 0.8,
    line_alpha = 1,
    x_breaks = NULL,
    y_breaks = NULL,
    x_limits = NULL,
    y_limits = NULL,
    x_expand = ggplot2::expansion(mult = c(0.01, 0.01)),
    y_expand = ggplot2::waiver(),
    x_label = NULL,
    y_label = "Predicted",
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    theme = NULL,
    base_size = 10,
    text_face = "bold",
    axis_title_size = 12,
    strip_background = ggplot2::element_rect(fill = "white", colour = "black"),
    legend_position = "none",
    line_args = list(),
    ribbon_args = list(),
    add = NULL,
    ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required by `plot_scenarios()`.", call. = FALSE)
  }

  scalar_name <- function(z) {
    is.character(z) && length(z) == 1L && !is.na(z) && nzchar(z)
  }
  scalar_flag <- function(z) {
    is.logical(z) && length(z) == 1L && !is.na(z)
  }
  whole_scalar <- function(z) {
    is.numeric(z) && length(z) == 1L && !is.na(z) &&
      is.finite(z) && z == as.integer(z)
  }
  validate_optional_integer <- function(z, argument) {
    if (!is.null(z) && (!whole_scalar(z) || z < 1L)) {
      stop("`", argument, "` must be NULL or a positive integer.", call. = FALSE)
    }
    if (is.null(z)) NULL else as.integer(z)
  }
  validate_numeric <- function(z, argument, length_two = FALSE) {
    if (is.null(z) || inherits(z, "waiver")) return(z)
    if (!is.numeric(z) || !length(z) || anyNA(z) || any(!is.finite(z))) {
      stop("`", argument, "` must contain finite numeric values.", call. = FALSE)
    }
    if (length_two && length(z) != 2L) {
      stop("`", argument, "` must contain exactly two values.", call. = FALSE)
    }
    z
  }
  validate_named_list <- function(z, argument) {
    if (!is.list(z)) stop("`", argument, "` must be a list.", call. = FALSE)
    if (length(z) && (is.null(names(z)) || anyNA(names(z)) ||
                      any(!nzchar(names(z))))) {
      stop("Every entry in `", argument, "` must be named.", call. = FALSE)
    }
    z
  }
  add_components <- function(plot_object, components) {
    if (is.null(components)) return(plot_object)
    if (!is.list(components) ||
        inherits(components, c("gg", "ggplot", "Scale", "theme", "Coord"))) {
      components <- list(components)
    }
    for (component in components) plot_object <- plot_object + component
    plot_object
  }

  if (!is.data.frame(data) || !nrow(data)) {
    stop("`data` must be a non-empty data frame returned by `simulate_scenarios()`.",
         call. = FALSE)
  }
  for (argument in c("x", "y", "sample", "lower", "upper")) {
    if (!scalar_name(get(argument, inherits = FALSE))) {
      stop("`", argument, "` must be one non-empty column name.", call. = FALSE)
    }
  }
  if (!is.null(facet) &&
      (!is.character(facet) || !length(facet) || anyNA(facet) ||
       any(!nzchar(facet)) || anyDuplicated(facet))) {
    stop("`facet` must be NULL or unique non-empty column names.", call. = FALSE)
  }

  mode <- match.arg(mode)
  if (mode == "auto") mode <- if (sample %in% names(data)) "samples" else "summary"

  required <- unique(c(x, y, facet, if (mode == "samples") sample))
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Column(s) missing from `data`: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }

  numeric_required <- unique(c(x, y))
  invalid <- numeric_required[!vapply(data[numeric_required], function(z) {
    is.numeric(z) && !anyNA(z) && all(is.finite(z))
  }, logical(1))]
  if (length(invalid)) {
    stop("`x` and `y` must identify finite numeric columns: ",
         paste(invalid, collapse = ", "), ".", call. = FALSE)
  }

  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      is.na(scale_factor) || !is.finite(scale_factor)) {
    stop("`scale_factor` must be one finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(ribbon_alpha) || length(ribbon_alpha) != 1L ||
      is.na(ribbon_alpha) || !is.finite(ribbon_alpha) ||
      ribbon_alpha < 0 || ribbon_alpha > 1) {
    stop("`ribbon_alpha` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(line_alpha) || length(line_alpha) != 1L ||
      is.na(line_alpha) || !is.finite(line_alpha) ||
      line_alpha < 0 || line_alpha > 1) {
    stop("`line_alpha` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(line_width) || length(line_width) != 1L ||
      is.na(line_width) || !is.finite(line_width) || line_width <= 0) {
    stop("`line_width` must be one positive finite number.", call. = FALSE)
  }
  for (argument in c("show_ribbon", "smooth", "smooth_se")) {
    if (!scalar_flag(get(argument, inherits = FALSE))) {
      stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  facet_ncol <- validate_optional_integer(facet_ncol, "facet_ncol")
  facet_nrow <- validate_optional_integer(facet_nrow, "facet_nrow")
  x_breaks <- validate_numeric(x_breaks, "x_breaks")
  y_breaks <- validate_numeric(y_breaks, "y_breaks")
  x_limits <- validate_numeric(x_limits, "x_limits", TRUE)
  y_limits <- validate_numeric(y_limits, "y_limits", TRUE)
  x_expand <- validate_numeric(x_expand, "x_expand")
  y_expand <- validate_numeric(y_expand, "y_expand")
  line_args <- validate_named_list(line_args, "line_args")
  ribbon_args <- validate_named_list(ribbon_args, "ribbon_args")
  dots <- validate_named_list(list(...), "...")

  if (!is.null(color_scale) && !inherits(color_scale, "Scale")) {
    stop("`color_scale` must be NULL or a ggplot2 color scale.", call. = FALSE)
  }

  plot_data <- data
  plot_data$.plot_x <- plot_data[[x]]
  plot_data$.plot_y <- plot_data[[y]] * scale_factor

  if (mode == "samples") {
    sample_values <- plot_data[[sample]]
    if (is.numeric(sample_values)) {
      plot_data$.plot_sample <- as.numeric(sample_values)
    } else {
      sample_numeric <- suppressWarnings(as.numeric(as.character(sample_values)))
      if (anyNA(sample_numeric)) {
        sample_numeric <- as.numeric(factor(sample_values, levels = unique(sample_values)))
      }
      plot_data$.plot_sample <- sample_numeric
    }
    if (length(facet)) {

      facet_id <- do.call(
        interaction,
        c(
          lapply(
            facet,
            function(v) plot_data[[v]]
          ),
          drop = TRUE
        )
      )

      plot_data$.plot_group <- interaction(
        plot_data[[sample]],
        facet_id,
        drop = TRUE
      )

    } else {

      plot_data$.plot_group <- plot_data[[sample]]

    }
  }

  has_interval <- mode == "summary" && show_ribbon &&
    all(c(lower, upper) %in% names(plot_data))
  if (has_interval) {
    if (!is.numeric(plot_data[[lower]]) || !is.numeric(plot_data[[upper]]) ||
        anyNA(plot_data[[lower]]) || anyNA(plot_data[[upper]]) ||
        any(!is.finite(plot_data[[lower]])) ||
        any(!is.finite(plot_data[[upper]]))) {
      stop("The selected interval columns must contain finite numeric values.",
           call. = FALSE)
    }
    plot_data$.plot_lower <- plot_data[[lower]] * scale_factor
    plot_data$.plot_upper <- plot_data[[upper]] * scale_factor
    if (any(plot_data$.plot_lower > plot_data$.plot_upper)) {
      stop("Every lower interval value must be less than or equal to upper.",
           call. = FALSE)
    }
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$.plot_x))

  if (has_interval) {
    ribbon_defaults <- list(
      mapping = ggplot2::aes(
        ymin = .data$.plot_lower,
        ymax = .data$.plot_upper
      ),
      fill = ribbon_fill,
      alpha = ribbon_alpha,
      inherit.aes = TRUE
    )
    ribbon_defaults[names(ribbon_args)] <- ribbon_args
    p <- p + do.call(ggplot2::geom_ribbon, ribbon_defaults)
  }

  layer_args <- dots
  layer_args[names(line_args)] <- line_args

  if (mode == "samples") {
    mapping <- ggplot2::aes(
      y = .data$.plot_y,
      group = .data$.plot_group,
      colour = .data$.plot_sample
    )
    layer_args$mapping <- mapping
    layer_args$linewidth <- line_width
    layer_args$alpha <- line_alpha
    layer_args$inherit.aes <- TRUE

    if (smooth) {
      layer_args$method <- smooth_method
      layer_args$se <- smooth_se
      if (!is.null(smooth_formula)) layer_args$formula <- smooth_formula
      p <- p + do.call(ggplot2::geom_smooth, layer_args)
    } else {
      p <- p + do.call(ggplot2::geom_line, layer_args)
    }

    if (is.null(color_scale)) {
      color_scale <- ggplot2::scale_color_viridis_c(name = color_label)
    }
    p <- p + color_scale
  } else {
    mapping <- ggplot2::aes(y = .data$.plot_y)
    layer_args$mapping <- mapping
    layer_args$colour <- line_color
    layer_args$linewidth <- line_width
    layer_args$alpha <- line_alpha
    layer_args$inherit.aes <- TRUE

    if (smooth) {
      layer_args$method <- smooth_method
      layer_args$se <- smooth_se
      if (!is.null(smooth_formula)) layer_args$formula <- smooth_formula
      p <- p + do.call(ggplot2::geom_smooth, layer_args)
    } else {
      p <- p + do.call(ggplot2::geom_line, layer_args)
    }
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = x_limits,
      expand = x_expand
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      limits = y_limits,
      expand = y_expand
    )

  if (!is.null(facet)) {
    if (!requireNamespace("rlang", quietly = TRUE)) {
      stop("Package 'rlang' is required when facets are used.", call. = FALSE)
    }
    facet_formula <- as.formula(
      paste(
        "~",
        paste(facet, collapse = " + ")
      )
    )

    p <- p + ggplot2::facet_wrap(
      facet_formula,
      ncol = facet_ncol,
      nrow = facet_nrow
    )
  }

  if (is.null(x_label)) x_label <- tools::toTitleCase(gsub("_", " ", x))
  if (is.null(theme)) {

    theme <-
      ggplot2::theme_bw(
        base_size = base_size
      ) +
      ggplot2::theme(
        strip.background =
          ggplot2::element_rect(
            fill = "white",
            colour = "black"
          )
      )

  }

  p <- p +
    theme +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(face = text_face),
      axis.title = ggplot2::element_text(
        size = axis_title_size,
        face = text_face
      ),
      strip.background = strip_background,
      legend.position = legend_position
    )

  add_components(p, add)
}
