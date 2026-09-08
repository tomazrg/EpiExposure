#' Plot simulated yield and economic losses
#'
#' Creates coordinated raster plots of yield and economic losses from an object
#' returned by `simulate_losses()`. By default, yield loss is shown above
#' economic loss, both use `ggplot2::scale_fill_viridis_b()`, and scenario
#' variables can be displayed as facets.
#'
#' @param data A non-empty data frame returned by `simulate_losses()`.
#' @param x Character scalar naming the variable displayed on the x-axis.
#' @param y Character scalar naming the variable displayed on the y-axis.
#'   Default is `"attainable_yield"`.
#' @param yield_loss Character scalar naming the yield-loss column.
#' @param economic_loss Character scalar naming the economic-loss column.
#' @param facet Optional character vector naming facet variables.
#' @param plots Character vector selecting `"yield"`, `"economic"`, or both.
#' @param ncol,nrow Optional positive integers controlling the arrangement of
#'   the main loss plots. When both are `NULL` and both plots are requested,
#'   `nrow = 2` is used.
#' @param facet_ncol,facet_nrow Optional positive integers controlling the
#'   arrangement of facets within each plot.
#' @param yield_fill_scale,economic_fill_scale Optional continuous ggplot2 fill
#'   scales. When `NULL`, `ggplot2::scale_fill_viridis_b()` is used. Any
#'   compatible scale may be supplied, including `scale_fill_gradient()`,
#'   `scale_fill_gradientn()`, `scale_fill_distiller()`, or a customized
#'   `scale_fill_viridis_b()`.
#' @param yield_fill_label,economic_fill_label Legend titles used by the default
#'   fill scales. Defaults are `kg ha^-1` and `USD ha^-1`, respectively.
#'   Character values and plotmath expressions are accepted.
#' @param x_breaks,y_breaks Optional numeric axis breaks.
#' @param x_expand,y_expand Numeric expansion vectors passed to the axes.
#' @param x_limits,y_limits Optional numeric vectors of length two.
#' @param x_label,y_label Axis labels. If `x_label = NULL`, a title-case label
#'   is generated from `x`.
#' @param yield_title,economic_title Optional panel titles.
#' @param show_yield_x_text Logical. Show x-axis text and ticks in the yield
#'   panel. Default is `FALSE`.
#' @param show_yield_x_label Logical. Show the x-axis title in the yield panel.
#'   Default is `FALSE`.
#' @param show_economic_strips Logical. Show facet strips in the economic panel.
#'   Default is `FALSE`.
#' @param raster_interpolate Logical passed to `ggplot2::geom_raster()`.
#' @param raster_args Named list of arguments passed to both raster layers.
#' @param yield_raster_args,economic_raster_args Named lists of panel-specific
#'   raster arguments. Panel-specific arguments override `raster_args` and `...`.
#' @param theme A ggplot2 theme applied to both plots. The default is
#'   `ggplot2::theme_bw(base_size = base_size)`.
#' @param base_size Base font size used by the default theme.
#' @param text_face Font face used for common plot text.
#' @param axis_title_size Axis-title size.
#' @param strip_background Facet-strip background in the yield panel.
#' @param legend_position Legend position used in both panels.
#' @param tag_levels Tag sequence passed to `patchwork::plot_annotation()`.
#' @param tag_prefix,tag_suffix Prefix and suffix for plot tags.
#' @param tag_size,tag_face,tag_position Appearance and position of plot tags.
#' @param collect_guides Logical. Collect compatible guides with patchwork.
#' @param yield_add,economic_add Optional ggplot2 component or list of
#'   components added to the corresponding panel.
#' @param combined_add Optional patchwork-compatible component added to the
#'   combined plot.
#' @param ... Additional named arguments passed to both raster layers.
#'
#' @return A `patchwork` object when both plots are requested, or a `ggplot`
#'   object when only one plot is requested. Individual panels are retained as
#'   `"yield_plot"` and `"economic_plot"` attributes.
#'
#' @details
#' Custom fill scales are used exactly as supplied. Therefore, when passing a
#' custom scale, define its `name` argument when a custom legend title is needed.
#' The default scales receive `yield_fill_label` and `economic_fill_label`
#' automatically.
#'
#' @examples
#' \dontrun{
#' plot_losses(
#'   losses,
#'   x = "rain",
#'   facet = c("tmean", "wetness"),
#'   x_breaks = seq(0, 15, by = 3),
#'   y_breaks = seq(4000, 12000, by = 2000),
#'   x_label = "Daily precipitation (mm)",
#'   nrow = 2,
#'   facet_ncol = 4
#' )
#'
#' plot_losses(
#'   losses,
#'   x = "rain",
#'   facet = c("tmean", "wetness"),
#'   yield_fill_scale = ggplot2::scale_fill_gradientn(
#'     colours = c("white", "gold", "firebrick"),
#'     name = expression(kg ~ ha^{-1})
#'   ),
#'   economic_fill_scale = ggplot2::scale_fill_viridis_b(
#'     option = "magma",
#'     name = expression(USD ~ ha^{-1})
#'   )
#' )
#' }
#'
#' @export
plot_losses <- function(
    data,
    x,
    y = "attainable_yield",
    yield_loss = "yield_loss",
    economic_loss = "economic_loss",
    facet = NULL,
    plots = c("yield", "economic"),
    ncol = NULL,
    nrow = NULL,
    facet_ncol = 4L,
    facet_nrow = NULL,
    yield_fill_scale = NULL,
    economic_fill_scale = NULL,
    yield_fill_label = expression(bold(kg ~ ha^{-1})),
    economic_fill_label = expression(bold(USD ~ ha^{-1})),
    x_breaks = NULL,
    y_breaks = NULL,
    x_expand = c(0.02, 0.02),
    y_expand = c(0.02, 0.02),
    x_limits = NULL,
    y_limits = NULL,
    x_label = NULL,
    y_label = "Attainable yield (kg/ha)",
    yield_title = NULL,
    economic_title = NULL,
    show_yield_x_text = FALSE,
    show_yield_x_label = FALSE,
    show_economic_strips = FALSE,
    raster_interpolate = FALSE,
    raster_args = list(),
    yield_raster_args = list(),
    economic_raster_args = list(),
    theme = NULL,
    base_size = 10,
    text_face = "bold",
    axis_title_size = 12,
    strip_background = ggplot2::element_rect(fill = "white", colour = "black"),
    legend_position = "right",
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")",
    tag_size = 12,
    tag_face = "bold",
    tag_position = c(0.01, 1.02),
    collect_guides = FALSE,
    yield_add = NULL,
    economic_add = NULL,
    combined_add = NULL,
    ...
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required by `plot_losses()`.", call. = FALSE)
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
    if (is.null(z)) return(NULL)
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
  validate_fill_scale <- function(z, argument) {
    if (!is.null(z) && !inherits(z, "Scale")) {
      stop("`", argument, "` must be NULL or a ggplot2 fill scale.",
           call. = FALSE)
    }
    z
  }
  merge_args <- function(common, dots, specific) {
    out <- common
    out[names(dots)] <- dots
    out[names(specific)] <- specific
    out$interpolate <- raster_interpolate
    out
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
    stop("`data` must be a non-empty data frame returned by `simulate_losses()`.",
         call. = FALSE)
  }
  for (argument in c("x", "y", "yield_loss", "economic_loss")) {
    if (!scalar_name(get(argument, inherits = FALSE))) {
      stop("`", argument, "` must be one non-empty column name.", call. = FALSE)
    }
  }
  if (!is.null(facet) &&
      (!is.character(facet) || !length(facet) || anyNA(facet) ||
       any(!nzchar(facet)) || anyDuplicated(facet))) {
    stop("`facet` must be NULL or unique non-empty column names.", call. = FALSE)
  }
  
  plots <- match.arg(plots, c("yield", "economic"), several.ok = TRUE)
  if (anyDuplicated(plots)) stop("`plots` must not contain duplicates.", call. = FALSE)
  if (length(plots) > 1L && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required when both plots are requested.",
         call. = FALSE)
  }
  if (!is.null(facet) && !requireNamespace("rlang", quietly = TRUE)) {
    stop("Package 'rlang' is required when facets are requested.", call. = FALSE)
  }
  
  required <- unique(c(
    x, y, facet,
    if ("yield" %in% plots) yield_loss,
    if ("economic" %in% plots) economic_loss
  ))
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Column(s) missing from `data`: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  
  numeric_columns <- unique(c(
    x, y,
    if ("yield" %in% plots) yield_loss,
    if ("economic" %in% plots) economic_loss
  ))
  invalid <- numeric_columns[!vapply(data[numeric_columns], function(z) {
    is.numeric(z) && !anyNA(z) && all(is.finite(z))
  }, logical(1))]
  if (length(invalid)) {
    stop("Plotting columns must contain finite numeric values: ",
         paste(invalid, collapse = ", "), ".", call. = FALSE)
  }
  
  ncol <- validate_optional_integer(ncol, "ncol")
  nrow <- validate_optional_integer(nrow, "nrow")
  facet_ncol <- validate_optional_integer(facet_ncol, "facet_ncol")
  facet_nrow <- validate_optional_integer(facet_nrow, "facet_nrow")
  x_breaks <- validate_numeric(x_breaks, "x_breaks")
  y_breaks <- validate_numeric(y_breaks, "y_breaks")
  x_expand <- validate_numeric(x_expand, "x_expand")
  y_expand <- validate_numeric(y_expand, "y_expand")
  x_limits <- validate_numeric(x_limits, "x_limits", TRUE)
  y_limits <- validate_numeric(y_limits, "y_limits", TRUE)
  
  for (argument in c(
    "show_yield_x_text", "show_yield_x_label", "show_economic_strips",
    "raster_interpolate", "collect_guides"
  )) {
    if (!scalar_flag(get(argument, inherits = FALSE))) {
      stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  
  raster_args <- validate_named_list(raster_args, "raster_args")
  yield_raster_args <- validate_named_list(yield_raster_args, "yield_raster_args")
  economic_raster_args <- validate_named_list(
    economic_raster_args, "economic_raster_args"
  )
  dots <- validate_named_list(list(...), "...")
  
  yield_fill_scale <- validate_fill_scale(yield_fill_scale, "yield_fill_scale")
  economic_fill_scale <- validate_fill_scale(
    economic_fill_scale, "economic_fill_scale"
  )
  if (is.null(yield_fill_scale)) {
    yield_fill_scale <- ggplot2::scale_fill_viridis_b(name = yield_fill_label)
  }
  if (is.null(economic_fill_scale)) {
    economic_fill_scale <- ggplot2::scale_fill_viridis_b(
      name = economic_fill_label
    )
  }
  
  if (is.null(x_label)) x_label <- tools::toTitleCase(gsub("_", " ", x))
  if (is.null(theme)) theme <- ggplot2::theme_bw(base_size = base_size)
  
  facet_layer <- if (is.null(facet)) NULL else ggplot2::facet_wrap(
    ggplot2::vars(!!!rlang::syms(facet)),
    ncol = facet_ncol,
    nrow = facet_nrow
  )
  
  yield_geom_args <- merge_args(raster_args, dots, yield_raster_args)
  economic_geom_args <- merge_args(raster_args, dots, economic_raster_args)
  
  make_plot <- function(fill_column, fill_scale, title, geom_args) {
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[x]], y = .data[[y]], fill = .data[[fill_column]]
      )
    ) +
      do.call(ggplot2::geom_raster, geom_args) +
      ggplot2::scale_x_continuous(
        breaks = x_breaks, limits = x_limits, expand = x_expand
      ) +
      ggplot2::scale_y_continuous(
        breaks = y_breaks, limits = y_limits, expand = y_expand
      ) +
      fill_scale +
      theme +
      ggplot2::labs(x = x_label, y = y_label, title = title) +
      ggplot2::theme(
        text = ggplot2::element_text(face = text_face),
        axis.title = ggplot2::element_text(
          size = axis_title_size, face = text_face
        ),
        legend.position = legend_position
      )
    if (!is.null(facet_layer)) p <- p + facet_layer
    p
  }
  
  yield_plot <- NULL
  economic_plot <- NULL
  
  if ("yield" %in% plots) {
    yield_plot <- make_plot(
      yield_loss, yield_fill_scale, yield_title, yield_geom_args
    ) + ggplot2::theme(
      strip.background = strip_background,
      axis.text.x = if (show_yield_x_text) {
        ggplot2::element_text()
      } else ggplot2::element_blank(),
      axis.ticks.x = if (show_yield_x_text) {
        ggplot2::element_line()
      } else ggplot2::element_blank(),
      axis.title.x = if (show_yield_x_label) {
        ggplot2::element_text(size = axis_title_size, face = text_face)
      } else ggplot2::element_blank()
    )
    yield_plot <- add_components(yield_plot, yield_add)
  }
  
  if ("economic" %in% plots) {
    economic_plot <- make_plot(
      economic_loss, economic_fill_scale, economic_title, economic_geom_args
    )
    if (!show_economic_strips) {
      economic_plot <- economic_plot + ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank()
      )
    }
    economic_plot <- add_components(economic_plot, economic_add)
  }
  
  selected <- list()
  if ("yield" %in% plots) selected$yield <- yield_plot
  if ("economic" %in% plots) selected$economic <- economic_plot
  
  if (length(selected) == 1L) {
    result <- selected[[1L]]
    attr(result, "yield_plot") <- yield_plot
    attr(result, "economic_plot") <- economic_plot
    return(result)
  }
  
  if (is.null(nrow) && is.null(ncol)) nrow <- 2L
  result <- patchwork::wrap_plots(
    selected,
    nrow = nrow,
    ncol = ncol,
    guides = if (collect_guides) "collect" else "keep"
  ) +
    patchwork::plot_annotation(
      tag_levels = tag_levels,
      tag_prefix = tag_prefix,
      tag_suffix = tag_suffix
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = tag_size, face = tag_face),
      plot.tag.position = tag_position
    )
  
  if (!is.null(combined_add)) result <- result + combined_add
  attr(result, "yield_plot") <- yield_plot
  attr(result, "economic_plot") <- economic_plot
  result
}
