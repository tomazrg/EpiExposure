#' Plot simulated yield and economic losses
#'
#' Creates either coordinated heatmaps of yield and economic losses or a
#' regression-style visualization from an object returned by `simulate_losses()`.
#'
#' With `type = "heatmap"` (default), the function preserves the original
#' raster-plot workflow. Yield loss and economic loss can be displayed alone or
#' together, and scenario variables can be shown as facets.
#'
#' With `type = "regression"`, `plots` is not used and a single regression plot
#' is returned. The response plotted on the y-axis is selected with `reg_y`.
#' Optionally, `att_yield` filters the data to one attainable-yield scenario
#' before plotting. A `yield_model` facet is generated automatically from
#' `intercept_sim` and `slope_sim`.
#'
#' Regression input is recognized as draw-level when a `sample` column is
#' present. In that case, one `geom_smooth(se = FALSE)` curve is drawn for each
#' sample using `group = sample` and `color = sample`. Otherwise, the data are
#' treated as summary output. The central response is drawn as a solid black
#' smooth and, when `<reg_y>_lower` and `<reg_y>_upper` are present, the lower
#' and upper curves are drawn as dashed black smooths.
#'
#' @param data A non-empty data frame returned by `simulate_losses()`.
#' @param x Character scalar naming the variable displayed on the x-axis.
#' @param type Character scalar selecting `"heatmap"` or `"regression"`.
#'   Default is `"heatmap"`.
#' @param y Character scalar naming the y-axis variable for heatmaps.
#'   Default is `"att_yield"`, the shortened output name returned by the current
#'   `simulate_losses()`.
#' @param yield_loss Character scalar naming the yield-loss column used by the
#'   heatmap. Default is `"yield_loss"`.
#' @param economic_loss Character scalar naming the economic-loss column used by
#'   the heatmap. Default is `"econ_loss"`.
#' @param reg_y Character scalar naming the y-axis variable when
#'   `type = "regression"`. Default is `"ap_yield"`.
#' @param att_yield Optional finite numeric scalar used only with
#'   `type = "regression"`. When supplied, rows are filtered to the selected
#'   `att_yield` value before plotting.
#' @param facet Optional character vector naming scenario variables used as
#'   facets. For heatmaps, `facet_wrap()` is used as in the original function.
#'   For regression plots, the first facet variable is placed in rows and any
#'   remaining facet variables are placed in columns together with the
#'   automatically generated `yield_model` facet. Thus
#'   `facet = c("tmean", "wetness")` reproduces
#'   `facet_grid(rows = vars(tmean), cols = vars(wetness, yield_model))`.
#' @param plots Character vector selecting `"yield"`, `"economic"`, or both.
#'   Used only when `type = "heatmap"` and ignored for regression plots.
#' @param ncol,nrow Optional positive integers controlling the arrangement of
#'   the main heatmap plots. When both are `NULL` and both plots are requested,
#'   `nrow = 2` is used.
#' @param facet_ncol,facet_nrow Optional positive integers controlling the
#'   arrangement of heatmap facets.
#' @param yield_fill_scale,economic_fill_scale Optional continuous ggplot2 fill
#'   scales for heatmaps. When `NULL`, `ggplot2::scale_fill_viridis_b()` is used.
#' @param yield_fill_label,economic_fill_label Legend titles used by the default
#'   heatmap fill scales.
#' @param x_breaks,y_breaks Optional numeric axis breaks.
#' @param x_expand,y_expand Numeric expansion vectors passed to the axes.
#' @param x_limits,y_limits Optional numeric vectors of length two.
#' @param x_label Axis label. If `NULL`, a title-case label is generated from
#'   `x`.
#' @param y_label Y-axis label used for heatmaps.
#' @param reg_y_label Optional y-axis label used for regression plots. If `NULL`,
#'   a label is generated from `reg_y`, with informative defaults for standard
#'   `simulate_losses()` output columns.
#' @param yield_title,economic_title Optional heatmap panel titles.
#' @param show_yield_x_text Logical. Show x-axis text and ticks in the yield
#'   heatmap panel. Default is `FALSE`.
#' @param show_yield_x_label Logical. Show the x-axis title in the yield heatmap
#'   panel. Default is `FALSE`.
#' @param show_economic_strips Logical. Show facet strips in the economic
#'   heatmap panel. Default is `FALSE`.
#' @param raster_interpolate Logical passed to `ggplot2::geom_raster()`.
#' @param raster_args Named list of arguments passed to both raster layers.
#' @param yield_raster_args,economic_raster_args Named lists of panel-specific
#'   raster arguments. Panel-specific arguments override `raster_args` and `...`.
#' @param theme A ggplot2 theme. The default is
#'   `ggplot2::theme_bw(base_size = base_size)`.
#' @param base_size Base font size used by the default theme.
#' @param text_face Font face used for common heatmap text.
#' @param axis_title_size Axis-title size.
#' @param strip_background Facet-strip background in heatmaps and regression
#'   plots.
#' @param legend_position Legend position used in heatmaps. Regression plots
#'   intentionally suppress the legend, matching the requested display.
#' @param tag_levels Tag sequence passed to `patchwork::plot_annotation()`.
#' @param tag_prefix,tag_suffix Prefix and suffix for plot tags.
#' @param tag_size,tag_face,tag_position Appearance and position of plot tags.
#' @param collect_guides Logical. Collect compatible guides with patchwork.
#' @param yield_add,economic_add Optional ggplot2 component or list of
#'   components added to the corresponding heatmap panel.
#' @param combined_add Optional patchwork-compatible component added to the
#'   combined heatmap.
#' @param ... Additional named arguments passed to both heatmap raster layers.
#'
#' @return With `type = "regression"`, a `ggplot` object. With
#'   `type = "heatmap"`, a `patchwork` object when both heatmaps are requested or
#'   a `ggplot` object when only one is requested. Heatmap panels are retained as
#'   `"yield_plot"` and `"economic_plot"` attributes.
#'
#' @details
#' ## Regression mode
#'
#' For summary-style data, the function plots the central column selected by
#' `reg_y`. If columns named `<reg_y>_lower` and `<reg_y>_upper` are available,
#' they are plotted using dashed black `geom_smooth(se = FALSE)` curves. These
#' bounds are used exactly as returned by `simulate_losses()`; no interval is
#' recalculated by `plot_losses()`.
#'
#' For draw-level data containing `sample`, the function plots one smooth curve
#' per sample and maps both `group` and `color` to `sample`. No summary interval
#' is calculated.
#'
#' @examples
#' \dontrun{
#' # Heatmaps using the shortened simulate_losses() output names
#' plot_losses(
#'   losses,
#'   x = "rain",
#'   type = "heatmap",
#'   facet = c("tmean", "wetness"),
#'   x_breaks = seq(0, 15, by = 3),
#'   y_breaks = seq(4000, 12000, by = 2000),
#'   x_label = "Daily precipitation (mm)",
#'   nrow = 2,
#'   facet_ncol = 4
#' )
#'
#' # Regression from summary-style output
#' plot_losses(
#'   losses,
#'   x = "rain",
#'   type = "regression",
#'   reg_y = "ap_yield",
#'   att_yield = 9000,
#'   facet = c("tmean", "wetness")
#' )
#'
#' # Regression from draw-level output containing `sample`
#' plot_losses(
#'   losses_samples,
#'   x = "rain",
#'   type = "regression",
#'   reg_y = "ap_yield",
#'   att_yield = 9000,
#'   facet = c("tmean", "wetness")
#' )
#' }
#'
#' @importFrom rlang .data
#' @export
plot_losses <- function(
    data,
    x,
    type = c("heatmap", "regression"),
    y = "att_yield",
    yield_loss = "yield_loss",
    economic_loss = "econ_loss",
    reg_y = "ap_yield",
    att_yield = NULL,
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
    reg_y_label = NULL,
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

  # HELPERS

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required by `plot_losses()`.", call. = FALSE)
  }

  scalar_name <- function(z) {
    is.character(z) &&
      length(z) == 1L &&
      !is.na(z) &&
      nzchar(z)
  }

  scalar_flag <- function(z) {
    is.logical(z) &&
      length(z) == 1L &&
      !is.na(z)
  }

  scalar_number <- function(z) {
    is.numeric(z) &&
      length(z) == 1L &&
      !is.na(z) &&
      is.finite(z)
  }

  whole_scalar <- function(z) {
    scalar_number(z) &&
      z == as.integer(z)
  }

  validate_optional_integer <- function(z, argument) {
    if (!is.null(z) &&
        (!whole_scalar(z) || z < 1L)) {
      stop(
        "`", argument, "` must be NULL or a positive integer.",
        call. = FALSE
      )
    }

    if (is.null(z)) NULL else as.integer(z)
  }

  validate_numeric <- function(
    z,
    argument,
    length_two = FALSE
  ) {
    if (is.null(z)) {
      return(NULL)
    }

    if (!is.numeric(z) ||
        !length(z) ||
        anyNA(z) ||
        any(!is.finite(z))) {
      stop(
        "`", argument, "` must contain finite numeric values.",
        call. = FALSE
      )
    }

    if (length_two &&
        length(z) != 2L) {
      stop(
        "`", argument, "` must contain exactly two values.",
        call. = FALSE
      )
    }

    z
  }

  validate_named_list <- function(z, argument) {
    if (!is.list(z)) {
      stop(
        "`", argument, "` must be a list.",
        call. = FALSE
      )
    }

    if (length(z) &&
        (is.null(names(z)) ||
         anyNA(names(z)) ||
         any(!nzchar(names(z))))) {
      stop(
        "Every entry in `", argument, "` must be named.",
        call. = FALSE
      )
    }

    z
  }

  validate_fill_scale <- function(z, argument) {
    if (!is.null(z) &&
        !inherits(z, "Scale")) {
      stop(
        "`", argument, "` must be NULL or a ggplot2 fill scale.",
        call. = FALSE
      )
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
    if (is.null(components)) {
      return(plot_object)
    }

    if (!is.list(components) ||
        inherits(
          components,
          c("gg", "ggplot", "Scale", "theme", "Coord")
        )) {
      components <- list(components)
    }

    for (component in components) {
      plot_object <- plot_object + component
    }

    plot_object
  }

  automatic_regression_label <- function(column) {
    known <- c(
      ap_yield = "Attainable predicted yield (kg/ha)",
      pred_yield = "Predicted yield",
      yield_loss = "Yield loss (kg/ha)",
      econ_loss = "Economic loss (USD/ha)",
      yl_prop_raw = "Raw yield-loss proportion",
      yl_prop = "Yield-loss proportion",
      yl_pct = "Yield loss (%)",
      rl_prop = "Relative-yield proportion",
      rl_pct = "Relative yield (%)"
    )

    if (column %in% names(known)) {
      return(unname(known[[column]]))
    }

    tools::toTitleCase(
      gsub("_", " ", column)
    )
  }

  make_regression_labeller <- function(facet_names) {
    if (!length(facet_names)) {
      return(ggplot2::labeller())
    }

    labellers <- lapply(
      facet_names,
      function(variable) {
        force(variable)
        function(values) {
          paste0(variable, " = ", values)
        }
      }
    )
    names(labellers) <- facet_names

    do.call(
      ggplot2::labeller,
      labellers
    )
  }

  # GENERAL VALIDATION

  if (!is.data.frame(data) ||
      !nrow(data)) {
    stop(
      "`data` must be a non-empty data frame returned by `simulate_losses()`.",
      call. = FALSE
    )
  }

  if (!scalar_name(x)) {
    stop(
      "`x` must be one non-empty column name.",
      call. = FALSE
    )
  }

  type <- match.arg(type)

  if (!is.null(facet) &&
      (!is.character(facet) ||
       !length(facet) ||
       anyNA(facet) ||
       any(!nzchar(facet)) ||
       anyDuplicated(facet))) {
    stop(
      "`facet` must be NULL or unique non-empty column names.",
      call. = FALSE
    )
  }

  if (!is.null(facet) &&
      !requireNamespace("rlang", quietly = TRUE)) {
    stop(
      "Package 'rlang' is required when facets are requested.",
      call. = FALSE
    )
  }

  x_breaks <- validate_numeric(
    x_breaks,
    "x_breaks"
  )
  y_breaks <- validate_numeric(
    y_breaks,
    "y_breaks"
  )
  x_expand <- validate_numeric(
    x_expand,
    "x_expand"
  )
  y_expand <- validate_numeric(
    y_expand,
    "y_expand"
  )
  x_limits <- validate_numeric(
    x_limits,
    "x_limits",
    TRUE
  )
  y_limits <- validate_numeric(
    y_limits,
    "y_limits",
    TRUE
  )

  if (!scalar_number(base_size) ||
      base_size <= 0) {
    stop(
      "`base_size` must be one positive finite number.",
      call. = FALSE
    )
  }

  if (!scalar_number(axis_title_size) ||
      axis_title_size <= 0) {
    stop(
      "`axis_title_size` must be one positive finite number.",
      call. = FALSE
    )
  }

  if (is.null(x_label)) {
    x_label <- tools::toTitleCase(
      gsub("_", " ", x)
    )
  }

  if (is.null(theme)) {
    theme <- ggplot2::theme_bw(
      base_size = base_size
    )
  }

  # REGRESSION MODE

  if (identical(type, "regression")) {

    if (!requireNamespace("rlang", quietly = TRUE)) {
      stop(
        "Package 'rlang' is required for regression faceting.",
        call. = FALSE
      )
    }

    if (!scalar_name(reg_y)) {
      stop(
        "`reg_y` must be one non-empty column name when ",
        "`type = 'regression'`.",
        call. = FALSE
      )
    }

    required_regression <- unique(
      c(
        x,
        reg_y,
        "att_yield",
        "intercept_sim",
        "slope_sim",
        facet
      )
    )

    missing_regression <- setdiff(
      required_regression,
      names(data)
    )

    if (length(missing_regression)) {
      stop(
        "Column(s) missing from `data` for regression plotting: ",
        paste(
          missing_regression,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    numeric_regression <- unique(
      c(
        x,
        reg_y,
        "att_yield",
        "intercept_sim",
        "slope_sim"
      )
    )

    invalid_regression <- numeric_regression[
      !vapply(
        data[numeric_regression],
        function(z) {
          is.numeric(z) &&
            !anyNA(z) &&
            all(is.finite(z))
        },
        logical(1)
      )
    ]

    if (length(invalid_regression)) {
      stop(
        "Regression plotting columns must contain finite numeric values: ",
        paste(
          invalid_regression,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    if (!is.null(att_yield) &&
        !scalar_number(att_yield)) {
      stop(
        "`att_yield` must be NULL or one finite numeric value.",
        call. = FALSE
      )
    }

    regression_data <- data

    if (!is.null(att_yield)) {
      tolerance <- sqrt(.Machine$double.eps) *
        pmax(
          1,
          abs(regression_data$att_yield),
          abs(att_yield)
        )

      keep <- abs(
        regression_data$att_yield -
          att_yield
      ) <= tolerance

      if (!any(keep)) {
        stop(
          "No rows were found for `att_yield = ",
          format(att_yield, digits = 10),
          "`.",
          call. = FALSE
        )
      }

      regression_data <- regression_data[
        keep,
        ,
        drop = FALSE
      ]
    }

    regression_data$yield_model <- paste0(
      "B0 = ",
      regression_data$intercept_sim,
      " | B1 = ",
      regression_data$slope_sim
    )

    regression_is_samples <-
      "sample" %in% names(regression_data)

    lower_column <- paste0(
      reg_y,
      "_lower"
    )
    upper_column <- paste0(
      reg_y,
      "_upper"
    )

    has_lower <- lower_column %in%
      names(regression_data)
    has_upper <- upper_column %in%
      names(regression_data)

    if (!regression_is_samples &&
        xor(has_lower, has_upper)) {
      stop(
        "Summary regression data contain only one interval endpoint for `",
        reg_y,
        "`. Both `",
        lower_column,
        "` and `",
        upper_column,
        "` are required to plot interval curves.",
        call. = FALSE
      )
    }

    if (!regression_is_samples &&
        has_lower &&
        has_upper) {
      invalid_bounds <- c(
        lower_column,
        upper_column
      )[
        !vapply(
          regression_data[
            c(
              lower_column,
              upper_column
            )
          ],
          function(z) {
            is.numeric(z) &&
              !anyNA(z) &&
              all(is.finite(z))
          },
          logical(1)
        )
      ]

      if (length(invalid_bounds)) {
        stop(
          "Regression interval columns must contain finite numeric values: ",
          paste(
            invalid_bounds,
            collapse = ", "
          ),
          ".",
          call. = FALSE
        )
      }
    }

    if (regression_is_samples) {
      if (anyNA(regression_data$sample)) {
        stop(
          "`sample` cannot contain missing values in regression draw-level data.",
          call. = FALSE
        )
      }

      regression_plot <- ggplot2::ggplot(
        regression_data
      ) +
        ggplot2::geom_smooth(
          ggplot2::aes(
            x = .data[[x]],
            y = .data[[reg_y]],
            group = .data$sample,
            color = .data$sample
          ),
          se = FALSE
        )

      if (is.numeric(regression_data$sample)) {
        regression_plot <-
          regression_plot +
          ggplot2::scale_color_viridis_c()
      } else {
        regression_plot <-
          regression_plot +
          ggplot2::scale_color_viridis_d()
      }

    } else {

      regression_plot <- ggplot2::ggplot(
        regression_data
      ) +
        ggplot2::geom_smooth(
          ggplot2::aes(
            x = .data[[x]],
            y = .data[[reg_y]]
          ),
          se = FALSE,
          color = "black"
        )

      if (has_lower &&
          has_upper) {
        regression_plot <-
          regression_plot +
          ggplot2::geom_smooth(
            ggplot2::aes(
              x = .data[[x]],
              y = .data[[lower_column]]
            ),
            se = FALSE,
            color = "black",
            linetype = 2,
            linewidth = 0.5
          ) +
          ggplot2::geom_smooth(
            ggplot2::aes(
              x = .data[[x]],
              y = .data[[upper_column]]
            ),
            se = FALSE,
            color = "black",
            linetype = 2,
            linewidth = 0.5
          )
      }
    }

    row_facets <- if (length(facet)) {
      facet[1L]
    } else {
      character(0)
    }

    column_facets <- c(
      if (length(facet) > 1L) {
        facet[-1L]
      } else {
        character(0)
      },
      "yield_model"
    )

    rows_spec <- if (length(row_facets)) {
      ggplot2::vars(
        !!!rlang::syms(
          row_facets
        )
      )
    } else {
      ggplot2::vars()
    }

    cols_spec <- ggplot2::vars(
      !!!rlang::syms(
        column_facets
      )
    )

    regression_plot <-
      regression_plot +
      ggplot2::facet_grid(
        rows = rows_spec,
        cols = cols_spec,
        labeller = make_regression_labeller(
          facet
        )
      ) +
      ggplot2::scale_x_continuous(
        breaks = if (is.null(x_breaks)) {ggplot2::waiver()} else {x_breaks},
        limits = x_limits,
        expand = x_expand
      ) +
      ggplot2::scale_y_continuous(
        breaks = if(is.null(y_breaks)) {ggplot2::waiver()} else {y_breaks},
        limits = y_limits,
        expand = y_expand
      )

    if (is.null(reg_y_label)) {
      reg_y_label <-
        automatic_regression_label(
          reg_y
        )
    }

    regression_plot <-
      regression_plot +
      theme +
      ggplot2::labs(
        x = x_label,
        y = reg_y_label
      ) +
      ggplot2::theme(
        legend.position = "none",
        strip.background = strip_background,
        strip.text = ggplot2::element_text(
          face = "bold"
        )
      )

    attr(
      regression_plot,
      "epiexposure_plot_type"
    ) <- "regression"

    attr(
      regression_plot,
      "epiexposure_regression_output"
    ) <- if (regression_is_samples) {
      "samples"
    } else {
      "summary"
    }

    attr(
      regression_plot,
      "epiexposure_regression_y"
    ) <- reg_y

    attr(
      regression_plot,
      "epiexposure_att_yield_filter"
    ) <- att_yield

    return(
      regression_plot
    )
  }

  # HEATMAP MODE

  for (argument in c(
    "y",
    "yield_loss",
    "economic_loss"
  )) {
    if (!scalar_name(
      get(
        argument,
        inherits = FALSE
      )
    )) {
      stop(
        "`", argument,
        "` must be one non-empty column name.",
        call. = FALSE
      )
    }
  }

  plots <- match.arg(
    plots,
    c(
      "yield",
      "economic"
    ),
    several.ok = TRUE
  )

  if (anyDuplicated(plots)) {
    stop(
      "`plots` must not contain duplicates.",
      call. = FALSE
    )
  }

  if (length(plots) > 1L &&
      !requireNamespace(
        "patchwork",
        quietly = TRUE
      )) {
    stop(
      "Package 'patchwork' is required when both plots are requested.",
      call. = FALSE
    )
  }

  required <- unique(
    c(
      x,
      y,
      facet,
      if ("yield" %in% plots) {
        yield_loss
      },
      if ("economic" %in% plots) {
        economic_loss
      }
    )
  )

  missing <- setdiff(
    required,
    names(data)
  )

  if (length(missing)) {
    stop(
      "Column(s) missing from `data`: ",
      paste(
        missing,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  numeric_columns <- unique(
    c(
      x,
      y,
      if ("yield" %in% plots) {
        yield_loss
      },
      if ("economic" %in% plots) {
        economic_loss
      }
    )
  )

  invalid <- numeric_columns[
    !vapply(
      data[numeric_columns],
      function(z) {
        is.numeric(z) &&
          !anyNA(z) &&
          all(is.finite(z))
      },
      logical(1)
    )
  ]

  if (length(invalid)) {
    stop(
      "Plotting columns must contain finite numeric values: ",
      paste(
        invalid,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  ncol <- validate_optional_integer(
    ncol,
    "ncol"
  )
  nrow <- validate_optional_integer(
    nrow,
    "nrow"
  )
  facet_ncol <- validate_optional_integer(
    facet_ncol,
    "facet_ncol"
  )
  facet_nrow <- validate_optional_integer(
    facet_nrow,
    "facet_nrow"
  )

  for (argument in c(
    "show_yield_x_text",
    "show_yield_x_label",
    "show_economic_strips",
    "raster_interpolate",
    "collect_guides"
  )) {
    if (!scalar_flag(
      get(
        argument,
        inherits = FALSE
      )
    )) {
      stop(
        "`", argument,
        "` must be TRUE or FALSE.",
        call. = FALSE
      )
    }
  }

  raster_args <- validate_named_list(
    raster_args,
    "raster_args"
  )
  yield_raster_args <- validate_named_list(
    yield_raster_args,
    "yield_raster_args"
  )
  economic_raster_args <- validate_named_list(
    economic_raster_args,
    "economic_raster_args"
  )
  dots <- validate_named_list(
    list(...),
    "..."
  )

  yield_fill_scale <- validate_fill_scale(
    yield_fill_scale,
    "yield_fill_scale"
  )
  economic_fill_scale <- validate_fill_scale(
    economic_fill_scale,
    "economic_fill_scale"
  )

  if (is.null(yield_fill_scale)) {
    yield_fill_scale <-
      ggplot2::scale_fill_viridis_b(
        name = yield_fill_label
      )
  }

  if (is.null(economic_fill_scale)) {
    economic_fill_scale <-
      ggplot2::scale_fill_viridis_b(
        name = economic_fill_label
      )
  }

  facet_layer <- if (is.null(facet)) {
    NULL
  } else {
    ggplot2::facet_wrap(
      ggplot2::vars(
        !!!rlang::syms(
          facet
        )
      ),
      ncol = facet_ncol,
      nrow = facet_nrow
    )
  }

  yield_geom_args <- merge_args(
    raster_args,
    dots,
    yield_raster_args
  )
  economic_geom_args <- merge_args(
    raster_args,
    dots,
    economic_raster_args
  )

  make_plot <- function(
    fill_column,
    fill_scale,
    title,
    geom_args
  ) {
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[x]],
        y = .data[[y]],
        fill = .data[[fill_column]]
      )
    ) +
      do.call(
        ggplot2::geom_raster,
        geom_args
      ) +
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        limits = x_limits,
        expand = x_expand
      ) +
      ggplot2::scale_y_continuous(
        breaks = y_breaks,
        limits = y_limits,
        expand = y_expand
      ) +
      fill_scale +
      theme +
      ggplot2::labs(
        x = x_label,
        y = y_label,
        title = title
      ) +
      ggplot2::theme(
        text = ggplot2::element_text(
          face = text_face
        ),
        axis.title = ggplot2::element_text(
          size = axis_title_size,
          face = text_face
        ),
        legend.position = legend_position
      )

    if (!is.null(facet_layer)) {
      p <- p +
        facet_layer
    }

    p
  }

  yield_plot <- NULL
  economic_plot <- NULL

  if ("yield" %in% plots) {
    yield_plot <- make_plot(
      yield_loss,
      yield_fill_scale,
      yield_title,
      yield_geom_args
    ) +
      ggplot2::theme(
        strip.background = strip_background,
        axis.text.x = if (show_yield_x_text) {
          ggplot2::element_text()
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.x = if (show_yield_x_text) {
          ggplot2::element_line()
        } else {
          ggplot2::element_blank()
        },
        axis.title.x = if (show_yield_x_label) {
          ggplot2::element_text(
            size = axis_title_size,
            face = text_face
          )
        } else {
          ggplot2::element_blank()
        }
      )

    yield_plot <- add_components(
      yield_plot,
      yield_add
    )
  }

  if ("economic" %in% plots) {
    economic_plot <- make_plot(
      economic_loss,
      economic_fill_scale,
      economic_title,
      economic_geom_args
    )

    if (!show_economic_strips) {
      economic_plot <-
        economic_plot +
        ggplot2::theme(
          strip.background =
            ggplot2::element_blank(),
          strip.text =
            ggplot2::element_blank()
        )
    }

    economic_plot <- add_components(
      economic_plot,
      economic_add
    )
  }

  selected <- list()

  if ("yield" %in% plots) {
    selected$yield <- yield_plot
  }

  if ("economic" %in% plots) {
    selected$economic <- economic_plot
  }

  if (length(selected) == 1L) {
    result <- selected[[1L]]

    attr(
      result,
      "yield_plot"
    ) <- yield_plot

    attr(
      result,
      "economic_plot"
    ) <- economic_plot

    attr(
      result,
      "epiexposure_plot_type"
    ) <- "heatmap"

    return(
      result
    )
  }

  if (is.null(nrow) &&
      is.null(ncol)) {
    nrow <- 2L
  }

  result <- patchwork::wrap_plots(
    selected,
    nrow = nrow,
    ncol = ncol,
    guides = if (collect_guides) {
      "collect"
    } else {
      "keep"
    }
  ) +
    patchwork::plot_annotation(
      tag_levels = tag_levels,
      tag_prefix = tag_prefix,
      tag_suffix = tag_suffix
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        size = tag_size,
        face = tag_face
      ),
      plot.tag.position = tag_position
    )

  if (!is.null(combined_add)) {
    result <- result +
      combined_add
  }

  attr(
    result,
    "yield_plot"
  ) <- yield_plot

  attr(
    result,
    "economic_plot"
  ) <- economic_plot

  attr(
    result,
    "epiexposure_plot_type"
  ) <- "heatmap"

  result
}
