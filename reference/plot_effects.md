# Plot lag-specific or period-specific effects

Plot lag-specific or period-specific effects

## Usage

``` r
plot_effects(
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
)
```

## Arguments

- data:

  Data frame returned by summarise_effects().

- scale:

  Plot scale. Use only \\lag\\ or \\period\\.

- vars:

  Character vector with variables to plot. NULL uses all variables.

- metric:

  For lag plots, one or both of c("effect", "delta").

- delta_multiplier:

  Multiplier applied to delta before plotting.

- metric_labels:

  Named labels for the lag color legends.

- ylab:

  Lag y-axis labels. Prefer a named vector keyed by variable.

- metric_ncol:

  Number of metric blocks per row. Default 1 places effect above delta,
  matching the manual figure.

- vars_ncol:

  Number of variable panels per metric block. NULL uses all selected
  variables in one row.

- effect_palette:

  Three colors: low, mid, high.

- delta_palette:

  NULL uses viridis; otherwise a vector of \>= 2 colors.

- delta_viridis_option:

  Viridis option used when delta_palette is NULL.

- clip_quantiles:

  Quantiles used to clip lag fill values.

- effect_breaks:

  Contour breaks for effect.

- delta_bins:

  Number of contour bins for delta.

- lag_reverse:

  Reverse the lag axis.

- lag_xlab:

  X-axis label for lag plots.

- legend_position:

  Named positions for effect and delta legends.

- lag_theme:

  One of "bw", "minimal", "classic", or a ggplot2 theme.

- lag_base_size:

  Base size for lag theme.

- axis_title_size:

  Axis-title size for lag plots.

- axis_text_size:

  Axis-text size for lag plots.

- contour_colour:

  Contour color.

- contour_linewidth:

  Contour line width.

- zero_linewidth:

  Delta zero-contour line width.

- effect_zero_text_size:

  Effect zero-contour label size.

- delta_zero_text_size:

  Delta zero-contour label size.

- interpolate:

  Passed to geom_raster().

- period_order:

  Period facet order. NULL reproduces reversed PeriodN order.

- vars_order:

  Period variable order. NULL follows \`vars\`.

- period_response:

  Response column for period plots; default "eta".

- period_exclude_samples:

  Samples removed from period plot. Default 10 reproduces the supplied
  manual plot; use NULL to keep every sample.

- period_palette:

  NULL uses viridis; otherwise vector of \>= 2 colors.

- period_viridis_option:

  Viridis option for period curves.

- period_xlab:

  X-axis label for period plots.

- period_ylab:

  Y-axis label for period plots.

- period_theme:

  One of "bw", "minimal", "classic", or a ggplot2 theme.

- period_base_size:

  Base size for period theme. Default 11 matches theme_bw().

- period_linewidth:

  Line width passed to geom_smooth().

- period_show_legend:

  Whether to show the period color legend.

- panel_labels:

  Character vector used to label the lag-specific metric blocks in the
  combined plot. Default is \`c("(a)", "(b)")\`, so the first selected
  metric block is labelled "(a)" and the second "(b)". Labels are
  applied only when \`scale = "lag"\` and \`output = "plot"\`; \`output
  = "list"\` returns the original unlabelled component plots.

- label_size:

  Positive finite size used for \`panel_labels\`. Default is \`12\`.

- output:

  Either "plot" or "list". For lag, "list" returns every panel.

## Value

A ggplot/cowplot object, or a list of ggplots when output = "list".
