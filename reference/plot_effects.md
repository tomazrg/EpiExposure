# Plot lag-specific or period-specific effects

Plots lag-specific effect surfaces or period-specific exposure-response
curves from an object returned by \`summarise_effects()\`.

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
  period_interval_linetype = 2,
  period_interval_linewidth = 0.5,
  period_show_legend = FALSE,
  panel_labels = c("(a)", "(b)"),
  label_size = 12,
  output = c("plot", "list")
)
```

## Arguments

- data:

  Data frame returned by \`summarise_effects()\`.

- scale:

  Plot scale. Use only \`"lag"\` or \`"period"\`.

- vars:

  Character vector containing variables to plot. \`NULL\` uses all
  available variables.

- metric:

  For lag plots, one or both of \`c("effect", "delta")\`.

- delta_multiplier:

  Multiplier applied to \`delta\` before plotting.

- metric_labels:

  Named labels for lag color legends.

- ylab:

  Lag y-axis labels. Prefer a named vector keyed by variable.

- metric_ncol:

  Number of metric blocks per row. Default \`1\` places effect above
  delta.

- vars_ncol:

  Number of variable panels per metric block. \`NULL\` uses all selected
  variables in one row.

- effect_palette:

  Three colors representing low, middle, and high values.

- delta_palette:

  \`NULL\` uses viridis; otherwise a vector of at least two colors.

- delta_viridis_option:

  Viridis option used when \`delta_palette = NULL\`.

- clip_quantiles:

  Quantiles used to clip lag fill values.

- effect_breaks:

  Contour breaks for \`effect\`.

- delta_bins:

  Number of contour bins for \`delta\`.

- lag_reverse:

  Logical; reverse the lag axis.

- lag_xlab:

  X-axis label for lag plots.

- legend_position:

  Named positions for effect and delta legends.

- lag_theme:

  One of \`"bw"\`, \`"minimal"\`, \`"classic"\`, or a ggplot2 theme.

- lag_base_size:

  Base font size for the lag theme.

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

  Passed to \`ggplot2::geom_raster()\`.

- period_order:

  Period facet order. \`NULL\` reverses names following the \`PeriodN\`
  convention.

- vars_order:

  Period variable order. \`NULL\` follows \`vars\`.

- period_response:

  Response displayed in period plots. Supported values are \`"eta"\` or
  \`"linear"\` for the linear-predictor contrast; \`"effect"\` for the
  effect column as stored; \`"exponentiated"\` (the alias
  \`"exponential"\` is also accepted) for an effect generated with
  \`effect_measure = "exponentiated"\`; \`"percent"\` for an effect
  generated with \`effect_measure = "percent"\`; and \`"baseline"\`,
  \`"predicted"\`, or \`"delta"\` for response-scale quantities. For
  summary input, corresponding columns named \`\<response\>\_lower\` and
  \`\<response\>\_upper\`, when present, are drawn as dashed uncertainty
  curves. Because exponentiated and percent results are both stored in
  \`effect\`, their aliases validate the \`epiexposure_effect_measure\`
  attribute created by \`summarise_effects()\`.

- period_exclude_samples:

  Samples removed from period draw-level plots. Default \`10\` preserves
  the historical display; use \`NULL\` to retain all samples.

- period_palette:

  \`NULL\` uses viridis; otherwise a vector of at least two colors.

- period_viridis_option:

  Viridis option for period curves.

- period_xlab:

  X-axis label for period plots.

- period_ylab:

  Y-axis label for period plots.

- period_theme:

  One of \`"bw"\`, \`"minimal"\`, \`"classic"\`, or a ggplot2 theme.

- period_base_size:

  Base size for the period theme.

- period_linewidth:

  Line width passed to \`ggplot2::geom_smooth()\`.

- period_interval_linetype:

  Line type used for lower and upper summary uncertainty curves.

- period_interval_linewidth:

  Line width used for lower and upper summary uncertainty curves.

- period_show_legend:

  Logical; show the period color legend.

- panel_labels:

  Character vector used to label lag-specific metric blocks.

- label_size:

  Positive finite size used for \`panel_labels\`.

- output:

  Either \`"plot"\` or \`"list"\`. For lag plots, \`"list"\` returns
  every component panel; for period plots, it returns \`list(period =
  plot)\`.

## Value

A ggplot/cowplot object, or a list of ggplots when \`output = "list"\`.

## Details

For period plots, draw-level input is recognized by the presence of a
\`sample\` column. One smoothed curve is then drawn for each parameter
draw. Otherwise, the data are treated as deterministic or summary
output. The selected central response is drawn as a solid black curve
and, when matching lower and upper columns are available, the
uncertainty limits are drawn as dashed black curves. The bounds are used
exactly as returned by \`summarise_effects()\`; \`plot_effects()\` does
not recalculate uncertainty.
