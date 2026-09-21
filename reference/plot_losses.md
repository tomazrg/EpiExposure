# Plot simulated yield and economic losses

Creates either coordinated heatmaps of yield and economic losses or a
regression-style visualization from an object returned by
\`simulate_losses()\`.

## Usage

``` r
plot_losses(
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
  yield_fill_label = expression(bold(kg ~ ha^{
-1
 })),
  economic_fill_label = expression(bold(USD ~ ha^{
-1
 })),
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
)
```

## Arguments

- data:

  A non-empty data frame returned by \`simulate_losses()\`.

- x:

  Character scalar naming the variable displayed on the x-axis.

- type:

  Character scalar selecting \`"heatmap"\` or \`"regression"\`. Default
  is \`"heatmap"\`.

- y:

  Character scalar naming the y-axis variable for heatmaps. Default is
  \`"att_yield"\`, the shortened output name returned by the current
  \`simulate_losses()\`.

- yield_loss:

  Character scalar naming the yield-loss column used by the heatmap.
  Default is \`"yield_loss"\`.

- economic_loss:

  Character scalar naming the economic-loss column used by the heatmap.
  Default is \`"econ_loss"\`.

- reg_y:

  Character scalar naming the y-axis variable when \`type =
  "regression"\`. Default is \`"ap_yield"\`.

- att_yield:

  Optional finite numeric scalar used only with \`type = "regression"\`.
  When supplied, rows are filtered to the selected \`att_yield\` value
  before plotting.

- facet:

  Optional character vector naming scenario variables used as facets.
  For heatmaps, \`facet_wrap()\` is used as in the historical function.
  For regression plots, the first facet variable is placed in rows and
  any remaining facet variables are placed in columns together with the
  automatically generated \`yield_model\` facet. Thus \`facet =
  c("tmean", "wetness")\` reproduces \`facet_grid(rows = vars(tmean),
  cols = vars(wetness, yield_model))\`.

- plots:

  Character vector selecting \`"yield"\`, \`"economic"\`, or both. Used
  only when \`type = "heatmap"\` and ignored for regression plots.

- ncol, nrow:

  Optional positive integers controlling the arrangement of the main
  heatmap plots. When both are \`NULL\` and both plots are requested,
  \`nrow = 2\` is used.

- facet_ncol, facet_nrow:

  Optional positive integers controlling the arrangement of heatmap
  facets.

- yield_fill_scale, economic_fill_scale:

  Optional continuous ggplot2 fill scales for heatmaps. When \`NULL\`,
  \`ggplot2::scale_fill_viridis_b()\` is used.

- yield_fill_label, economic_fill_label:

  Legend titles used by the default heatmap fill scales.

- x_breaks, y_breaks:

  Optional numeric axis breaks.

- x_expand, y_expand:

  Numeric expansion vectors passed to the axes.

- x_limits, y_limits:

  Optional numeric vectors of length two.

- x_label:

  Axis label. If \`NULL\`, a title-case label is generated from \`x\`.

- y_label:

  Y-axis label used for heatmaps.

- reg_y_label:

  Optional y-axis label used for regression plots. If \`NULL\`, a label
  is generated from \`reg_y\`, with informative defaults for standard
  \`simulate_losses()\` output columns.

- yield_title, economic_title:

  Optional heatmap panel titles.

- show_yield_x_text:

  Logical. Show x-axis text and ticks in the yield heatmap panel.
  Default is \`FALSE\`.

- show_yield_x_label:

  Logical. Show the x-axis title in the yield heatmap panel. Default is
  \`FALSE\`.

- show_economic_strips:

  Logical. Show facet strips in the economic heatmap panel. Default is
  \`FALSE\`.

- raster_interpolate:

  Logical passed to \`ggplot2::geom_raster()\`.

- raster_args:

  Named list of arguments passed to both raster layers.

- yield_raster_args, economic_raster_args:

  Named lists of panel-specific raster arguments. Panel-specific
  arguments override \`raster_args\` and \`...\`.

- theme:

  A ggplot2 theme. The default is \`ggplot2::theme_bw(base_size =
  base_size)\`.

- base_size:

  Base font size used by the default theme.

- text_face:

  Font face used for common heatmap text.

- axis_title_size:

  Axis-title size.

- strip_background:

  Facet-strip background in heatmaps and regression plots.

- legend_position:

  Legend position used in heatmaps. Regression plots intentionally
  suppress the legend, matching the requested display.

- tag_levels:

  Tag sequence passed to \`patchwork::plot_annotation()\`.

- tag_prefix, tag_suffix:

  Prefix and suffix for plot tags.

- tag_size, tag_face, tag_position:

  Appearance and position of plot tags.

- collect_guides:

  Logical. Collect compatible guides with patchwork.

- yield_add, economic_add:

  Optional ggplot2 component or list of components added to the
  corresponding heatmap panel.

- combined_add:

  Optional patchwork-compatible component added to the combined heatmap.

- ...:

  Additional named arguments passed to both heatmap raster layers.

## Value

With \`type = "regression"\`, a \`ggplot\` object. With \`type =
"heatmap"\`, a \`patchwork\` object when both heatmaps are requested or
a \`ggplot\` object when only one is requested. Heatmap panels are
retained as \`"yield_plot"\` and \`"economic_plot"\` attributes.

## Details

With \`type = "heatmap"\` (default), the function preserves the
historical raster-plot workflow. Yield loss and economic loss can be
displayed alone or together, and scenario variables can be shown as
facets.

With \`type = "regression"\`, \`plots\` is not used and a single
regression plot is returned. The response plotted on the y-axis is
selected with \`reg_y\`. Optionally, \`att_yield\` filters the data to
one attainable-yield scenario before plotting. A \`yield_model\` facet
is generated automatically from \`intercept_sim\` and \`slope_sim\`.

Regression input is recognized as draw-level when a \`sample\` column is
present. In that case, one \`geom_smooth(se = FALSE)\` curve is drawn
for each sample using \`group = sample\` and \`color = sample\`.
Otherwise, the data are treated as summary output. The central response
is drawn as a solid black smooth and, when \`\<reg_y\>\_lower\` and
\`\<reg_y\>\_upper\` are present, the lower and upper curves are drawn
as dashed black smooths.

\## Regression mode

For summary-style data, the function plots the central column selected
by \`reg_y\`. If columns named \`\<reg_y\>\_lower\` and
\`\<reg_y\>\_upper\` are available, they are plotted using dashed black
\`geom_smooth(se = FALSE)\` curves. These bounds are used exactly as
returned by \`simulate_losses()\`; no interval is recalculated by
\`plot_losses()\`.

For draw-level data containing \`sample\`, the function plots one smooth
curve per sample and maps both \`group\` and \`color\` to \`sample\`. No
summary interval is calculated.

## Examples

``` r
if (FALSE) { # \dontrun{
# Heatmaps using the shortened simulate_losses() output names
plot_losses(
  losses,
  x = "rain",
  type = "heatmap",
  facet = c("tmean", "wetness"),
  x_breaks = seq(0, 15, by = 3),
  y_breaks = seq(4000, 12000, by = 2000),
  x_label = "Daily precipitation (mm)",
  nrow = 2,
  facet_ncol = 4
)

# Regression from summary-style output
plot_losses(
  losses,
  x = "rain",
  type = "regression",
  reg_y = "ap_yield",
  att_yield = 9000,
  facet = c("tmean", "wetness")
)

# Regression from draw-level output containing `sample`
plot_losses(
  losses_samples,
  x = "rain",
  type = "regression",
  reg_y = "ap_yield",
  att_yield = 9000,
  facet = c("tmean", "wetness")
)
} # }
```
