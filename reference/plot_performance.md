# Plot observed versus predicted model performance

Creates observed-versus-predicted performance plots for models evaluated
by \`find_bestfit()\`. The function retrieves out-of-fold predictions
from the \`"predictions_by_model"\` attribute and model-level validation
metrics from the main \`find_bestfit()\` result.

## Usage

``` r
plot_performance(
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
)
```

## Arguments

- object:

  A data frame returned by \`find_bestfit()\`. The object must contain a
  \`model_id\` column and the requested validation-metric columns, and
  must retain the \`"predictions_by_model"\` attribute.

- model_id:

  Optional vector of model identifiers to plot. If \`NULL\`, all models
  retained in \`object\` that have prediction data are plotted. Models
  are displayed in the order supplied by \`model_id\`, or in the current
  order of \`object\` when \`model_id = NULL\`.

- metrics:

  Character vector naming validation metrics stored as columns in
  \`object\`, such as \`"CCC"\`, \`"Cb"\`, \`"rho"\`, \`"RMSE"\`,
  \`"MAE"\`, \`"ROC_AUC"\`, or \`"Brier"\`. Metrics are printed in the
  supplied order.

- x:

  Character scalar naming the observed-value column in the prediction
  data. Default is \`"observed"\`.

- y:

  Character scalar naming the predicted-value column in the prediction
  data. Default is \`"predicted"\`.

- scale_factor:

  Finite numeric scalar multiplying both plotted axes. Use \`100\` to
  display proportions as percentages and \`1\` to retain the original
  response scale.

- ncol:

  Optional positive integer passed to \`ggplot2::facet_wrap()\`.

- metric_digits:

  Non-negative integer number of decimal places used in metric labels.

- x_lab:

  Optional x-axis label. If \`NULL\`, the name supplied to \`x\` is
  converted to title case.

- y_lab:

  Optional y-axis label. If \`NULL\`, the name supplied to \`y\` is
  converted to title case.

- show_fit:

  Logical. If \`TRUE\`, adds a fitted regression line using
  \`ggplot2::geom_smooth(method = "lm")\`.

- fit_se:

  Logical. Display the uncertainty band around the fitted line.

- fit_color:

  Color of the fitted regression line.

- show_identity:

  Logical. If \`TRUE\`, adds the 1:1 identity line.

- identity_color:

  Color of the identity line.

- identity_linetype:

  Line type of the identity line.

- label_x, label_y:

  Numeric positions for metric labels. Defaults are \`-Inf\` and
  \`Inf\`, placing labels at the upper-left corner of each panel.

- label_hjust, label_vjust:

  Horizontal and vertical justification for metric labels.

- label_size:

  Numeric text size for metric labels.

- base_size:

  Base font size passed to \`ggplot2::theme_bw()\`.

- point_args:

  Optional named list of additional arguments passed to
  \`ggplot2::geom_point()\`. This is useful when point customization
  should be supplied programmatically.

- ...:

  Additional arguments passed directly to \`ggplot2::geom_point()\`,
  such as \`size\`, \`shape\`, \`colour\`, \`alpha\`, \`fill\`, or
  \`stroke\`.

## Value

A \`ggplot\` object. The plot can be extended with additional ggplot2
layers, scales, labels, or themes using \`+\`.

## Details

One panel is produced per selected model. Validation metrics are
displayed inside each panel in the order supplied through \`metrics\`.
For example, \`metrics = c("CCC", "Cb", "rho")\` displays CCC first,
followed by Cb and rho on separate lines.

\`plot_performance()\` does not recalculate predictions or performance
metrics. Observed and predicted values are obtained from:

\`attr(object, "predictions_by_model")\`

The metric annotations are obtained from columns of \`object\`.
Consequently, the function preserves the exact out-of-fold values and
validation metrics calculated by \`find_bestfit()\`.

The default axes use the original response scale. For Beta or Binomial
outcomes stored as proportions, use \`scale_factor = 100\`, together
with axis labels such as \`x_lab = "Observed (

## Examples

``` r
if (FALSE) { # \dontrun{
plot_performance(
  best_models,
  model_id = best_models$model_id[1:6],
  metrics = c("CCC", "Cb", "rho"),
  scale_factor = 100,
  ncol = 3,
  x_lab = "Observed (%)",
  y_lab = "Predicted (%)",
  size = 2.5,
  alpha = 0.8
)

plot_performance(
  best_models,
  metrics = c("RMSE", "MAE"),
  point_args = list(shape = 21, fill = "steelblue", colour = "black")
) +
  ggplot2::theme(text = ggplot2::element_text(face = "bold"))
} # }
```
