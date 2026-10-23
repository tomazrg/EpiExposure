# Plot out-of-fold model performance

Creates performance plots for models evaluated by \`find_bestfit()\`.
For non-binary outcomes, the default is an observed-versus-predicted
plot. For binomial outcomes, the default is a receiver operating
characteristic (ROC) curve.

## Usage

``` r
plot_performance(
  object,
  model_id = NULL,
  metrics = "CCC",
  type = c("auto", "observed_predicted", "roc", "threshold", "confusion"),
  threshold = NULL,
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

  A non-empty data frame returned by \`find_bestfit()\`. It must contain
  a unique \`model_id\` column, the requested metric columns, and
  retained out-of-fold predictions.

- model_id:

  Optional vector of model identifiers to plot. If \`NULL\`, every
  retained model with available predictions is plotted.

- metrics:

  Unique character vector naming validation metrics stored as numeric
  columns of \`object\`. Metrics are displayed in the supplied order.
  Common non-binary metrics include \`"CCC"\`, \`"Cb"\`, \`"rho"\`,
  \`"RMSE"\`, and \`"MAE"\`. Binomial metrics include \`"ROC_AUC"\`,
  \`"Brier"\`, \`"LogLoss"\`, \`"Accuracy"\`, \`"Balanced_Accuracy"\`,
  \`"Sensitivity"\`, \`"Specificity"\`, \`"Precision"\`, \`"F1"\`, and
  \`"MCC"\`.

- type:

  Character. Plot type selected from \`"auto"\`,
  \`"observed_predicted"\`, \`"roc"\`, \`"threshold"\`, or
  \`"confusion"\`. \`"auto"\` uses an ROC curve for binomial outcomes
  and an observed-versus-predicted plot otherwise.

- threshold:

  Optional probability threshold strictly between 0 and 1. If \`NULL\`,
  \`attr(object, "threshold")\` is used when available; otherwise
  \`0.5\` is used. It is used by \`"threshold"\` and \`"confusion"\`
  plots and identifies the selected operating point on ROC curves.

- x:

  Character scalar naming the observed-value column in prediction data.
  Default is \`"observed"\`.

- y:

  Character scalar naming the predicted-value column in prediction data.
  For binomial outcomes, this column must contain probabilities in
  \`\[0, 1\]\`. Default is \`"predicted"\`.

- scale_factor:

  Finite non-zero numeric scalar multiplying both axes only for \`type =
  "observed_predicted"\`. Use \`100\` to display proportions as
  percentages.

- ncol:

  Optional positive integer passed to \`ggplot2::facet_wrap()\`.

- metric_digits:

  Non-negative integer controlling metric-label precision.

- x_lab, y_lab:

  Optional axis labels. Plot-type-specific labels are used when either
  argument is \`NULL\`.

- show_fit:

  Logical. For observed-versus-predicted plots, add a linear fitted
  line.

- fit_se:

  Logical. Display an uncertainty band around the fitted line.

- fit_color:

  Color of the fitted regression line.

- show_identity:

  Logical. For observed-versus-predicted plots, add the 1:1 line. For
  ROC curves, add the no-discrimination diagonal.

- identity_color:

  Color of the identity or no-discrimination line.

- identity_linetype:

  Line type of the identity or no-discrimination line.

- label_x, label_y:

  Numeric positions for metric labels. Infinite defaults place labels
  close to a panel corner where supported.

- label_hjust, label_vjust:

  Horizontal and vertical justification for metric labels.

- label_size:

  Numeric text size for metric labels and confusion-matrix counts.

- base_size:

  Base font size passed to \`ggplot2::theme_bw()\`.

- point_args:

  Optional named list of additional arguments passed to
  \`ggplot2::geom_point()\` for observed-versus-predicted plots.

- ...:

  Additional named arguments passed to \`ggplot2::geom_point()\` for
  observed-versus-predicted plots.

## Value

A \`ggplot\` object that can be extended with additional ggplot2 layers,
scales, labels, or themes.

## Details

The function retrieves out-of-fold predictions from the
\`"predictions_by_model"\` attribute, or from the \`"predictions"\`
attribute when the model-specific list is unavailable. Model-level
validation metrics are retrieved from columns of the main
\`find_bestfit()\` result.

\`plot_performance()\` does not refit models or recalculate their stored
validation metrics. Predictions are obtained from:

\`attr(object, "predictions_by_model")\`

or, when that list is unavailable:

\`attr(object, "predictions")\`

For binomial outcomes, ROC curves are reconstructed directly from the
observed classes and out-of-fold predicted probabilities. Sensitivity
and specificity are calculated over the distinct probability thresholds
represented in the selected predictions.

The ROC curve uses false-positive rate (\`1 - specificity\`) on the
x-axis and sensitivity on the y-axis. A highlighted point identifies
performance at \`threshold\`.

Threshold plots display sensitivity and specificity across probability
thresholds and add a vertical line at the selected threshold.

Confusion matrices classify probabilities using \`predicted \>=
threshold\` and report classification counts.

The ROC and threshold plots require both observed classes. The observed
values must be coded as \`0\` and \`1\`, and predicted values must be
finite probabilities in \`\[0, 1\]\`.

## Examples

``` r
performance_results <- data.frame(
  model_id = c("model_1", "model_2"),
  CCC = c(0.91, 0.86),
  Cb = c(0.96, 0.93),
  rho = c(0.95, 0.92),
  RMSE = c(0.07, 0.10),
  MAE = c(0.05, 0.08)
)

continuous_predictions <- list(
  model_1 = data.frame(
    observed = c(0.05, 0.12, 0.20, 0.31, 0.43, 0.55),
    predicted = c(0.06, 0.14, 0.18, 0.29, 0.46, 0.53)
  ),
  model_2 = data.frame(
    observed = c(0.05, 0.12, 0.20, 0.31, 0.43, 0.55),
    predicted = c(0.09, 0.16, 0.24, 0.27, 0.39, 0.49)
  )
)

attr(
  performance_results,
  "predictions_by_model"
) <- continuous_predictions

continuous_plot <- plot_performance(
  performance_results,
  metrics = c("CCC", "Cb", "rho"),
  scale_factor = 100,
  x_lab = "Observed (%)",
  y_lab = "Predicted (%)",
  size = 2.5,
  alpha = 0.8
)

binary_results <- data.frame(
  model_id = c("model_1", "model_2"),
  ROC_AUC = c(0.88, 0.82),
  Brier = c(0.14, 0.18)
)

binary_predictions <- list(
  model_1 = data.frame(
    observed = c(0, 0, 0, 1, 1, 1),
    predicted = c(0.08, 0.21, 0.39, 0.58, 0.77, 0.91)
  ),
  model_2 = data.frame(
    observed = c(0, 0, 0, 1, 1, 1),
    predicted = c(0.18, 0.31, 0.48, 0.52, 0.68, 0.79)
  )
)

attr(binary_results, "outcome_type") <- "binary"
attr(binary_results, "threshold") <- 0.5
attr(binary_results, "predictions_by_model") <- binary_predictions

roc_plot <- plot_performance(
  binary_results,
  metrics = c("ROC_AUC", "Brier")
)

threshold_plot <- plot_performance(
  binary_results,
  model_id = "model_1",
  metrics = "ROC_AUC",
  type = "threshold"
)

confusion_plot <- plot_performance(
  binary_results,
  model_id = "model_1",
  metrics = c("ROC_AUC", "Brier"),
  type = "confusion"
)

stopifnot(
  inherits(continuous_plot, "ggplot"),
  inherits(roc_plot, "ggplot"),
  inherits(threshold_plot, "ggplot"),
  inherits(confusion_plot, "ggplot")
)
```
