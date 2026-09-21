# Compare period-specific DLNM effects

Compares period-specific DLNM associations returned by
\`summarise_effects(scale = "period")\` against one selected reference
period.

## Usage

``` r
compare_periods(
  period_df,
  period_ref = 1,
  periods = NULL,
  output = c("summary", "samples"),
  interval_probs = NULL,
  eps = 1e-10
)
```

## Arguments

- period_df:

  Output from the current \`summarise_effects(scale = "period")\`. A
  data frame is expected. If an \`"epiexposure_effects"\` object
  produced with diagnostics is supplied, its \`\$effects\` component is
  used. The object must carry the current EpiExposure temporal metadata:
  one common fitted \`max_lag\`, exact \`history_length = max_lag + 1\`,
  and the exact-history contract.

  Current EpiExposure period output must contain \`var\`, \`value\`,
  \`period\`, \`scale\`, \`eta\`, \`effect\`, \`baseline\`,
  \`predicted\`, and \`delta\`.

  If \`summarise_effects(..., uncertainty = TRUE)\` was used,
  \`period_df\` must be the \*\*sample-level\*\* output produced with
  \`output = "samples"\`. Comparing already-summarized medians and
  interval bounds would not preserve covariance among period effects.

- period_ref:

  Reference period. Either:

  \- one positive integer selecting the period by its defined order;
  or - one exact non-empty character period label, for example \`"W1"\`.

  Numeric references do not require period labels to end in numbers. If
  \`periods\` is supplied, its row order defines the numeric period
  order. Otherwise the first-occurrence order in \`period_df\$period\`
  is used.

- periods:

  Optional period-definition data frame, typically returned by
  \`define_periods()\`, containing \`period\`, \`lag_start\`, and
  \`lag_end\`.

  When supplied, period labels are validated against \`period_df\`,
  every period boundary must lie within the inherited common fitted lag
  range \`0:max_lag\`, and \`lag_start\`, \`lag_end\`, and \`n_lags\`
  are appended to the comparison output. The number of lags is

  \`n_lags = lag_end - lag_start + 1\`.

  Supplying \`periods\` is recommended when the compared periods have
  unequal widths because cumulative DLNM effects are sums across
  included lags.

- output:

  Character. \`"summary"\` or \`"samples"\`.

  For deterministic input (\`uncertainty = FALSE\`), only \`"summary"\`
  is available and contains the direct period comparisons.

  For sample-level uncertain input, \`"samples"\` returns draw-by-draw
  comparisons and \`"summary"\` returns medians, SDs, and empirical
  intervals calculated \*\*after\*\* the draw-by-draw comparison.

- interval_probs:

  Optional numeric vector of length two defining the empirical interval
  for uncertain comparisons. If \`NULL\`, the interval stored by
  \`summarise_effects()\` is reused; otherwise the default \`c(0.025,
  0.975)\` is used when no stored interval is available.

- eps:

  Positive finite numerical tolerance. It is used only for consistency
  checks of values that should be equal under the EpiExposure contract;
  it is not used to alter effects or to force a denominator away from
  zero.

## Value

A data frame containing the original period-specific quantities and the
corresponding reference-period quantities. Principal comparison columns
are:

- \`eta_diff\`:

  Difference in additive DLNM contrasts:

  \`eta - ref_eta\`.

  This is the canonical comparison on the fitted link scale.

- \`diff\`:

  Arithmetic difference in the displayed \`effect\`:

  \`effect - ref_effect\`.

  Its units depend on the \`effect_measure\` used by
  \`summarise_effects()\`.

- \`ratio\`:

  For \`log\` links, the ratio of the two period-specific response
  ratios. For \`logit\` links, the ratio of the two period-specific odds
  ratios:

  \`exp(eta - ref_eta)\`.

  For other links, no general multiplicative interpretation is defined,
  so \`ratio\` and \`ratio_percent\` are \`NA\`.

- \`ratio_percent\`:

  \`100 \* (ratio - 1)\` when \`ratio\` is defined.

- \`predicted_diff\`:

  Difference between the period-specific population expected responses:

  \`predicted - ref_predicted\`.

- \`delta_diff\`:

  Difference between period-specific absolute response-scale changes:

  \`delta - ref_delta\`.

  Because the joint \`baseline\` is identical across compared periods,
  \`delta_diff\` and \`predicted_diff\` should be numerically equal.

With uncertain sample-level input and \`output = "summary"\`, each
numeric quantity is summarized by its median plus \`\_sd\`, \`\_lower\`,
and \`\_upper\` columns.

Output attributes record the reference period, fitted link, effect
measure, prediction contract, uncertainty contract, ratio
interpretation, common fitted \`max_lag\`, exact history length, and the
inherited EpiExposure exact-history contract.

## Details

The comparison preserves the EpiExposure effect contract:

\- \`eta\` is the fundamental additive DLNM contrast on the fitted link
scale; - \`effect\` is the user-requested representation of \`eta\`
(\`"linear"\`, \`"exponentiated"\`, or \`"percent"\`); - \`predicted\`
is the population-level expected outcome for that period; - \`delta =
predicted - baseline\` is the absolute response-scale change from the
joint exposure reference profile.

Period comparisons are always made at the same exposure variable and the
same exposure evaluation value.

\## Exact common lag/history contract

\`compare_periods()\` does not receive the original long-format exposure
histories, so it cannot recount time rows itself. Instead, it requires
the temporal metadata propagated by the current \`summarise_effects()\`
output:

\$\$ history\\length = max\\lag + 1. \$\$

The stored contract must state that all fitted exposure variables used
the same exact lag window. This prevents period comparisons from
silently mixing effect objects created under incompatible temporal
histories.

When \`periods\` is supplied, its lag boundaries are additionally
checked against the inherited common fitted \`max_lag\`.
\`compare_periods()\` does not require a manually supplied \`periods\`
object to cover every fitted lag, because analysts may intentionally
compare a subset of defined periods; it only requires every supplied
period to remain inside the fitted lag window.

\## What is being compared?

\`summarise_effects(scale = "period")\` first sums lag-specific DLNM
contributions on the additive linear-predictor scale within each
requested period. \`compare_periods()\` then compares those
period-specific cumulative contrasts at the \*\*same exposure value\*\*.

If period \\j\\ has contrast \\\eta_j(x)\\ and the selected reference
period has contrast \\\eta_r(x)\\, the fundamental comparison is

\$\$\Delta\eta\_{j:r}(x) = \eta_j(x) - \eta_r(x).\$\$

Under a log link,

\$\$\exp\\\Delta\eta\_{j:r}(x)\\ = RR_j(x) / RR_r(x),\$\$

while under a logit link,

\$\$\exp\\\Delta\eta\_{j:r}(x)\\ = OR_j(x) / OR_r(x).\$\$

Therefore the returned multiplicative \`ratio\` is derived from
\`eta_diff\`, not from blindly dividing whichever transformed \`effect\`
happens to be displayed.

This distinction matters for \`effect_measure = "percent"\`. For
example, a +20 percent relative effect and a +10 percent relative effect
do \*\*not\*\* imply a period ratio of \`20 / 10 = 2\`. Their underlying
multiplicative effects are 1.20 and 1.10, and the appropriate relative
ratio is \`1.20 / 1.10\`.

\## Response-scale comparison

\`predicted_diff\` compares expected outcomes under the same joint
reference profile while changing the focal exposure over different lag
periods. For Beta and Binomial outcomes, multiplying \`predicted_diff\`
by 100 gives the difference in percentage points between the two
period-specific expected proportions/probabilities.

Since \`delta = predicted - baseline\` and the same baseline is used for
every period comparison at a given parameter draw, \`delta_diff\` must
equal \`predicted_diff\` up to numerical tolerance. The function
validates this identity.

\## Unequal period widths

Period effects are cumulative sums over lags. Consequently, comparing a
21-lag period with a 10-lag period compares the cumulative associations
as defined; the contrast can reflect both the lag-response pattern and
the different number of included lags. \`compare_periods()\` does not
divide cumulative effects by period width. Use period definitions with
the scientific interpretation intended by the analysis.

\## Uncertainty

When \`summarise_effects()\` uses uncertainty, all requested effects are
based on the same joint parameter draw. \`compare_periods()\` preserves
this covariance by matching the target and reference periods \*\*within
the same draw\*\*, calculating \`eta_diff\`, \`diff\`, \`ratio\`,
\`predicted_diff\`, and \`delta_diff\` draw by draw, and only then
summarizing those comparisons.

Already summarized uncertain output is therefore rejected explicitly.
Subtracting medians or combining marginal interval endpoints would not
recover the distribution of a difference or ratio.
