# Compare predicted outcomes between multiple exposure scenarios

Compares predicted outcomes across two or more exposure scenarios, using
\`predict_outcome()\` as the computational backend.

## Usage

``` r
compare_predictions(
  fit,
  profiles = NULL,
  profiles1 = NULL,
  profiles2 = NULL,
  re = c("population", "conditional"),
  id = NULL,
  allow_new_levels = FALSE,
  type = c("response", "link", "conditional"),
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000,
  eps = 1e-12
)
```

## Arguments

- fit:

  Fitted model (output of \`fit_epidlnm()\`).

- profiles:

  Can be: - a named list of profiles (recommended), e.g.:
  \`list(scenario1 = ..., scenario2 = ...)\`, or - \`NULL\` when using
  \`profiles1\` and \`profiles2\` for backward compatibility.

- profiles1:

  (legacy) First profile (used only if \`profiles\` is NULL).

- profiles2:

  (legacy) Second profile (used only if \`profiles\` is NULL).

- re:

  Character. Prediction level: - \`"population"\`: excludes random
  effects (default). - \`"conditional"\`: includes random effects where
  supported.

- id:

  Optional character string indicating the column used as identifier.

- allow_new_levels:

  Logical. Passed to \`predict_outcome()\`.

- type:

  Character. Scale of prediction: \`"response"\` (default), \`"link"\`,
  or \`"conditional"\`.

- uncertainty:

  Logical. If \`TRUE\`, uncertainty is propagated using sample-based
  predictions.

- output:

  Character. Output type when \`uncertainty = TRUE\`: - \`"summary"\`:
  returns median-based summaries (default) - \`"samples"\`: returns all
  simulated samples

- n_samples:

  Integer. Number of samples used for uncertainty propagation.

- eps:

  Small positive constant used for numerical stability.

## Value

A data.frame with pairwise comparisons including: - \`scenario1\`,
\`scenario2\` - \`pred1\`, \`pred2\` - \`diff\` - \`percent_change\` -
\`ratio\`

When \`uncertainty = TRUE\`: - \`"samples"\`: returns sample-level
comparisons - \`"summary"\`: returns median-based estimates, standard
deviation, and empirical interval limits

## Details

This function supports both deterministic comparisons and uncertainty
propagation. When \`uncertainty = TRUE\`, predictions are computed at
the sample level and comparisons are derived from these simulated
values.

If \`output = "summary"\`, the central estimate is computed as the
median of the sample-based distributions, and interval limits are
derived from empirical quantiles (default: 2.5

\*\*Important:\*\* although traditional terminology might suggest
"mean", all central estimates in summary outputs correspond to the
\*median\* when uncertainty is propagated, ensuring robustness under
asymmetric distributions.

When \`uncertainty = TRUE\`, this function always operates on simulated
prediction samples obtained from \`predict_outcome(output =
"samples")\`.

Summaries are then computed as: - central estimate: median - uncertainty
intervals: empirical quantiles

This approach ensures coherent uncertainty propagation for both
frequentist (simulation-based) and Bayesian (posterior-based) models,
avoiding incorrect analytic variance approximations.
