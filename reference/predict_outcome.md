# Predict outcome under user-defined exposure-lag profile(s)

Predicts the outcome associated with one or more user-defined
exposure-lag profiles, using the DLNM specification and template data
stored in the fitted model.

## Usage

``` r
predict_outcome(
  fit,
  profiles,
  re = c("population", "conditional"),
  id = NULL,
  allow_new_levels = FALSE,
  type = c("response", "link", "conditional"),
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000
)
```

## Arguments

- fit:

  Fitted model returned by \`fit_epidlnm()\`.

- profiles:

  Numeric vector (single exposure) or named list of numeric vectors. If
  the fitted model contains multiple exposure variables, \`profiles\`
  must be a named list matching those variables.

- re:

  Character. Prediction level: - \`"population"\`: excludes random
  effects - \`"conditional"\`: includes random effects where supported

- id:

  Optional vector of grouping levels (e.g., \`epi_id\`).

- allow_new_levels:

  Logical. Allow unseen grouping levels for mixed models.

- type:

  Prediction scale: \`"response"\`, \`"link"\`, or \`"conditional"\`.

- uncertainty:

  Logical. If \`TRUE\`, quantify uncertainty using simulated or
  posterior coefficient draws.

- output:

  Character. \`"summary"\` returns aggregated predictions; \`"samples"\`
  returns one row per sample (and per \`id\`, if provided).

- n_samples:

  Integer. Number of samples used for uncertainty quantification.

## Value

A data.frame.

\- If \`uncertainty = FALSE\`: same as before, with column
\`prediction\`.

\- If \`uncertainty = TRUE\` and \`output = "summary"\`: returns
columns: - \`prediction\` (median-based central estimate) - \`sd\` -
\`lower\` - \`upper\`

\- If \`uncertainty = TRUE\` and \`output = "samples"\`: returns one row
per sample (and per \`id\`, if provided), with columns: - \`sample\` -
\`prediction\`

## Details

This function supports both deterministic prediction and uncertainty
propagation. When \`uncertainty = TRUE\`, predictions are generated from
simulated or posterior draws of the model coefficients. If \`output =
"summary"\`, the central estimate is computed as the median of the
simulated predictions, while interval limits are obtained from empirical
quantiles (default: 2.5

\*\*Important:\*\* when \`uncertainty = TRUE\` and \`output =
"summary"\`, the column \`prediction\` represents the \*central
estimate\*, computed as the median of the predictive distribution.

Uncertainty is propagated using model-consistent sampling: - Bayesian
models (e.g., \`brms\`, \`INLA\`, \`bdlnm\`) use posterior draws -
Frequentist models use simulation from the asymptotic coefficient
distribution

For summary outputs under uncertainty, the median is used instead of the
mean to provide a more robust central estimate under asymmetric
predictive distributions.
