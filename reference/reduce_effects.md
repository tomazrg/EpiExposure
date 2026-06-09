# Reduce DLNM effects to one dimension (article-consistent)

Reduces a DLNM exposure-lag-response surface into a one-dimensional
summary using \`dlnm::crossreduce()\`, with optional uncertainty
propagation.

## Usage

``` r
reduce_effects(
  fit,
  wx_long,
  var,
  lag_max = NULL,
  df_var = NULL,
  df_lag = NULL,
  fun_var = NULL,
  fun_lag = NULL,
  type = c("overall", "lag", "var"),
  value = NULL,
  scale = c("percent", "response", "link"),
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000
)
```

## Arguments

- fit:

  Fitted model object.

- wx_long:

  Long-format weather data.

- var:

  Exposure variable (e.g. \`"tmax"\`).

- lag_max:

  Optional maximum lag. Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- df_var:

  Optional degrees of freedom (exposure). Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- df_lag:

  Optional degrees of freedom (lag). Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- fun_var:

  Optional basis (\`"ns"\`, \`"bs"\`, \`"poly"\`, \`"lin"\`). Ignored if
  \`fit\` contains \`epiexposure_spec\`.

- fun_lag:

  Optional basis (\`"ns"\`, \`"ps"\`, \`"lin"\`). Ignored if \`fit\`
  contains \`epiexposure_spec\`.

- type:

  Reduction type: \`"overall"\`, \`"lag"\`, or \`"var"\`.

- value:

  Required for \`type = "lag"\` or \`type = "var"\`.

- scale:

  Output scale: \`"link"\`, \`"response"\`, or \`"percent"\`.

- uncertainty:

  Logical. If \`TRUE\`, propagate uncertainty using simulated or
  posterior draws of the model coefficients.

- output:

  Character. \`"summary"\` returns aggregated estimates; \`"samples"\`
  returns all simulated values.

- n_samples:

  Integer. Number of samples used for uncertainty propagation.

## Value

A data.frame containing reduced effects. When \`uncertainty = TRUE\` and
\`output = "summary"\`, the result includes: - central estimate (median;
stored in \`eta\`) - \`eta_sd\`: standard deviation of simulated
values - \`low\` / \`high\`: empirical interval limits (quantiles)

## Details

When \`uncertainty = TRUE\`, summaries are computed from
simulated/posterior samples of the model coefficients. The central
estimate is obtained as the median of the simulated effects, while
uncertainty intervals are derived from empirical quantiles (default: 2.5

\*\*Important:\*\* although the output element is named \`mean\` for
backward compatibility, it represents the \*central estimate\*, computed
as the median when uncertainty is propagated.

Uncertainty is propagated using model-consistent sampling: - Bayesian
models (e.g., \`brms\`, \`INLA\`, \`bdlnm\`) use posterior draws -
Frequentist models use a normal approximation of the coefficient
distribution

The use of the median as the central estimate improves robustness under
non-normal or asymmetric effect distributions, which commonly arise in
DLNM applications.
