# Summarise DLNM effects

Summarises DLNM effects on either the daily or accumulated scale, using
exposure grids derived from the observed data and the DLNM specification
stored in the fitted model.

## Usage

``` r
summarise_effects(
  fit,
  wx_long,
  var = NULL,
  lag_max = NULL,
  df_var = NULL,
  df_lag = NULL,
  fun_var = NULL,
  fun_lag = NULL,
  scale = c("daily", "accumulated"),
  lag_periods = NULL,
  probs = seq(0.05, 0.95, by = 0.01),
  ref = list(method = "median", value = NULL),
  effect_measure = c("percent", "ratio", "linear"),
  incremental = FALSE,
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000
)
```

## Arguments

- fit:

  Fitted model.

- wx_long:

  Long-format weather data.

- var:

  Exposure variable name (character scalar), character vector of
  variables, or \`NULL\`. If \`NULL\`, all variables stored in fit
  metadata are used.

- lag_max:

  Optional maximum lag. Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- df_var:

  Optional degrees of freedom for exposure. Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- df_lag:

  Optional degrees of freedom for lag. Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- fun_var:

  Optional exposure basis function. Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- fun_lag:

  Optional lag basis function. Ignored if \`fit\` contains
  \`epiexposure_spec\`.

- scale:

  Character. \`"daily"\` or \`"accumulated"\`.

- lag_periods:

  Optional lag-period table used when \`scale = "accumulated"\` and
  \`incremental = FALSE\`. The table must contain columns \`period\`,
  \`lag_start\`, and \`lag_end\`.

- probs:

  Quantiles used to define the exposure grid.

- ref:

  Reference definition list.

- effect_measure:

  Character. \`"percent"\`, \`"ratio"\`, or \`"linear"\`.

- incremental:

  Logical. If \`TRUE\` and \`scale = "accumulated"\`, returns cumulative
  effects from lag 0 up to each lag.

- uncertainty:

  Logical. If \`TRUE\`, quantify uncertainty using simulated or
  posterior coefficient draws.

- output:

  Character. \`"summary"\` or \`"samples"\`.

- n_samples:

  Integer. Number of samples used for uncertainty.

## Value

A data.frame.

If \`uncertainty = FALSE\`, returns deterministic summaries of DLNM
effects.

If \`uncertainty = TRUE\` and \`output = "summary"\`, returns one row
per grid value / lag / period combination (depending on \`scale\`) with
median-based central estimates and empirical interval limits.

If \`uncertainty = TRUE\` and \`output = "samples"\`, returns one row
per simulated sample.

## Details

This function supports both deterministic summaries and uncertainty
propagation. When \`uncertainty = TRUE\`, summaries are computed from
simulated or posterior draws of the model coefficients. If \`output =
"summary"\`, the central estimate is computed as the median of the
simulated effects, while interval limits are obtained from empirical
quantiles (default: 2.5

\*\*Important:\*\* when \`uncertainty = TRUE\` and \`output =
"summary"\`, columns such as \`eta\`, \`effect\`, \`delta\`,
\`delta_pp\`, \`baseline\`, and \`predicted\` represent the \*central
estimate\*, computed as the median of the simulated distribution.

Uncertainty is propagated using model-consistent sampling: - Bayesian
models (e.g., \`brms\`, \`INLA\`, \`bdlnm\`) use posterior draws -
Frequentist models use simulation from the asymptotic coefficient
distribution

For summary outputs under uncertainty, the median is used instead of the
mean to provide a more robust central estimate under asymmetric effect
distributions, which are common in nonlinear DLNM settings.
