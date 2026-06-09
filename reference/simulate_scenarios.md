# Simulate epidemiological DLNM scenarios

Generates predicted outcomes for one or more user-defined
epidemiological scenarios, using exposure profiles assembled across lag
periods.

## Usage

``` r
simulate_scenarios(
  fit,
  scenarios,
  wx_long = NULL,
  periods = NULL,
  lag_max = NULL,
  df_var = NULL,
  df_lag = NULL,
  fun_var = NULL,
  fun_lag = NULL,
  ref_vals = NULL,
  pop_level = TRUE,
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000
)
```

## Arguments

- fit:

  Fitted model returned by \`fit_epidlnm()\`.

- scenarios:

  Output from \`simulate_range()\` or a named list of scenarios.
  Recommended structured object: \`list(scenarios = ..., info = ...,
  periods = ...)\`.

- wx_long:

  Optional long-format weather data. Used only as fallback to derive
  default reference values when \`ref_vals\` are not supplied and the
  fitted model does not store centering information.

- periods:

  Optional period table. If \`NULL\`, the function will try to read it
  from \`scenarios\$periods\`.

- lag_max:

  Optional maximum lag. Used only if \`fit\` does not store lag_max.

- df_var:

  Optional degrees of freedom (exposure). Ignored if \`fit\` contains
  \`epiexposure_spec\`. Kept only for fallback compatibility.

- df_lag:

  Optional degrees of freedom (lag). Ignored if \`fit\` contains
  \`epiexposure_spec\`. Kept only for fallback compatibility.

- fun_var:

  Optional basis function for exposure. Ignored if \`fit\` contains
  \`epiexposure_spec\`. Kept only for fallback compatibility.

- fun_lag:

  Optional basis function for lag. Ignored if \`fit\` contains
  \`epiexposure_spec\`. Kept only for fallback compatibility.

- ref_vals:

  Optional named list of reference values for each variable. If
  \`NULL\`, the function tries in order: 1. median from \`wx_long\` 2.
  centering value stored in \`fit\` spec (\`argvar\$cen\`)

- pop_level:

  Logical. If \`TRUE\`, predictions exclude random effects where
  supported.

- uncertainty:

  Logical. If \`TRUE\`, quantify uncertainty.

- output:

  Character. \`"summary"\` or \`"samples"\`.

- n_samples:

  Integer. Number of samples used for uncertainty quantification.

## Value

A data.frame.

If \`uncertainty = FALSE\`, returns one row per scenario point with
column: - \`prediction\`

If \`uncertainty = TRUE\` and \`output = "summary"\`, returns one row
per scenario point with columns: - \`prediction\` (median-based central
estimate) - \`sd\` - \`lower\` - \`upper\`

If \`uncertainty = TRUE\` and \`output = "samples"\`, returns one row
per sample with columns: - \`sample\` - \`prediction\`

Additional metadata columns from \`scenarios\$info\` are preserved.

## Details

Scenarios can be supplied either as: - a structured object returned by
\`simulate_range()\`, or - a named list of scenario definitions.

When \`uncertainty = TRUE\`, uncertainty is propagated through
\`predict_outcome()\`. If \`output = "summary"\`, the central estimate
is computed as the median of simulated predictions, and interval limits
are obtained from empirical quantiles (default: 2.5

\*\*Important:\*\* in the summary output, the column \`prediction\`
represents the central estimate, computed as the median when uncertainty
is propagated.

Uncertainty is propagated through \`predict_outcome()\`, which uses
model-consistent sampling: - Bayesian models use posterior draws -
Frequentist models use simulation from the asymptotic coefficient
distribution

For summary outputs under uncertainty, the median is used instead of the
mean to provide a more robust central estimate under asymmetric
predictive distributions.
