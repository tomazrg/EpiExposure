# Recenter DLNM effects using a new reference value

Recomputes the DLNM exposure-lag-response surface using a different
centering (reference) value, without refitting the model.

## Usage

``` r
recenter_effects(
  fit,
  wx_long,
  var,
  lag_max,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns",
  ref = list(method = "median", value = NULL),
  probs = seq(0.05, 0.95, by = 0.01),
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000
)
```

## Arguments

- fit:

  Fitted model from \`fit_epidlnm()\`.

- wx_long:

  Long-format weather data.

- var:

  Exposure variable (e.g. \`"tmax"\`).

- lag_max:

  Maximum lag.

- df_var:

  Degrees of freedom (exposure).

- df_lag:

  Degrees of freedom (lag).

- fun_var:

  Basis (\`"ns"\`, \`"bs"\`, \`"poly"\`, \`"lin"\`).

- fun_lag:

  Basis (\`"ns"\`, \`"ps"\`, \`"lin"\`).

- ref:

  New reference definition (centering value).

- probs:

  Quantiles used to define the exposure grid.

- uncertainty:

  Logical. If \`TRUE\`, quantify uncertainty using simulated or
  posterior draws of the model coefficients.

- output:

  Character. \`"summary"\` returns aggregated surfaces; \`"samples"\`
  returns all simulated surfaces.

- n_samples:

  Integer. Number of samples used for uncertainty propagation.

## Value

\- If \`uncertainty = FALSE\`: a \`crosspred\` object.

\- If \`uncertainty = TRUE\` and \`output = "summary"\`: a list with: -
\`mean\`: central estimate surface (median-based) - \`lower\`: lower
interval surface (quantile-based) - \`upper\`: upper interval surface
(quantile-based)

\- If \`uncertainty = TRUE\` and \`output = "samples"\`: a list with: -
\`mean\`: central estimate surface (median-based) - \`lower\`: lower
interval surface - \`upper\`: upper interval surface - \`samples\`: list
of \`crosspred\` objects for each simulation

## Details

When \`uncertainty = TRUE\`, the effect surface is recomputed across
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

The use of the median as the central estimate improves robustness to
asymmetry and non-normality in DLNM effect distributions, which commonly
arise from nonlinear exposure-lag-response relationships.
