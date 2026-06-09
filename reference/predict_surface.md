# Predict full DLNM exposure-lag-response surface Computes the full DLNM exposure–lag–response surface using the fitted model,#' without refitting. The surface is evaluated over a grid of exposure values and lags defined by the model specification or user input.

This function supports both deterministic predictions and uncertainty
propagation. When \`uncertainty = TRUE\`, the surface is recomputed
across simulated or posterior draws of the model coefficients.

## Usage

``` r
predict_surface(
  fit,
  wx_long,
  var,
  lag_max = NULL,
  df_var = NULL,
  df_lag = NULL,
  fun_var = NULL,
  fun_lag = NULL,
  ref = list(method = "median", value = NULL),
  probs = seq(0.05, 0.95, by = 0.01),
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

  Exposure variable name.

- lag_max:

  Optional maximum lag. Ignored if fit contains \`epiexposure_spec\`.

- df_var:

  Optional degrees of freedom (exposure). Ignored if fit contains
  metadata.

- df_lag:

  Optional degrees of freedom (lag). Ignored if fit contains metadata.

- fun_var:

  Optional basis function for exposure.

- fun_lag:

  Optional basis function for lag.

- ref:

  Reference exposure definition.

- probs:

  Quantiles used for the exposure grid.

- uncertainty:

  Logical. If \`TRUE\`, quantify uncertainty.

- output:

  Character. \`"summary"\` or \`"samples"\`.

- n_samples:

  Number of simulations.

## Value

\- If \`uncertainty = FALSE\`: a \`crosspred\` object.

\- If \`uncertainty = TRUE\` and \`output = "summary"\`: a list with: -
\`mean\`: central surface (median-based) - \`lower\`: lower surface
(quantile-based) - \`upper\`: upper surface (quantile-based)

\- If \`uncertainty = TRUE\` and \`output = "samples"\`: a list of
\`crosspred\` objects (one per simulation).

## Details

If \`output = "summary"\`, the central surface is computed as the median
of the simulated surfaces, while interval limits are obtained from
empirical quantiles (default: 2.5

\*\*Important:\*\* although the returned object uses the name \`mean\`
for backward compatibility, it represents the \*central estimate\*,
computed as the median when uncertainty is propagated.

Uncertainty is propagated using model-consistent sampling: - Bayesian
models (e.g., \`brms\`, \`INLA\`, \`bdlnm\`) use posterior draws -
Frequentist models use simulation from the asymptotic coefficient
distribution

Using the median as the central estimate improves robustness to
asymmetric distributions commonly observed in DLNM surfaces.
