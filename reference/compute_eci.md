# Compute Exposure Cumulative Impact (ECI - the exposure profile with the fitted model coefficients)

This function supports both deterministic estimation and uncertainty
propagation. When \`uncertainty = TRUE\`, the weighted ECI is recomputed
across simulated or posterior draws of the model coefficients.

## Usage

``` r
compute_eci(
  profile,
  fit = NULL,
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000
)
```

## Arguments

- profile:

  Numeric vector of exposure values representing a lag profile. Its
  length must be equal to \`lag_max + 1\` for the exposure variable
  stored in the fitted model.

- fit:

  Fitted model from \`fit_epidlnm()\`. Required to compute
  \`ECI_weighted\`. If \`NULL\`, only \`ECI_raw\` is returned.

- uncertainty:

  Logical. If \`TRUE\`, quantify uncertainty.

- output:

  Character. \`"summary"\` or \`"samples"\`.

- n_samples:

  Integer. Number of samples used for uncertainty quantification.

## Value

A data.frame.

\- If \`fit = NULL\`, returns: - \`ECI_raw\` - \`ECI_weighted = NA\`

\- If \`uncertainty = FALSE\`, returns: - \`ECI_raw\` - \`ECI_weighted\`

\- If \`uncertainty = TRUE\` and \`output = "summary"\`, returns: -
\`ECI_raw\` - \`ECI_weighted\` (median-based central estimate) -
\`sd\` - \`lower\` - \`upper\`

\- If \`uncertainty = TRUE\` and \`output = "samples"\`, returns: -
\`sample\` - \`ECI_raw\` - \`ECI_weighted\`

## Details

If \`output = "summary"\`, the central estimate is computed as the
median of the simulated ECI values, while interval limits are obtained
from empirical quantiles (default: 2.5

\*\*Important:\*\* when \`uncertainty = TRUE\` and \`output =
"summary"\`, \`ECI_weighted\` represents the \*central estimate\*,
computed as the median of the simulated distribution.

Uncertainty is propagated using model-consistent sampling: - Bayesian
models (e.g., \`brms\`, \`INLA\`, \`bdlnm\`) use posterior draws -
Frequentist models use simulation from the asymptotic coefficient
distribution

The use of the median as the central estimate improves robustness under
asymmetric or non-normal simulated ECI distributions.
