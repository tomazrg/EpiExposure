# Compute lag-specific decomposition of Exposure Cumulative Impact (ECI)

Decomposes the model-weighted Exposure Cumulative Impact (ECI) of one or
more exposure histories into exact lag-specific contributions on the
linear-predictor scale.

## Usage

``` r
compute_ecilag(
  profile = NULL,
  data = NULL,
  group = NULL,
  time = "time",
  group_level = NULL,
  fit,
  vars = NULL,
  ref = list(method = "median", value = NULL),
  absolute = TRUE,
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000,
  interval_probs = c(0.025, 0.975),
  seed = NULL,
  extrapolation = c("error", "warn", "allow")
)
```

## Arguments

- profile:

  Optional exposure-history input. Supply exactly one of \`profile\` or
  \`data\`.

  For a single focal exposure, \`profile\` may be one finite numeric
  vector. For multiple focal exposures, supply a list containing one
  numeric vector per variable. A named list is recommended; if unnamed,
  its order must match \`vars\`.

  Exposure histories must be chronological from the oldest observation
  to the most recent observation. Under the EpiExposure exact-history
  contract, every supplied profile must contain exactly \`max_lag + 1\`
  observations. The first value corresponds to lag \`max_lag\` and the
  final value to lag 0. Histories are never truncated, padded, or
  silently realigned.

- data:

  Optional non-empty long-format data frame containing observed exposure
  histories. Supply exactly one of \`profile\` or \`data\`.

  \`data\` must contain \`group\`, \`time\`, and every exposure
  requested in \`vars\`. Each evaluated group must contain exactly
  \`max_lag + 1\` observations, with unique and equally spaced time
  values. All evaluated groups must use the same temporal spacing.

  When a method-based reference is requested, \`data\` must additionally
  contain every exposure included in the fitted model because the
  complete joint reference profile is constructed from those variables.

- group:

  Character scalar naming the independent-history grouping column in
  \`data\`. Required when \`data\` is supplied.

- time:

  Character scalar naming the chronological numeric time column in
  \`data\`. Default is \`"time"\`.

- group_level:

  Optional grouping value or vector of grouping values to evaluate. If
  \`NULL\`, all groups are evaluated in first-occurrence order.

- fit:

  Fitted EpiExposure model returned by the current \`fit_epidlnm()\`
  implementation. The fitted model must satisfy the current
  population-level expected-response, central-parameter, uncertainty,
  and exact-history metadata contracts.

- vars:

  Optional character vector naming fitted exposure variables to
  decompose.

  If the fitted model contains only one exposure, \`vars = NULL\` uses
  that exposure automatically. For multivariable fits, \`vars\` must be
  supplied explicitly.

  Each requested variable is decomposed separately. When one focal
  variable is evaluated, all other fitted exposures remain fixed at
  their joint reference values; contributions from different exposure
  variables are not pooled into a single lag decomposition.

- ref:

  Reference exposure specification.

  Two forms are supported:

  - Method-based: \`list(method = "median", value = NULL)\`,
    \`list(method = "percentile", value = p)\`, or \`list(method =
    "fixed", value = x)\`.

  - Exposure-specific: a named list containing exactly one finite
    reference value for every fitted exposure, for example \`list(tmean
    = 25, rain = 5, wetness = 70)\`.

  Method-based references require \`data\`. \`"median"\` uses each
  fitted exposure's median, \`"percentile"\` applies the same percentile
  probability to each fitted exposure, and \`"fixed"\` applies the same
  numeric value to every fitted exposure and should therefore be used
  cautiously when exposures have different units.

  The default is \`list(method = "median", value = NULL)\`.

- absolute:

  Logical. If \`TRUE\`, add the absolute lag contribution and, when
  uncertainty is requested, its uncertainty summary. The signed
  \`ECI_weighted\` column is always retained.

  \`ECI_percent\` is always based on the absolute magnitude of each lag
  contribution:

  \$\$ 100 \frac{\|C_l\|} {\sum_j \|C_j\|}. \$\$

  Consequently, positive and negative lag effects do not cancel when
  relative lag importance is calculated. If the total absolute
  contribution is zero, percentage contributions are undefined and
  returned as \`NA\`.

- uncertainty:

  Logical. If \`FALSE\`, lag contributions are calculated from the
  harmonized central fixed/population parameter estimate. If \`TRUE\`,
  joint parameter uncertainty is propagated draw by draw through the
  complete lag decomposition.

- output:

  Character. \`"summary"\` returns deterministic lag contributions when
  \`uncertainty = FALSE\`, or medians, standard deviations, and
  empirical intervals when \`uncertainty = TRUE\`.

  \`"samples"\` additionally returns the lag-specific and total weighted
  ECI values for every parameter draw and requires \`uncertainty =
  TRUE\`.

- n_samples:

  Positive integer number of parameter draws used when \`uncertainty =
  TRUE\`. At least two are required. Default is 1000.

- interval_probs:

  Numeric vector of length two specifying the empirical uncertainty
  interval. The default \`c(0.025, 0.975)\` gives a 95 percent interval.

- seed:

  Optional finite integer used for parameter-draw generation or
  sampling. When supplied, the caller's random-number state is restored
  when the function exits.

- extrapolation:

  Character controlling exposure values outside the range recorded by
  the fitted cross-basis: \`"error"\` (default), \`"warn"\`, or
  \`"allow"\`. The fitted basis specification is never re-estimated from
  the ECI-lag input data.

## Value

A list containing cumulative and lag-specific ECI results.

For a single exposure profile, the principal components are:

- \`ECI_raw\`:

  Descriptive sum of the focal exposure history.

- \`ECI_raw_centered\`:

  Descriptive sum of exposure deviations from the focal reference value.

- \`ECI_weighted\`:

  Centered cumulative exposure contribution on the linear-predictor
  scale.

- \`reference_value\`:

  Reference value for the focal exposure.

- \`max_lag\`:

  Common fitted maximum lag.

- \`n_exposure_values\`:

  Number of observations in the history; always \`max_lag + 1\`.

- \`by_lag\`:

  Data frame containing the exact contribution assigned to each lag.

With \`uncertainty = TRUE\`, \`ECI_weighted\` is summarized by its
median, standard deviation, and empirical lower and upper quantiles. The
\`by_lag\` table similarly contains draw-based summaries for signed
contributions, percentage contributions, and, when \`absolute = TRUE\`,
absolute contributions.

With \`output = "samples"\`, \`by_lag_samples\` and
\`ECI_weighted_samples\` contain the corresponding draw-level
quantities.

When \`data\` or multiple focal exposure histories are evaluated,
cumulative results are returned in \`eci_summary\` and lag-specific
results in \`by_lag\`, with the grouping column included when
applicable.

The result also stores attributes describing the exact ECI-lag
decomposition, common fitted maximum lag, history length, reference
values, and uncertainty setting.

## Details

Each evaluated exposure history is compared with an explicit joint
reference profile. The focal exposure is allowed to vary through time
while every other fitted exposure remains fixed at its reference value.
Lag-specific contributions are obtained by changing one chronological
exposure position at a time relative to that same reference profile and
reconstructing the corresponding fitted cross-basis through the standard
EpiExposure prediction pipeline.

The resulting lag contributions are additive: their sum is required to
equal the centered cumulative ECI for the complete focal exposure
trajectory. \`compute_ecilag()\` therefore performs an exact
decomposition of the fitted DLNM contrast. It does not calculate
numerical derivatives, marginal effects, or local sensitivity measures.

\## Centered cumulative impact

Let \\X\_{ref}\\ denote the fixed/population design generated when every
fitted exposure follows its reference profile. For focal exposure \\v\\,
let \\X_v\\ denote the design in which only that exposure follows the
evaluated complete trajectory. The model-weighted cumulative ECI is

\$\$ ECI\_{weighted,v} = (X_v-X\_{ref})\beta. \$\$

This quantity is a centered contrast on the linear-predictor scale. The
intercept and all non-focal exposure terms cancel because they are
identical in the target and reference designs.

\## Exact lag decomposition

For each chronological position corresponding to lag \\l\\, the function
creates an isolated target profile in which only that exposure value
differs from the reference trajectory. The resulting design contrast is

\$\$ D_l = X_l-X\_{ref}. \$\$

Its lag-specific contribution is

\$\$ C_l = D_l\beta. \$\$

The implementation explicitly verifies

\$\$ \sum\_{l=0}^{L} D_l = X_v-X\_{ref} \$\$

and therefore

\$\$ \sum\_{l=0}^{L} C_l = ECI\_{weighted,v}. \$\$

If this additivity condition is not satisfied within numerical
tolerance, the function stops rather than returning an approximate
decomposition.

This is fundamentally different from a derivative-based sensitivity
analysis. \`compute_ecilag()\` assigns the complete centered fitted
contrast to exact lag positions; it does not estimate
\\\partial\eta/\partial x_l\\.

\## Chronological and lag indexing

Input histories are chronological. For fitted maximum lag \\L\\, a
history of length \\L+1\\ is interpreted as

\$\$ (x_L, x\_{L-1}, \ldots, x_1, x_0), \$\$

where the first supplied observation corresponds to lag \\L\\ and the
last to lag 0. Returned \`by_lag\` tables are ordered by increasing lag.

\## Joint reference profile

Even when only one exposure is decomposed, the reference condition is
defined jointly for all fitted exposures. This keeps the decomposition
consistent with the fitted multivariable model. Non-focal exposures
remain fixed at their own reference values throughout the entire lag
window.

\## Parameter uncertainty

With \`uncertainty = TRUE\`, the same joint fixed/population parameter
draw is used for all lag contributions within a trajectory. For draw
\\s\\,

\$\$ C_l^{(s)} = D_l\beta^{(s)} \$\$

and

\$\$ ECI\_{weighted}^{(s)} = \sum_l C_l^{(s)}. \$\$

Additivity is verified for every draw. Summary statistics are calculated
only after the draw-specific contributions and percentage contributions
have been obtained.

Frequentist engines use joint draws derived from the fitted fixed-effect
covariance matrix. Bayesian engines use joint posterior or
approximate-posterior fixed/population coefficient draws through the
centralized EpiExposure prediction helpers.

No residual, observation, process, dispersion, or posterior-predictive
noise is added. Uncertainty therefore represents parameter uncertainty
in the fitted exposure-lag association.

## See also

\`compute_eci()\`, \`summarise_effects()\`, \`predict_outcomes()\`
