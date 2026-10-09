# Compute Exposure Cumulative Impact

Computes a descriptive raw cumulative exposure (\`ECI_raw\`) and, when a
fitted EpiExposure model is supplied, a model-weighted cumulative
exposure impact (\`ECI_weighted\`) for one or more exposure variables.

## Usage

``` r
compute_eci(
  profile = NULL,
  data = NULL,
  group = NULL,
  time = "time",
  group_level = NULL,
  fit = NULL,
  vars = NULL,
  ref = list(method = "median", value = NULL),
  scale = c("link", "response", "percent"),
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

  Optional exposure profile input. Supply exactly one of \`profile\` or
  \`data\`.

  When one exposure variable is evaluated, \`profile\` may be a finite
  numeric vector. For multiple variables, supply a list with one numeric
  vector per requested variable. A named list is recommended. If the
  list is unnamed, its order must match \`vars\`.

  Profiles must be chronological from oldest to most recent. When
  \`fit\` is supplied, every focal exposure profile must contain
  \*\*exactly\*\* \`max_lag + 1\` observations for that fitted exposure.
  Profiles that are shorter or longer are rejected; \`compute_eci()\`
  never selects a temporal window silently from a longer profile.

  When \`fit = NULL\`, no fitted \`max_lag\` exists. In that raw-only
  mode, \`profile\` must be one finite numeric vector and the function
  simply returns its descriptive cumulative sum.

- data:

  Optional non-empty long-format data frame containing observed
  chronological exposure profiles. Supply exactly one of \`profile\` or
  \`data\`.

  \`data\` must contain the grouping column named by \`group\`, the time
  column named by \`time\`, and every exposure requested in \`vars\`.
  Every evaluated group must contain exactly \`max_lag + 1\` equally
  spaced observations for the requested focal exposure(s). Groups with
  either fewer or more rows are rejected. If a method-based reference is
  used, \`data\` must additionally contain every exposure in the fitted
  model because those variables are needed to construct the joint
  reference profile.

- group:

  Character scalar naming the independent-profile grouping column in
  \`data\`. Required when \`data\` is supplied.

- time:

  Character scalar naming the chronological numeric time column in
  \`data\`. Default is \`"time"\`. Times must be unique and equally
  spaced within each evaluated group, and the temporal spacing must be
  common across groups.

- group_level:

  Optional grouping value or vector of values to evaluate. If \`NULL\`,
  all group levels are evaluated in first-occurrence order.

- fit:

  Optional fitted model returned by the current \`fit_epidlnm()\`. If
  \`NULL\`, only \`ECI_raw\` is calculated and model-weighted ECI,
  uncertainty, and response-scale quantities are unavailable.

- vars:

  Optional character vector naming fitted exposure variables to
  evaluate.

  If \`fit\` contains one exposure and \`vars = NULL\`, that exposure is
  used. For a multivariable fit, \`vars\` must be supplied explicitly.
  Each requested variable is evaluated separately against the same joint
  reference profile; the function does not sum impacts from different
  exposure variables.

- ref:

  Reference exposure specification used when \`fit\` is supplied.

  Two forms are supported:

  - Method-based: \`list(method = "median", value = NULL)\`,
    \`list(method = "percentile", value = p)\`, or \`list(method =
    "fixed", value = x)\`.

  - Exposure-specific: a named list containing exactly one finite
    numeric reference value for \*\*every fitted exposure\*\*, for
    example \`list(tmean = 25, rain = 5, wetness = 10)\`.

  Method-based references require \`data\` because the reference values
  are derived from the supplied exposure data. \`"median"\` uses each
  fitted exposure's own median. \`"percentile"\` uses the same
  percentile \`p\` for every fitted exposure. \`"fixed"\` applies the
  same numeric value to all fitted exposures and should be used
  cautiously when variables have different units.

  For direct \`profile\` input, use the exposure-specific named-list
  form.

  The default is \`list(method = "median", value = NULL)\`.

- scale:

  Character defining the interpretation and units of the single returned
  model-weighted impact column, \`ECI_weighted\`:

  \- \`"link"\`: the centered cumulative DLNM contrast \\\Delta\eta\\; -
  \`"response"\`: the absolute response-scale impact \\P - B\\, where
  \`P\` is the population expected response for the focal profile and
  \`B\` is the population expected response for the joint reference
  profile; - \`"percent"\`: \`100 \* (exp(Delta eta) - 1)\`, available
  only for log and logit links.

  Internally, \`compute_eci()\` calculates both the centered link-scale
  contrast \\\Delta\eta\\ and the response-scale difference \\\Delta =
  P-B\\. These internal quantities are intentionally \*\*not returned as
  separate \`eta\` and \`delta\` columns\*\* because they duplicate
  \`ECI_weighted\` on the corresponding scale: with \`scale = "link"\`,
  \`ECI_weighted = Delta eta\`; with \`scale = "response"\`,
  \`ECI_weighted = P - B\`. The internal quantities are still retained
  during computation because they are required for scale transformations
  and draw-by-draw uncertainty propagation.

  For a log link, \`"percent"\` is the percent relative change in the
  expected response. For a logit link, it is the percent change in odds,
  not the percentage-point change in probability.

- uncertainty:

  Logical. If \`FALSE\`, use the harmonized central fixed/population
  parameter estimate. If \`TRUE\`, propagate joint fixed/population
  parameter uncertainty draw by draw.

- output:

  Character. \`"summary"\` returns deterministic values when
  \`uncertainty = FALSE\`, or medians, SDs, and empirical intervals when
  \`uncertainty = TRUE\`. \`"samples"\` returns one row per parameter
  draw and requires \`uncertainty = TRUE\`.

- n_samples:

  Positive integer number of parameter draws used when \`uncertainty =
  TRUE\`. At least two are required. Default is 1000.

- interval_probs:

  Numeric vector of length two defining the empirical uncertainty
  interval. Default \`c(0.025, 0.975)\` gives a 95 percent interval.

- seed:

  Optional finite integer used for parameter-draw sampling. The caller's
  global random-number state is restored when the function exits.

- extrapolation:

  Character controlling exposure values outside the fitted cross-basis
  exposure range: \`"error"\` (default), \`"warn"\`, or \`"allow"\`. The
  fitted basis is never re-estimated from ECI input data.

## Value

A data frame.

With a fitted model, principal columns are:

- \`var\`:

  Focal exposure variable.

- \`reference_value\`:

  Reference exposure value for the focal variable.

- \`max_lag\`:

  Common fitted maximum lag.

- \`n_exposure_values\`:

  Number of exposure values used; always \`max_lag + 1\`.

- \`ECI_raw\`:

  Descriptive sum of the focal exposure values over the fitted lag
  window.

- \`ECI_raw_centered\`:

  Descriptive sum of exposure deviations from the focal reference value
  over the same lag window.

- \`baseline\`:

  Population-level expected response under the joint reference exposure
  profile, with fitted random effects excluded.

- \`predicted\`:

  Population-level expected response when the focal exposure follows the
  evaluated profile and every other fitted exposure remains at its
  reference value.

- \`ECI_weighted\`:

  The model-weighted exposure impact on the scale requested by
  \`scale\`. It equals the internally calculated \\\Delta\eta\\ on the
  link scale, the internally calculated \\P-B\\ on the response scale,
  or the permitted exponential percent transformation on the percent
  scale.

- \`scale\`:

  Requested weighted-ECI scale.

\`eta\` and \`delta\` are not returned as separate columns. Their
information is represented by \`ECI_weighted\` when \`scale = "link"\`
and \`scale = "response"\`, respectively.

With \`uncertainty = TRUE\` and \`output = "summary"\`, \`baseline\`,
\`predicted\`, and \`ECI_weighted\` are each summarized by their median
plus \`\_sd\`, \`\_lower\`, and \`\_upper\` columns. \`ECI_raw\` and
\`ECI_raw_centered\` are exposure-profile descriptors and therefore
remain deterministic.

With \`uncertainty = TRUE\` and \`output = "samples"\`, \`baseline\`,
\`predicted\`, and \`ECI_weighted\` are returned draw by draw with a
\`sample\` column. The internal link-scale and response-scale contrasts
are calculated for every draw but are not duplicated as output columns.

When \`fit = NULL\`, the result contains only \`ECI_raw\` and the number
of exposure values used.

The output attribute \`epiexposure_eci_history_contract\` records
\`"exact_max_lag_plus_one"\` whenever a fitted model is used.

## Details

Exposure profiles are interpreted chronologically, from the oldest
observation to the most recent observation. The most recent value is
associated internally with lag 0.

The model-weighted ECI is a \*\*contrast relative to an explicit joint
exposure reference profile\*\*. It is not the uncentered cross-basis
contribution \`cb(x) is not, by itself, an interpretable exposure
effect.

\## Raw ECI

For the discrete exposure profile used by the fitted lag window,

\$\$ ECI\_{raw} = \sum\_{j=0}^{L} x_j. \$\$

This is a discrete cumulative exposure sum, not a continuous-time
integral. It retains the units of the exposure multiplied by the number
of discrete lag positions. Comparisons are meaningful only for the same
exposure variable, units, lag window, and temporal resolution.

The centered descriptive counterpart is

\$\$ ECI\_{raw,centered} = \sum\_{j=0}^{L}(x_j-x\_{ref}). \$\$

Neither raw quantity incorporates the fitted exposure-response or
lag-response shape.

\## Exact exposure-profile length

When a fitted model is supplied, \`compute_eci()\` uses the EpiExposure
exact temporal-profile contract. For a fitted maximum lag \`L\`, every
evaluated focal exposure profile must contain exactly

\$\$ L + 1 \$\$

chronological observations: the first value represents lag \`L\` and the
last value represents lag 0. All fitted exposures share the same
\`max_lag\`. Longer profiles are not truncated with \`tail()\` and
shorter profiles are not padded. This prevents ambiguity about which
exposure window generated the cumulative impact and keeps ECI
calculations comparable across profiles, groups, and downstream lag
decompositions.

In raw-only mode (\`fit = NULL\`), this restriction cannot be inferred
because no fitted \`max_lag\` is available; the supplied vector is
summed as given.

\## Centered weighted ECI

A DLNM cross-basis is a basis expansion. The uncentered quantity

\$\$ cb(x)^T\beta \$\$

is a contribution to the linear predictor under the chosen basis
parameterization, but it is not an invariant exposure effect. In
particular, applying the inverse link to this contribution alone does
not produce a valid expected response because the model intercept and
the reference contribution of the other fitted exposures are absent.

EpiExposure therefore defines the model-weighted cumulative impact as a
contrast between the focal profile and an explicit joint reference
profile:

\$\$ \Delta\eta = \\cb(x)-cb(x\_{ref})\\^T\beta\_{focal}. \$\$

Every non-focal fitted exposure is held at its own reference value
throughout its complete lag window.

This centering makes \`ECI_weighted\` consistent with the DLNM contrast
logic used by \`summarise_effects()\`.

\## Why only \`ECI_weighted\` is returned

Three mathematical quantities are required internally:

\$\$ \Delta\eta = \eta\_{target} - \eta\_{reference}, \$\$

\$\$ B = g^{-1}(\eta\_{reference}), \qquad P = g^{-1}(\eta\_{target}),
\$\$

and

\$\$ \Delta = P-B. \$\$

Earlier output designs could expose \`eta\`, \`delta\`, and
\`ECI_weighted\` simultaneously. This is redundant because the selected
scale already defines which model-weighted impact is being reported. The
current output therefore uses one canonical column:

- \`scale = "link"\`: \`ECI_weighted\` is \\\Delta\eta\\;

- \`scale = "response"\`: \`ECI_weighted\` is \\\Delta = P-B\\;

- \`scale = "percent"\`: \`ECI_weighted\` is the permitted percent
  transformation of \\\Delta\eta\\.

\`baseline\` and \`predicted\` remain in the output because they are
distinct expected-response quantities that provide the response-scale
context for the contrast. The intermediate \`eta\` and \`delta\`
quantities continue to be computed internally, including for every
uncertainty draw, but are not duplicated as public output columns.

\## Response-scale ECI

Let \\X\_{ref}\\ denote the complete fixed/population design under the
joint reference exposure profile. The population expected-response
baseline is

\$\$ B = g^{-1}(X\_{ref}\beta). \$\$

Let \\X\_{target}\\ replace only the focal exposure profile by the
evaluated profile. Then

\$\$ P = g^{-1}(X\_{target}\beta) \$\$

and

\$\$ \Delta = P-B. \$\$

With \`scale = "response"\`, \`ECI_weighted = Delta\`.

\## Percent scale

With a log link,

\$\$ 100\\\exp(\Delta\eta)-1\\ \$\$

is the percent relative change in the population expected response.

With a logit link, the same mathematical transformation is the percent
change in odds. It is not the percent change or percentage-point change
in predicted probability.

Other links do not receive this exponential interpretation, so \`scale =
"percent"\` is rejected explicitly.

\## Population-level prediction contract

Weighted ECI follows the same EpiExposure v1 prediction target as
\`predict_outcomes()\` and \`summarise_effects()\`:

\$\$ E(Y\mid X,\theta) \$\$

using only the fixed/population component. Fitted group-specific random
effects are set to zero/excluded. No residual, observation, process,
dispersion, or posterior-predictive noise is added.

With \`uncertainty = FALSE\`, frequentist fits use fitted fixed
coefficients and Bayesian fits use posterior-mean fixed/population
coefficients.

\## Uncertainty

With \`uncertainty = TRUE\`, one coherent fixed/population parameter
draw is applied jointly to the reference design and every target ECI
design. For draw \\s\\,

\$\$ \Delta\eta_i^{(s)} = (X_i-X\_{ref})\beta^{(s)}, \$\$

\$\$ B^{(s)} = g^{-1}(X\_{ref}\beta^{(s)}), \$\$

\$\$ P_i^{(s)} = g^{-1}(X_i\beta^{(s)}), \$\$

and

\$\$ \Delta_i^{(s)} = P_i^{(s)}-B^{(s)}. \$\$

The requested \`ECI_weighted\` transformation is then applied \*\*within
each draw\*\* before medians, SDs, and empirical quantiles are
calculated. The same parameter draw determines the baseline and every
target row, preserving covariance among variables, groups, baseline,
target predictions, and the weighted ECI.

Frequentist engines use the joint asymptotic fixed-effect covariance
matrix. Bayesian engines use joint posterior or approximate-posterior
fixed-effect draws through the centralized EpiExposure prediction
helpers. The function does not independently reconstruct engine-specific
uncertainty and does not add residual or posterior-predictive outcome
noise.

\## Method-based reference values

Method-based references are resolved once from the complete supplied
\`data\`, not independently within each group. This ensures every
evaluated group is compared with the same joint exposure reference
condition.
