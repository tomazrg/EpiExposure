# Summarise DLNM exposure-lag effects on link and response scales

Summarises lag-specific or period-specific distributed lag nonlinear
model (DLNM) associations for models fitted with \`fit_epidlnm()\`. The
function preserves the centering logic of \`dlnm::crosspred()\`: each
exposure-lag association is expressed relative to a reference exposure
value (\`cen\`) for the focal exposure. EpiExposure then extends that
DLNM contrast by anchoring it to a joint reference exposure profile,
allowing the same association to be reported as \`baseline\`,
\`predicted\`, and \`delta\` on the response scale.

## Usage

``` r
summarise_effects(
  fit,
  data,
  vars = NULL,
  scale = c("lag", "period"),
  lag_periods = NULL,
  probs = seq(0.05, 0.95, by = 0.01),
  at = NULL,
  ref = list(method = "median", value = NULL),
  effect_measure = c("linear", "exponentiated", "percent"),
  incremental = FALSE,
  uncertainty = FALSE,
  output = c("summary", "samples"),
  n_samples = 1000,
  interval_probs = c(0.025, 0.975),
  diagnostics = FALSE,
  seed = NULL,
  extrapolation = c("error", "warn", "allow")
)
```

## Arguments

- fit:

  Fitted model returned by the current \`fit_epidlnm()\`.

- data:

  Long-format exposure data containing \`epi_id\`, \`time\`, and every
  exposure variable used in the fitted model. \`data\` is used to define
  default exposure grids and reference values and to validate temporal
  coverage. Under the EpiExposure exact-profile contract, every
  \`epi_id\` must contain exactly the common fitted \`max_lag + 1\` time
  points, and all fitted exposures therefore share the same
  temporal-profile length. Longer profiles are not truncated and shorter
  profiles are not padded. \`data\` is \*\*not\*\* used to re-estimate
  spline knots, boundary knots, lag bases, or any other fitted
  cross-basis component.

- vars:

  Exposure-variable name or names to summarise, or \`NULL\` to summarise
  all fitted exposures.

- scale:

  Character. \`"lag"\` returns lag-specific associations. \`"period"\`
  aggregates lag-specific linear-predictor contrasts over intervals.

- lag_periods:

  Optional data frame containing \`period\`, \`lag_start\`, and
  \`lag_end\`. Required when \`scale = "period"\` and \`incremental =
  FALSE\`. Both \`lag_start\` and \`lag_end\` are interpreted on the
  retrospective lag scale, where lag 0 is the most recent exposure and
  increasing lag values correspond to older exposure observations.
  Period bounds are inclusive. Periods may be non-overlapping or
  overlapping; each period is interpreted independently.

- probs:

  Unique probabilities used only to construct the default exposure
  evaluation grid when \`at\` is not supplied for a variable. Defaults
  to \`seq(0.05, 0.95, by = 0.01)\`.

- at:

  Optional exposure values at which associations are evaluated. \`NULL\`
  uses the quantiles specified by \`probs\`. A numeric vector can be
  used when exactly one exposure is requested. For multiple exposures,
  supply a named list, for example \`list(tmean = seq(20, 35, by =
  0.25), rain = seq(0, 40, by = 1))\`. A named list may be partial;
  omitted requested variables use \`probs\`.

- ref:

  Reference exposure specification.

  Two forms are supported:

  - Method-based: \`list(method = "median", value = NULL)\`,
    \`list(method = "percentile", value = p)\`, or \`list(method =
    "fixed", value = x)\`.

  - Exposure-specific: a named list with exactly one finite reference
    value for every fitted exposure, for example \`list(tmean = 25, rain
    = 5, wetness = 10)\`.

  Method-based references are applied consistently to \*\*every\*\*
  fitted exposure. Thus \`method = "median"\` uses each exposure's own
  median; \`method = "percentile"\` uses the same percentile for each
  exposure; and \`method = "fixed"\` applies the same numeric value to
  every exposure. Because exposures often have different units, the
  exposure-specific form is recommended for multivariable models when
  fixed reference values are desired.

  The focal exposure's reference value is passed to the DLNM centering
  operation (\`cen\`). The complete set of reference values
  simultaneously defines the joint reference profile used for
  \`baseline\`.

- effect_measure:

  Character. \`"linear"\` returns the DLNM contrast on the
  linear-predictor scale. \`"exponentiated"\` returns \`exp(eta)\` and
  \`"percent"\` returns \`100 \* (exp(eta) - 1)\`.

  Exponentiation is allowed only for fitted \`"log"\` and \`"logit"\`
  links, matching the convention used by \`dlnm::crosspred()\`. With a
  log link, \`exp(eta)\` is a response-scale ratio and \`"percent"\` is
  the corresponding percent relative change in the expected response.
  With a logit link, \`exp(eta)\` is an odds ratio and \`"percent"\` is
  the percent change in odds; neither is a percentage-point change in
  probability. For links such as \`"identity"\`, \`"probit"\`,
  \`"cloglog"\`, or \`"inverse"\`, use \`effect_measure = "linear"\`.

- incremental:

  Logical. Used only with \`scale = "period"\`. If \`TRUE\`, return
  cumulative effects from lag 0 through each successive lag. If
  \`FALSE\`, aggregate over \`lag_periods\`.

- uncertainty:

  Logical. If \`FALSE\`, use the harmonized central fixed/population
  parameter estimate. If \`TRUE\`, propagate joint fixed/population
  parameter uncertainty draw by draw.

- output:

  Character. \`"summary"\` returns deterministic values when
  \`uncertainty = FALSE\`, or median, SD, and empirical intervals when
  \`uncertainty = TRUE\`. \`"samples"\` returns one row per parameter
  draw and exposure-lag/period combination and therefore requires
  \`uncertainty = TRUE\`.

- n_samples:

  Positive integer number of parameter draws used when \`uncertainty =
  TRUE\`. At least two are required.

- interval_probs:

  Numeric vector of length two defining the empirical uncertainty
  interval. The default \`c(0.025, 0.975)\` gives a 95 percent interval.

- diagnostics:

  Logical. If \`TRUE\`, additionally return lag-ranking and
  lag-contribution diagnostics. Diagnostics are available only for
  \`scale = "lag"\` and are always based on the additive
  linear-predictor contrast \`eta\`, regardless of \`effect_measure\`.
  This avoids treating the neutral value 1 of an exponentiated effect as
  if it were an additive contribution.

- seed:

  Optional finite integer used for parameter sampling. The caller's
  global random-number state is restored when the function exits.

- extrapolation:

  Character controlling values outside the exposure range stored with
  the fitted cross-basis: \`"error"\` (default), \`"warn"\`, or
  \`"allow"\`. This does not re-estimate the basis; it only controls
  whether extrapolation is rejected, warned about, or allowed.

## Value

A data frame, or when \`diagnostics = TRUE\`, an object of class
\`"epiexposure_effects"\` containing \`effects\` and \`diagnostics\`.

The principal columns have the following meanings:

- \`eta\`:

  DLNM association contrast on the fitted link/linear-predictor scale,
  centered at the focal exposure reference. For \`scale = "lag"\`, this
  is the contribution of the specified exposure value at that lag
  relative to the same lag at the reference exposure. For period
  summaries, lag-specific \`eta\` contrasts are summed on the additive
  linear-predictor scale.

- \`effect\`:

  \`eta\` itself for \`effect_measure = "linear"\`, or its permitted
  log/logit transformation for \`"exponentiated"\` or \`"percent"\`.

- \`baseline\`:

  Population-level expected outcome under the \*\*joint reference
  exposure profile\*\*: every fitted exposure is held at its own
  reference value at every fitted lag and fitted random effects are set
  to zero. This is an EpiExposure response-scale reference prediction;
  it is not the DLNM contrast itself and should not be confused with the
  neutral DLNM values \`eta = 0\` or \`exp(eta) = 1\`.

- \`predicted\`:

  Population-level expected outcome after applying the focal DLNM
  contrast to that joint baseline. For a lag-specific row, only the
  focal exposure at that lag is changed from its reference to \`value\`;
  all other lags and all other fitted exposures remain at their
  reference values. For a period row, the focal exposure is changed to
  \`value\` throughout the specified lag interval while all remaining
  exposure-lag positions remain at reference.

- \`delta\`:

  Absolute response-scale difference \`predicted - baseline\`. For
  Gaussian/Gamma-type mean responses this is a difference in expected
  means; for Poisson/negative-binomial outcomes it is a difference in
  expected counts; and for Beta/Binomial outcomes it is a difference in
  expected proportions/probabilities. Multiplying \`delta\` by 100 for
  Beta or Binomial models gives percentage-point differences.

## Details

Exposure profiles are supplied chronologically: the oldest exposure is
first and the most recent exposure is last. Internally, EpiExposure
converts these profiles to the retrospective lag scale used by DLNMs,
where lag 0 represents the most recent exposure observation and
increasing lag values represent progressively older exposure
observations. All lag-specific and period-specific summaries reported by
\`summarise_effects()\` are therefore indexed on this retrospective lag
scale.

\## DLNM contrast and the EpiExposure response-scale extension

\`eta\` is constructed using the centering mechanism of
\`dlnm::crosspred()\`. The fitted cross-basis is never redefined from
\`data\`: EpiExposure uses the stored \`argvar\`, \`arglag\`, maximum
lag, and, when available, the original stored \`crossbasis\` object.
Consequently, \`eta\` retains the standard DLNM interpretation as an
association relative to \`cen\`.

EpiExposure additionally defines a joint response-scale reference
prediction

\$\$B = g^{-1}(X\_{ref}\beta),\$\$

where every fitted exposure profile is held at its selected reference
value. A lag- or period-specific DLNM contrast \\\Delta\eta\_{x,l}\\ is
then translated to the response scale as

\$\$P\_{x,l} = g^{-1}\\g(B) + \Delta\eta\_{x,l}\\,\$\$

and

\$\$D\_{x,l} = P\_{x,l} - B.\$\$

Thus \`baseline\`, \`predicted\`, and \`delta\` do not replace the DLNM
contrast; they provide an additional outcome-scale interpretation
anchored to the same reference condition.

\## Uncertainty

With \`uncertainty = TRUE\`, one \*\*joint\*\* fixed/population
coefficient draw is used consistently across the complete requested
effect surface. For draw \\s\\,

\$\$B^{(s)} = g^{-1}(X\_{ref}\beta^{(s)}),\$\$

\$\$P^{(s)}\_{x,l} = g^{-1}\\X\_{ref}\beta^{(s)} +
\Delta\eta^{(s)}\_{x,l}\\,\$\$

and

\$\$D^{(s)}\_{x,l} = P^{(s)}\_{x,l} - B^{(s)}.\$\$

The same draw therefore determines \`baseline\`, \`eta\`, \`effect\`,
\`predicted\`, and \`delta\`, preserving their covariance. \`baseline\`
changes across parameter draws because the fitted parameters are
uncertain, but for a given draw it remains the same across all
exposure-lag/period rows that use the same joint reference profile. In
summary output, the same baseline uncertainty summary is therefore
repeated across those rows.

Parameter uncertainty includes the joint uncertainty of the fitted
fixed/population coefficients. Fitted group-specific random effects are
excluded because EpiExposure v1 reports population-level effects.
Residual, observation, process, dispersion, and posterior-predictive
noise are not added. The uncertainty interval therefore describes
uncertainty in the expected response and in the exposure-lag
association, not the dispersion of a future individual observation.

For uncertainty summaries, all response-scale transformations are
performed draw by draw before medians, SDs, and empirical quantiles are
calculated.

\## Exact common lag/profile contract

\`summarise_effects()\` requires the fitted model to use one common
\`max_lag\` across all fitted exposure variables. If that common maximum
lag is \`L\`, every epidemic profile supplied in \`data\` must contain
exactly

\$\$ L + 1 \$\$

equally spaced observations. Profiles with fewer or more observations
are rejected explicitly. Because all exposure variables are columns of
the same validated long-format profile and missing/non-finite exposure
values are not permitted, all fitted exposures necessarily use the same
time support.

\## Period effects

Period definitions always use the retrospective lag scale employed by
EpiExposure and DLNMs. Thus, lag 0 corresponds to the most recent
exposure observation and increasing lag values correspond to
progressively older exposure observations. For example, if \`max_lag =
85\`, a period defined as \`0–14\` represents the 15 most recent
exposure observations, whereas a period defined as \`71–85\` represents
the oldest portion of the fitted exposure profile.

DLNM contributions are additive on the linear-predictor scale. Period
summaries therefore sum lag-specific \`eta\` contrasts first and only
then transform the resulting contrast, if requested. With \`incremental
= TRUE\`, periods are \`0-0\`, \`0-1\`, ..., \`0-L\`. Otherwise each row
of \`lag_periods\` defines an independent lag interval.
