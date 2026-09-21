# Compute sensitivity of exposure-response curves

Estimates the local derivative of an already calculated
exposure-response curve. The function can work directly with supplied
curve values or first smooth the curve with a GAM and differentiate the
fitted smooth.

## Usage

``` r
epi_sensitivity(
  data,
  x = "value",
  y = "prediction",
  scenario_var = "scenario",
  k = 10,
  method = c("analytical", "finite"),
  smooth_basis = c("cs", "cr", "tp", "ts"),
  elasticity = TRUE,
  critical = TRUE,
  eps = 1e-06
)
```

## Arguments

- data:

  Non-empty data frame containing the curve to differentiate. At
  minimum, it must contain the columns selected by \`x\` and \`y\`.

  Within each group defined by \`scenario_var\`, \`x\` must contain
  unique finite numeric values. Duplicate \`x\` values are rejected
  rather than aggregated silently because a derivative requires one
  unambiguous curve value at each evaluation position.

- x:

  Character scalar naming the predictor/evaluation coordinate. Default
  is \`"value"\`.

- y:

  Character scalar naming the curve value to differentiate. Default is
  \`"prediction"\`.

- scenario_var:

  \`NULL\` or a character vector naming one or more grouping columns.
  The derivative is computed independently within every unique
  combination of these columns.

  The historical default is \`"scenario"\`. Set \`scenario_var = NULL\`
  to differentiate the entire supplied data set as one curve. If several
  exposure variables, scenarios, or parameter draws are present, include
  all columns needed to identify one unique curve, for example
  \`scenario_var = c("var", "scenario")\` or \`scenario_var =
  c("scenario", "sample")\`.

- k:

  Positive integer GAM basis dimension requested when \`method =
  "analytical"\`. Default is \`10\`.

  The effective value is reduced, when necessary, to the number of
  distinct \`x\` values in the current curve. It is never increased
  silently. A non-constant curve analyzed with \`method = "analytical"\`
  requires at least three distinct \`x\` values and \`k \>= 3\`.

- method:

  Character. One of:

  \- \`"analytical"\`: historical EpiExposure name for a GAM-based
  smooth derivative. A univariate GAM is fitted with \`mgcv::gam(...,
  method = "REML")\`, using the auxiliary smooth basis selected by
  \`smooth_basis\`, then its linear-predictor matrix is differentiated
  numerically using a small within-range perturbation; - \`"finite"\`:
  differentiates the supplied curve directly using local
  finite-difference formulas. With at least three points, second-order
  three-point formulas are used and support irregular \`x\` spacing.

  The label \`"analytical"\` is retained for backward compatibility, but
  this method is not a symbolic analytical derivative. It is a
  derivative of the fitted GAM smooth obtained by finite differencing
  its \`lpmatrix\`.

- smooth_basis:

  Character. Auxiliary \`mgcv\` smooth basis used only when \`method =
  "analytical"\`. One of:

  \- \`"cs"\` (default): shrinkage cubic regression spline; - \`"cr"\`:
  cubic regression spline without the additional shrinkage
  modification; - \`"tp"\`: thin plate regression spline; - \`"ts"\`:
  shrinkage thin plate regression spline.

  This choice belongs only to the \*\*post-processing GAM fitted inside
  \`epi_sensitivity()\`\*\*. It does \*\*not\*\* redefine, replace,
  inherit, or need to match the exposure-response or lag-response basis
  used by \`define_exposures()\`, \`fit_epidlnm()\`, or
  \`find_bestfit()\`. The DLNM has already been fitted before this
  function is called.

  \`"cs"\` remains the default because it provides a low-rank univariate
  cubic smooth with shrinkage, which is a useful conservative default
  when differentiating an already estimated curve. The alternatives are
  exposed to support transparent sensitivity/robustness analyses of this
  secondary smoothing step.

  \`smooth_basis\` is validated even when \`method = "finite"\` for API
  consistency, but it has no computational effect in that method and the
  output metadata records it as unused.

- elasticity:

  Logical scalar. If \`TRUE\`, compute local elasticity

  \$\$ E(x) = \frac{dy}{dx}\frac{x}{y}. \$\$

  For \`method = "analytical"\`, the denominator is the GAM-fitted curve
  value, so the derivative and elasticity refer to the same fitted
  function. For \`method = "finite"\`, the denominator is the supplied
  \`y\`.

  Elasticity is returned as \`NA\` when its denominator is numerically
  indistinguishable from zero.

- critical:

  Logical scalar. If \`TRUE\`, detect derivative sign changes. A sign
  change from positive to negative is classified as a local maximum; a
  change from negative to positive is classified as a local minimum.

  Because the true zero crossing can lie between evaluated \`x\` values,
  the nearest supplied/evaluated row is flagged and the interpolated
  crossing is returned in \`critical_x\`.

- eps:

  Positive finite numeric scalar smaller than 1. Default is \`1e-6\`.

  It is used as a relative numerical tolerance. For the GAM derivative,
  the perturbation is proportional to the observed \`x\` span and is
  constrained to remain inside that span. For elasticity and
  critical-point detection, scale-aware tolerances are derived from the
  corresponding curve/derivative magnitudes.

## Value

A data frame containing all original input columns, ordered within group
by increasing \`x\`, plus:

\- \`sensitivity\`: estimated local derivative \`dy/dx\`; -
\`sensitivity_curve_value\`: value of the curve whose derivative is
being reported. This equals the GAM-fitted value for \`method =
"analytical"\` and the supplied \`y\` for \`method = "finite"\`; -
\`elasticity\`, when \`elasticity = TRUE\`; - \`critical\`,
\`critical_type\`, and \`critical_x\`, when \`critical = TRUE\`.

\`critical_type\` is \`"local_maximum"\` or \`"local_minimum"\` for
detected derivative sign changes and \`NA\` otherwise. \`critical_x\` is
a linear interpolation of the derivative zero crossing between the
bracketing evaluation points.

Output attributes document the derivative method, grouping variables,
numerical contract, requested/effective GAM basis dimensions, and the
auxiliary \`smooth_basis\` used by the GAM derivative.

## Details

\`epi_sensitivity()\` is a \*\*post-processing\*\* function. It does not
fit an EpiExposure DLNM, reconstruct exposure histories, or redefine the
fitted lag window. Consequently, rows supplied to this function are
curve-evaluation points, not necessarily the \`max_lag + 1\`
observations that constituted an original exposure history.

\## GAM-based derivative

For \`method = "analytical"\`, EpiExposure fits the descriptive curve

\$\$ y = s(x) + \epsilon \$\$

with REML smoothing-parameter estimation and the auxiliary \`mgcv\`
basis selected through \`smooth_basis\`.

The supported bases are deliberately restricted to four standard
one-dimensional choices:

\- \`"cs"\`: shrinkage cubic regression spline; - \`"cr"\`: cubic
regression spline; - \`"tp"\`: thin plate regression spline; - \`"ts"\`:
shrinkage thin plate regression spline.

This auxiliary spline is \*\*not the DLNM spline\*\*. It is fitted only
after an exposure-response curve has already been produced. For example,
a DLNM may have been fitted with a natural spline exposure basis and the
resulting curve may still be differentiated here with \`smooth_basis =
"cs"\` or \`"tp"\`. There is no requirement that the two bases match
because they serve different statistical roles.

The DLNM basis determines the epidemiological model itself. In contrast,
\`smooth_basis\` controls only how the already evaluated curve is
optionally smoothed before numerical differentiation. Users interested
in robustness can therefore compare \`"cs"\`, \`"cr"\`, \`"tp"\`, and
\`"ts"\` without changing the previously fitted DLNM.

Let \\X_p(x)\\ denote the GAM linear-predictor matrix and \\\hat\beta\\
its fitted coefficient vector. The fitted curve is

\$\$ \hat f(x) = X_p(x)\hat\beta. \$\$

The derivative mapping is approximated from two nearby prediction
matrices:

\$\$ X'\_p(x) \approx \frac{X_p(x\_+) - X_p(x\_-)}{x\_+ - x\_-}, \$\$

and

\$\$ \hat f'(x) = X'\_p(x)\hat\beta. \$\$

Interior points use a centered perturbation. Boundary points use a
one-sided perturbation that remains inside the observed \`x\` range; the
function does not intentionally extrapolate beyond the supplied curve.

The external method name \`"analytical"\` is retained because it was
part of the previous EpiExposure API. Methodologically, however, this is
a GAM-smoothed numerical derivative, not symbolic differentiation.

\## Why \`"cs"\` remains the default

Derivatives can amplify small-scale wiggles in an estimated curve. A
shrinkage cubic regression spline is therefore retained as the default
auxiliary smoother because it offers a simple one-dimensional cubic
basis while allowing stronger penalization of weak structure. This is a
post-processing default, not a statement that \`"cs"\` is universally
optimal.

\`smooth_basis = "cr"\` is useful when the analyst wants the
corresponding cubic regression spline without the extra shrinkage
modification. \`"tp"\` provides the general thin plate regression spline
commonly used by \`mgcv\`, while \`"ts"\` adds shrinkage to that family.

Because the estimated derivative can depend on the secondary smoothing
choice, reporting or checking sensitivity across these bases can be
useful when derivative-based scientific conclusions depend on fine
features of the curve.

\## Direct finite differences

With \`method = "finite"\`, no smoother is fitted. For three or more
points, the derivative at an interior point is obtained from the
quadratic interpolant through the previous, current, and next
observations. This gives the standard centered three-point derivative
for equally spaced data and its corresponding unequal-spacing form for
irregular \`x\`.

Endpoints use the corresponding one-sided three-point formula. With
exactly two points, the secant slope is the only estimable derivative
and is assigned to both positions.

\`x\` values must be strictly unique after grouping. The function does
not average duplicate positions, discard rows, or replace a zero spacing
by \`eps\`.

\## Elasticity

Elasticity is dimensionless:

\$\$ E(x) = f'(x)x/f(x). \$\$

For the GAM method, both \\f'(x)\\ and \\f(x)\\ come from the same
fitted smooth. For the direct finite method, \\f(x)\\ is the supplied
\`y\`.

Elasticity can be numerically unstable near zero. EpiExposure therefore
returns \`NA\` when the relevant curve value is within a scale-aware
tolerance of zero. No arbitrary constant is added to the denominator.

\## Critical points

\`critical = TRUE\` searches for changes in the sign of the estimated
derivative. Near-zero derivative values are treated as a bridge between
the nearest non-zero derivative signs. Only genuine positive-to-negative
or negative-to-positive transitions are marked.

Therefore:

\- positive -\> negative: local maximum; - negative -\> positive: local
minimum.

Flat regions without a sign reversal are not automatically labelled as
extrema. Boundary extrema are also not inferred from a one-sided
derivative alone.

\## Relationship to the EpiExposure exact-history contract

\`epi_sensitivity()\` does not consume raw fitted exposure histories. It
works on an already evaluated curve, whose number of rows may
legitimately differ from the original history length. Consequently, it
does \*\*not\*\* require

\$\$ nrow(data) = max\\lag + 1. \$\$

If current EpiExposure temporal metadata are present on \`data\`,
however, the function validates and propagates them. Specifically,

\$\$ history\\length = max\\lag + 1 \$\$

and the stored history contract must be
\`"all_fitted_exposures_same_exact_max_lag_plus_one"\`.

This validates the provenance of an EpiExposure-derived curve without
incorrectly treating curve-evaluation rows as original exposure-history
observations.

\## Uncertainty

\`epi_sensitivity()\` does not generate coefficient/posterior draws. If
the supplied data contain draw-specific curves, include the draw
identifier in \`scenario_var\` so each draw is differentiated
independently. The function does not silently pool repeated \`x\`
positions across draws.
