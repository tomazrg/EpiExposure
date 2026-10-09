# Define DLNM exposure templates

Defines the training cross-basis parameterization for one or more
exposure variables from long-format epidemic profiles.

## Usage

``` r
define_exposures(
  data,
  vars,
  max_lag,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns"
)
```

## Arguments

- data:

  Non-empty long-format data frame containing \`epi_id\`, \`time\`, and
  every exposure named in \`vars\`.

  Rows do not need to be pre-sorted. They are ordered internally by
  \`epi_id\` and \`time\`.

  Each epidemic must form one complete exposure profile containing
  exactly \`max_lag + 1\` time points. When more than one time point is
  present, spacing must be constant within epidemics and identical
  across epidemics. Missing exposure values are not allowed in
  EpiExposure v1. Because all exposures are finite columns of the same
  validated rows, every fitted exposure uses the same temporal support
  and profile length.

- vars:

  Unique character vector naming exposure variables in \`data\`.

- max_lag:

  Non-negative integer maximum lag. The legacy form \`c(0, L)\` is also
  accepted, but EpiExposure v1 requires the lag range to start at 0.

- df_var:

  Positive integer controlling the exposure-response basis when
  \`fun_var\` is \`"ns"\` or \`"bs"\`. For \`fun_var = "poly"\`,
  \`df_var\` is used as the polynomial degree. It is ignored for
  \`fun_var = "lin"\`.

- df_lag:

  Positive integer controlling the lag-response basis when \`fun_lag\`
  is \`"ns"\` or \`"bs"\`. For \`fun_lag = "poly"\`, \`df_lag\` is used
  as the polynomial degree. It is ignored for \`fun_lag = "lin"\`.

- fun_var:

  Character exposure-basis function. Supported unpenalized EpiExposure
  v1 choices are \`"ns"\`, \`"bs"\`, \`"poly"\`, and \`"lin"\`.

- fun_lag:

  Character lag-basis function. Supported unpenalized EpiExposure v1
  choices are \`"ns"\`, \`"bs"\`, \`"poly"\`, and \`"lin"\`.

  Penalized dlnm bases such as \`"ps"\` are not supported by the current
  engine-agnostic EpiExposure fitting contract. A P-spline
  transformation is not, by itself, a penalized DLNM: the corresponding
  penalty must also be supplied during model fitting (for example
  through \`dlnm::cbPen()\` and \`mgcv::gam(..., paraPen = ...)\`, or
  through a dedicated cross-basis smooth). EpiExposure v1 therefore
  stops explicitly rather than fitting an unintentionally unpenalized
  \`"ps"\` basis.

## Value

A named list of \`dlnm::crossbasis\` training templates, one per
exposure. The list carries:

\- \`attr(x, "spec")\`: effective per-exposure basis specification; -
\`attr(x, "epiexposure_vars")\`: exposure names; - \`attr(x,
"epiexposure_max_lag")\`: the common fitted maximum lag; - \`attr(x,
"epiexposure_history_length")\`: the exact required exposure profile
length, \`max_lag + 1\`; - \`attr(x, "epiexposure_history_contract")\`:
\`"all_fitted_exposures_same_exact_max_lag_plus_one"\`; - \`attr(x,
"epiexposure_time_step")\`: the common temporal spacing, or \`NA\` when
\`max_lag = 0\` because one observation has no estimable spacing; -
\`attr(x, "epiexposure_template_contract")\`:
\`"training_grouped_crossbasis"\`; - \`attr(x,
"epiexposure_profile_order")\`: \`"chronological"\`.

Each individual cross-basis also stores:

\- \`attr(cb, "cb_colnames")\`: native dlnm cross-basis column names; -
\`attr(cb, "epiexposure_canonical_cb_colnames")\`: the canonical future
design names \`cb\_\<variable\>\_\<index\>\`; - \`attr(cb,
"epiexposure_variable")\`: exposure name.

## Details

\`define_exposures()\` is the first step in the standard EpiExposure
workflow:

“\` templates \<- define_exposures(...) design \<- build_design(data,
templates) fit \<- fit_epidlnm(...) “\`

The function estimates \*\*basis definitions\*\*, not disease effects.
For each exposure it uses all training epidemics jointly to determine
any data-dependent exposure-basis features (for example spline knots and
boundary knots), while \`dlnm::crossbasis(..., group = ...)\` keeps the
individual epidemic profiles as independent time series so lags never
cross from one epidemic into another.

\## Independent epidemic series

Earlier EpiExposure code concatenated epidemic profiles and inserted
\`max_lag\` missing values between them before calling \`crossbasis()\`.
The current implementation uses the native \`group\` argument of
\`dlnm::crossbasis()\` instead.

This is the intended dlnm representation for several independent time
series. After internal sorting, every epidemic remains consecutive and
complete, while exposure values from all training epidemics still
contribute to the common exposure-basis definition. Thus, for a spline
exposure basis, data-dependent knots are estimated from the pooled
\*\*training exposure distribution\*\*, but lagged profiles never cross
epidemic boundaries.

\## Temporal requirements

A vector passed to \`dlnm::crossbasis()\` represents one complete
exposure profile. \`define_exposures()\` therefore rejects:

\- duplicated time values within an epidemic; - irregular spacing within
an epidemic when more than one time point exists; - different time steps
across epidemics; - epidemics with either fewer \*\*or more\*\* than
\`max_lag + 1\` observations.

If the common maximum lag is \`L\`, every epidemic must contain exactly

\$\$ L + 1 \$\$

observations. Profiles are never truncated to a trailing window, padded,
or silently realigned. Because all fitted exposures are columns of the
same long-format rows, all variables necessarily have the same temporal
length and support within each epidemic.

For \`max_lag = 0\`, the exact profile length is one observation. In
that special case no time interval exists to estimate, so
\`epiexposure_time_step\` is stored as \`NA\`.

These are the same temporal assumptions enforced downstream by
\`build_design()\`.

\## Exposure and lag intercepts

The exposure basis is constructed without an intercept. This follows the
standard DLNM identifiability requirement and avoids rank deficiency
when a model intercept is present.

The lag basis includes an intercept. This is the standard cross-basis
parameterization for the lag dimension and allows an effect distributed
over the full fitted lag window.

\## Effective versus requested basis arguments

\`dlnm::crossbasis()\` calls \`onebasis()\` internally and may normalize
or modify basis arguments. Accordingly, EpiExposure stores the
\*\*effective\*\* \`argvar\`, \`arglag\`, and lag attributes returned by
the completed cross-basis, rather than assuming that the original
request is the final parameterization.

These effective attributes are the authoritative training template
reused by \`build_design()\`, \`predict_outcomes()\`, cross-validation,
scenario functions, and effect summaries. Prediction data never redefine
knots or boundary knots.

\## Scope of \`df_var\` and \`df_lag\`

For spline bases, the user-supplied degrees of freedom determine the
requested flexibility, but the returned \`spec\` records the effective
basis arguments. For \`"poly"\`, the corresponding \`df\_\*\` argument
is interpreted as polynomial degree. For \`"lin"\`, the corresponding
\`df\_\*\` value has no effect.

For predictors with highly skewed distributions or many repeated values,
\`splines::ns()\` may issue a knot-placement warning when interior knots
coincide with boundary values. This adjustment is handled automatically
and does not prevent model fitting.
