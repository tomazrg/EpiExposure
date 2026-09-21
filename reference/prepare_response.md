# Prepare and validate the epidemic-level response for modeling

Validates an epidemic-level response against the distribution family
selected by the user and creates the standardized \`y_model\` column
required by \`fit_epidlnm()\`.

## Usage

``` r
prepare_response(
  data,
  response = NULL,
  family,
  beta_scale = c("auto", "proportion", "percent"),
  y_var = NULL
)
```

## Arguments

- data:

  Non-empty data frame containing the response variable. In the standard
  workflow this is the output of \`build_design()\`.

- response:

  Character scalar naming the original response variable in \`data\`. In
  the standard workflow this is \`"y"\`.

  For compatibility with EpiExposure model-selection code, \`response\`
  may be omitted when \`y_var\` is supplied. If both are supplied they
  must identify the same column.

- family:

  Distribution family supplied as a supported character name or as a
  family-like object containing a valid \`family\` field.

  Supported EpiExposure v1 canonical families are:

  \- \`"beta"\`: continuous proportions strictly inside \`(0, 1)\`; -
  \`"binomial"\`: binary outcomes coded exactly \`0/1\`; -
  \`"poisson"\`: non-negative integer counts; - \`"negative_binomial"\`:
  NB2 non-negative integer counts; - \`"gaussian"\`: finite continuous
  numeric values; - \`"gamma"\`: strictly positive continuous values.

  Aliases such as \`"bernoulli"\`, \`"normal"\`, \`"negbin"\`,
  \`"nbinom"\`, and \`"nbinom2"\` are normalized to the corresponding
  canonical family.

  \`nbinom1\` / NB1 are rejected explicitly. EpiExposure v1 standardizes
  \`"negative_binomial"\` to the NB2 variance parameterization. Ordinal
  outcomes are not supported in EpiExposure v1.

- beta_scale:

  Character scalar controlling interpretation of a Beta response. One
  of:

  \- \`"auto"\` (default): values already within \`\[0, 1\]\` are
  interpreted as proportions; otherwise, if all values lie within \`\[0,
  100\]\` and at least one value exceeds 1, they are interpreted as
  percentages and divided by 100; - \`"proportion"\`: require the
  supplied response to be on the \`\[0, 1\]\` scale; - \`"percent"\`:
  require values on the \`\[0, 100\]\` scale and divide by 100.

  After scale handling, Beta regression still requires the \*\*open\*\*
  interval \`(0, 1)\`. Exact boundary values 0 or 1 are not altered
  silently.

- y_var:

  Optional compatibility alias for \`response\`. This is used by
  existing EpiExposure internal workflows. New user-facing code should
  prefer \`response\`.

## Value

A data frame containing all original columns plus \`y_model\`.

Existing attributes on \`data\` (including the design/template metadata
created by \`build_design()\`) are preserved. The following response
metadata are added:

\- \`"response_family_name"\`: canonical EpiExposure family; -
\`"response_name"\`: original response-column name; -
\`"response_original_scale"\`: response scale before preparation; -
\`"response_model_scale"\`: scale represented by \`y_model\`; -
\`"response_transform"\`: transformation applied, if any; -
\`"y_scale_mult"\`: multiplier that maps the Beta model scale back to
the original numerical scale (\`100\` after percent-to-proportion
conversion, otherwise \`1\`); - \`"response_nb_parameterization"\`:
\`"NB2"\` for negative-binomial responses, otherwise \`NULL\`; -
\`"response_contract"\`: \`"family_validated_response_v1"\`.

## Details

In the standard EpiExposure workflow:

“\` templates \<- define_exposures(...) design \<- build_design( data =
data, cb_templates = templates, include_response = TRUE ) prepared \<-
prepare_response( data = design, response = "y", family = "beta" ) fit
\<- fit_epidlnm( data = prepared, ... ) “\`

\`build_design(include_response = TRUE)\` and \`prepare_response()\`
have different responsibilities. \`build_design()\` only carries the
original epidemic-level outcome into the design matrix.
\`prepare_response()\` checks whether that outcome is compatible with
the requested statistical family and performs only the explicitly
documented scale conversion needed for Beta percentages.

\`prepare_response()\` does \*\*not\*\* choose a probability
distribution from the observed data and does not perform goodness-of-fit
model selection. The \`family\` is selected by the analyst and is then
validated here.

\## Why this function remains necessary after \`build_design()\`

With \`include_response = TRUE\`, \`build_design()\` ensures that one
original epidemic-level \`y\` value is carried into each design row. It
deliberately does not interpret that outcome statistically.

\`prepare_response()\` is the family-aware layer between the design
matrix and \`fit_epidlnm()\`. \`fit_epidlnm()\` expects a standardized
numeric \`y_model\` column and checks that the response-family metadata
agree with the family used for fitting.

\## Beta responses

Beta regression requires

\$\$0 \< y \< 1.\$\$

When \`beta_scale = "auto"\`, the function retains the previous
EpiExposure convenience of converting a clear 0–100 percentage scale to
proportions. The conversion is recorded in output metadata.

Earlier EpiExposure code silently replaced exact 0 and 1 by \`1e-5\` and
\`1 - 1e-5\`. That behavior is intentionally removed. Boundary
modification changes observed outcomes and the appropriate treatment
depends on the scientific data-generating process. Exact 0/1 values
therefore produce an explicit error rather than an undocumented
numerical adjustment.

\## Count responses

Poisson and negative-binomial responses must be non-negative integers.
Values are validated but never truncated or rounded to make them valid.

EpiExposure v1 uses NB2 whenever the canonical family is
\`"negative_binomial"\`. Explicit NB1 requests generate an error before
model fitting so different engines cannot silently use different
mean–variance relationships.

\## Binomial, Gaussian, and Gamma responses

Binomial responses must be numeric 0/1. Gaussian responses may be any
finite numeric values. Gamma responses must be strictly positive.

\## No distribution inference

This function validates a \*\*chosen\*\* family; it does not infer the
best family from histograms, moments, normality tests, overdispersion
tests, or other automatic rules. Distribution choice remains a modeling
decision based on the outcome definition and study design.
