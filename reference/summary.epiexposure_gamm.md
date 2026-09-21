# Summarize an EpiExposure GAMM fit

Summarizes a generalized additive mixed model fitted through
\`fit_epidlnm()\` with \`model_engine = "gamm"\`.

## Usage

``` r
# S3 method for class 'epiexposure_gamm'
summary(object, component = c("gam", "lme", "both"), ...)
```

## Arguments

- object:

  A fitted EpiExposure GAMM object.

- component:

  Character. Component to summarize:

  \- \`"gam"\` returns the fixed/population model summary; - \`"lme"\`
  returns the mixed-model summary; - \`"both"\` returns both summaries.

- ...:

  Additional arguments passed to the corresponding summary method.

## Value

A model summary object. When \`component = "both"\`, a list containing
the GAM and LME summaries is returned.

## Details

The GAMM object contains separate \`gam\` and \`lme\` components. By
default, the method returns the \`gam\` summary containing the
fixed/population coefficients and their inferential statistics.
