# Fit DLNM inferential model

Fit DLNM inferential model

## Usage

``` r
fit_epidlnm(
  dat,
  model_engine,
  family,
  random_effect = NULL,
  epiexposure_spec = NULL,
  basis_objects = NULL,
  ...
)
```

## Arguments

- dat:

  Design matrix containing y_model and cb\_\* terms

- model_engine:

  Modeling engine
  ("glm","glmmTMB","gam","gamm","gls","spamm","brms","inla","bdlnm")

- family:

  Distribution (engine-specific; can be character or family object)

- random_effect:

  Random effect variable (optional)

- epiexposure_spec:

  Optional list describing how crossbasis was built. Recommended
  structure (by variable name): list( tmean = list(lag_max=85,
  argvar=list(...), arglag=list(...)), vpd = list(lag_max=85,
  argvar=list(...), arglag=list(...)) ) If NULL, the model is fit
  normally but downstream prediction from profiles will not be
  available.

- basis_objects:

  Optional named list of original crossbasis/onebasis objects. Required
  when model_engine = "bdlnm", because bdlnm::bdlnm() needs the basis
  objects explicitly in the formula environment.

- ...:

  Additional arguments passed to the modeling engine

## Value

Fitted model object (same class as before), with extra attributes
