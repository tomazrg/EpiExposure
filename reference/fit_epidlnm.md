# Fit DLNM inferential model

Fit DLNM inferential model

## Usage

``` r
fit_epidlnm(dat, model_engine, family, random_effect = NULL, ...)
```

## Arguments

- dat:

  Design matrix containing y_model and cb\_\* terms

- model_engine:

  Modeling engine
  ("glm","glmmTMB","gam","gamm","gls","spamm","brms","inla")

- family:

  Distribution

- random_effect:

  Random effect variable (optional)

- ...:

  Additional arguments passed to the modeling engine

## Value

Fitted model object
