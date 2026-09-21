# Model Selection and Ensemble Forecasting

## Introduction

Distributed Lag Nonlinear Models (DLNMs) require several modeling
decisions, including:

- exposure variables;
- lag duration;
- exposure spline complexity;
- lag spline complexity;
- regression engine.

Different specifications may lead to different predictive performance.

`EpiExposure` provides tools for systematically exploring model
configurations, ranking competing models, and combining predictions
through ensemble forecasting.

This section demonstrates:

- automated DLNM model search;
- predictive validation;
- model ranking;
- ensemble forecasting;
- comparison of ensemble strategies.

## Load example data

``` r

library(EpiExposure)
library(dplyr)
library(ggplot2)
```

``` r

data("epi_data")

epi_data
```

    ## # A tibble: 44,720 × 6
    ##    epi_id  time tmean  rain wetness     y
    ##     <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
    ##  1      1     0  22.0  0       8.83 0.330
    ##  2      1     1  21.5  0      11.9  0.330
    ##  3      1     2  21.5  5.38   12.5  0.330
    ##  4      1     3  21.8 14.6    12.3  0.330
    ##  5      1     4  20.8  0      11.4  0.330
    ##  6      1     5  23.0  8.24    7.98 0.330
    ##  7      1     6  23.9  4.69   10.8  0.330
    ##  8      1     7  23.1 12.3     8.58 0.330
    ##  9      1     8  25.3  4.14   10.6  0.330
    ## 10      1     9  26.0  1.91   14.1  0.330
    ## # ℹ 44,710 more rows

## Why model selection matters

Exposure-lag-response relationships are highly sensitive to model
specification.

For example, increasing lag flexibility may improve fit but increase
overfitting.

Likewise, adding additional predictors may improve prediction in some
pathosystems while reducing generalizability in others.

`EpiExposure` can evaluate predefined combinations of candidate DLNM
specifications using a consistent fitting and validation workflow.

## Searching candidate models

``` r

library(future)

future::plan(
  future::multisession,
  workers = 4
)
```

The appropriate number of workers depends on the available hardware and
local computing policies. Users should avoid requesting more workers
than their system can support.

``` r

future::plan(future::sequential ) 
```

### Run the model search

``` r

best_models <- find_bestfit(
  data = epi_data,
  response = "y",
  group = "epi_id",
  time = "time",
  vars = c("tmean","rain","wetness"),
  max_lag = 85,
  df_var_grid = c(2,3,4),
  df_lag_grid = c(2,3,4),
  model_engine = "glmmTMB",
  family = "beta",
  keep_fits = TRUE)
```

Inspect results

``` r

head(best_models)
```

Depending on the selected validation settings, the resulting object can
include:

- a model identifier;
- the fitted predictor combination;
- exposure and lag spline settings;
- the modeling engine and response family;
- the Concordance Correlation Coefficient;
- the Root Mean Square Error;
- the Mean Absolute Error;
- bias-related metrics;
- the fitted model object when `keep_fits = TRUE`.

## Ranking models

### Top-performing models

The following example ranks models by decreasing CCC. This is one
possible ranking criterion, not a universal rule. Candidate models
should also be examined using error metrics, bias, convergence
diagnostics, model complexity, and the scientific plausibility of the
exposure-lag structure.

``` r

best_models %>%
  arrange(desc(CCC)) %>%
  select(
    model_id,
    vars,
    CCC,
    RMSE,
    MAE
  ) %>%
  head(10)
```

## Visualizing model performance

``` r

best_models %>%
  mutate(
    splines = paste(
      df_var,
      df_lag,
      sep = " × "
    )
  ) %>%
  ggplot(
    aes(
      Cb,
      CCC,
      color = splines
    )
  ) +
  geom_point(size = 3) +
  theme_bw()
```

## Understanding validation metrics

### Concordance Correlation Coefficient

CCC evaluates agreement between observed and predicted outcomes by
combining precision and accuracy.

Values closer to 1 indicate stronger concordance. However, CCC should be
interpreted together with error and bias metrics rather than used as the
only selection criterion.

### Root Mean Square Error

RMSE summarizes the magnitude of prediction errors while assigning
greater influence to larger errors because residuals are squared before
averaging.

Lower values indicate better predictive performance when models are
evaluated on the same response scale and validation data.

### Mean Absolute Error

MAE is the average absolute difference between observed and predicted
outcomes.

Lower values indicate better predictive performance. Because MAE remains
on the response scale, it is often easier to interpret directly than
RMSE.

## Inspecting specific models

Suppose the best model is:

``` r

top_model <- best_models |>
  dplyr::arrange(
    dplyr::desc(
      CCC
    )
  ) |>
  dplyr::slice(
    1L
  )
```

Retrieve the fitted model:

``` r

fit_best <- best_models$fit[[
  which.max(
    best_models$CCC
  )
]]

summary(
  fit_best
)
```

Review model diagnostics:

``` r

summary(fit_best)
```

## Why use ensembles?

No single model is guaranteed to be optimal.

Different DLNM structures frequently capture different aspects of the
underlying epidemiological process.

Combining predictions from multiple models often improves stability and
predictive robustness.

## Selecting ensemble candidates

The identifiers below are illustrative. Replace them with model
identifiers selected from the user’s own `best_models` object.

``` r

candidate_models <- c(
  63,
  23,
  44
)
```

Therefore, three candidate models (best ranked \[63\], and the two
lower-ranked models \[23, 44\]) were selected to be used on the ensemble
approaches, with CCC used as the performance metric for the weighted
ensemble.

## Unweighted ensemble

``` r

ens_unweighted <- ensemble_bestfit(
  bestfit = best_models,

  ensemble_scope = "model",

  method = "unweighted",

  model_id = candidate_models)
```

### Inspect results

``` r

ens_unweighted$ensemble_summary
```

## Weighted ensemble

``` r

ens_weighted <- ensemble_bestfit(
  bestfit = best_models,

  ensemble_scope = "model",

  method = "weighted",

  model_id = candidate_models,

  weight_metric = "CCC",

  weight_transform = "rank_inverse"
)
```

### Inspect results

``` r

ens_weighted$ensemble_summary
```

## Stacked ensemble

``` r

ens_stacked <- ensemble_bestfit(
  bestfit = best_models,

  ensemble_scope = "model",

  method = "stacked",

  model_id = candidate_models,

  stacking_model = "ridge",

  stack_objective = "regularized"
)
```

### Inspect results

``` r

ens_stacked$ensemble_summary
```

## Comparing ensemble strategies

``` r

ensemble_metrics <- dplyr::bind_rows(
  Unweighted = ens_unweighted$ensemble_summary,
  Weighted = ens_weighted$ensemble_summary,
  Stacked = ens_stacked$ensemble_summary,
  .id = "method"
) |>
  dplyr::select(
    method,
    CCC,
    RMSE,
    MAE
  )

ensemble_metrics
```

## Visualizing ensemble predictions

``` r

ensemble_all <- bind_rows(

  mutate(
    ens_unweighted$ensemble_predictions,
    method = "Unweighted"
  ),

  mutate(
    ens_weighted$ensemble_predictions,
    method = "Weighted"
  ),

  mutate(
    ens_stacked$ensemble_predictions,
    method = "Stacked"
  )
)
```

``` r

ggplot(
  ensemble_all,
  aes(
    observed,
    predicted_ensemble
  )
) +

  geom_point(alpha = 0.5) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red"
  ) +

  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "blue"
  ) +

  facet_wrap(
    ~method
  ) +

  theme_bw() +

  labs(
    x = "Observed",
    y = "Predicted"
  )
```

## Choosing an ensemble strategy

Each ensemble method has advantages.

**Unweighted ensemble**

- gives equal contribution to the selected models;
- is simple and transparent;
- avoids estimating additional model weights;
- can be affected by weak or redundant candidate models.

**Weighted ensemble**

- assigns different contributions to candidate models;
- can emphasize models with stronger validation performance;
- depends on the selected metric and weight transformation;
- does not guarantee improvement over the best individual model.

**Stacked ensemble**

- estimates combination coefficients from model predictions;
- can account for complementary predictive information;
- requires careful separation between weight estimation and final
  evaluation;
- can overfit when the validation sample is limited or when candidate
  predictions are highly similar.

Therefore, no ensemble method is universally superior.

An unweighted ensemble is useful as a transparent reference because
every candidate contributes equally.

A weighted ensemble is appropriate when the selected validation metric
provides a defensible basis for assigning model contributions.

A stacked ensemble can be useful when candidate models contain
complementary predictive information and the stacking weights are
estimated without using the same observations reserved for final
evaluation.

Regardless of the method, compare the ensemble against:

- the strongest individual candidate;
- a simple baseline model;
- alternative ensemble strategies;
- results on validation data not reused for weight estimation.

### Avoiding optimistic ensemble evaluation

Ensemble weights should be estimated from out-of-sample or
cross-validated predictions whenever possible. Estimating weights and
reporting performance on the same fitted predictions can produce
optimistic results.

The final evaluation should use observations that were not reused to
select candidate models or estimate ensemble weights. If a completely
independent test set is unavailable, nested or repeated cross-validation
may provide a more defensible assessment.

## Summary

This vignette introduced the `EpiExposure` model-selection and ensemble
framework.

The workflow includes:

- defining a candidate model space;
- fitting candidate DLNM specifications;
- evaluating predictive performance;
- ranking models using multiple criteria;
- retaining selected fitted models;
- combining predictions through unweighted, weighted, or stacked
  ensembles;
- comparing ensemble predictions with individual models and simple
  baselines.

Model selection and ensemble construction do not remove model
uncertainty. Their value depends on the candidate model space,
validation design, selection rules, and independence of the final
performance assessment.
