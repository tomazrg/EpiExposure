# Identify critical lags based on daily DLNM effects

Supports three input types: 1) Deterministic (no uncertainty) 2) Summary
(with sd / CI) 3) Samples (posterior / simulated draws)

## Usage

``` r
identify_critical_lags(daily_df, metric = c("max", "mean", "absmean"))
```

## Arguments

- daily_df:

  Output from summarise_effects(scale = "daily"). Must contain at least:
  lag, effect. Optional columns: sample, var, effect_sd, effect_lower,
  effect_upper

- metric:

  Summary metric used to quantify lag importance. One of "max", "mean",
  or "absmean".

## Value

A data.frame with lag importance ranking
