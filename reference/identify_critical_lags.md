# Identify critical lags based on daily DLNM effects

Identifies lags with the strongest disease association by summarising
daily DLNM effects across exposure levels.

## Usage

``` r
identify_critical_lags(daily_df, metric = c("max", "mean", "absmean"))
```

## Arguments

- daily_df:

  Output from summarise_effects(scale = "daily"). Must contain columns:
  lag, effect.

- metric:

  Summary metric used to quantify lag importance. One of "max", "mean",
  or "absmean".

## Value

A data.frame with one row per lag and summary statistics used to
identify critical lags.
