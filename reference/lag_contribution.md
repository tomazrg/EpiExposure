# Quantify lag contributions to the cumulative DLNM effect

Supports: 1) Deterministic input 2) Summary input (effect with sd/CI) 3)
Samples input (with 'sample' column)

## Usage

``` r
lag_contribution(daily_df, lag_window = NULL, absolute = TRUE)
```

## Arguments

- daily_df:

  Output from summarise_effects(scale = "daily")

- lag_window:

  Optional lag interval c(start, end)

- absolute:

  Logical (default TRUE)

## Value

data.frame
