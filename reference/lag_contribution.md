# Quantify lag contributions to the cumulative DLNM effect

Decomposes the cumulative DLNM effect into relative contributions of
individual lags based on daily lag-specific effects.

## Usage

``` r
lag_contribution(daily_df, lag_window = NULL, absolute = TRUE)
```

## Arguments

- daily_df:

  Output from summarise_effects(scale = "daily"). Must contain columns:
  lag, effect.

- lag_window:

  Optional vector of length 2 specifying the lag interval c(lag_start,
  lag_end). If NULL, the full lag range is used.

- absolute:

  Logical. If TRUE (default), contributions are calculated using
  absolute effects to avoid sign cancellation.

## Value

A data.frame with lag-specific contributions to the cumulative effect.
