# Compare accumulated DLNM effects between periods

Compares accumulated effects across epidemiological periods using a
reference period. Differences and ratios are computed relative to the
reference.

## Usage

``` r
compare_periods(accumulated_df, period_ref = 1, eps = 1e-10)
```

## Arguments

- accumulated_df:

  Output from summarise_effects(scale = "accumulated"). Must contain
  columns: value, period, effect.

- period_ref:

  Reference period (numeric index or character label, e.g. 1 or "W1")

- eps:

  Small constant to avoid division by zero in ratio calculation.

## Value

A data.frame with comparisons relative to the reference period.
