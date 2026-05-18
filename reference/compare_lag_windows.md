# Compare accumulated DLNM effects between lag windows

Compares accumulated effects across epidemiological lag windows using a
reference window. Differences and ratios are computed relative to the
reference.

## Usage

``` r
compare_lag_windows(accumulated_df, window_ref = 1, eps = 1e-10)
```

## Arguments

- accumulated_df:

  Output from summarise_effects(scale = "accumulated"). Must contain
  columns: value, window, effect.

- window_ref:

  Reference window (numeric index or character label, e.g. 1 or "W1")

- eps:

  Small constant to avoid division by zero in ratio calculation.

## Value

A data.frame with pairwise comparisons relative to the reference window.
