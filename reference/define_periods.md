# Define epidemiological lag periods

Creates lag periods to summarise cumulative DLNM effects.

## Usage

``` r
define_periods(lag_max, cuts = NULL, prefix = "W")
```

## Arguments

- lag_max:

  Maximum lag (integer)

- cuts:

  Optional numeric vector of cut points defining periods

- prefix:

  Character string used to label periods (default = "W")

## Value

data.frame with columns: - period - lag_start - lag_end
