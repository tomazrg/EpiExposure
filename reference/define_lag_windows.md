# Define epidemiological lag windows

Creates lag windows to summarise cumulative DLNM effects. If no cut
points are provided, a single window from lag 0 to lag_max is returned.

## Usage

``` r
define_lag_windows(lag_max, cuts = NULL, prefix = "W")
```

## Arguments

- lag_max:

  Maximum lag (integer)

- cuts:

  Optional numeric vector of cut points defining windows

- prefix:

  Character string used to label windows (default = "W")

## Value

data.frame with columns: - window_id (character labels: e.g., "W1",
"W2", ...) - lag_start - lag_end
