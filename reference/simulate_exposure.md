# Simulate exposure history across lags (full or patterned)

Creates one or multiple exposure-lag profiles of length (lag_max + 1)

## Usage

``` r
simulate_exposure(
  lag_max,
  n = 1,
  mode = c("profile", "pattern"),
  pattern = c("random_lags", "alternating_lags", "block_lags"),
  fixed_value = NULL,
  n_lags = NULL,
  lag_range = NULL,
  alternating_values = NULL,
  block_lags = NULL,
  background = list(dist = "normal", mean = 0, sd = 1),
  bounds = NULL,
  seed = NULL,
  cumulative = FALSE
)
```

## Arguments

- lag_max:

  Integer. Maximum lag (profile length = lag_max + 1)

- n:

  Integer. Number of simulations (default = 1)

- mode:

  Character. "profile" or "pattern"

- pattern:

  Character. Pattern type

- fixed_value:

  Numeric. Value for pattern-controlled lags

- n_lags:

  Integer for "random_lags"

- lag_range:

  Integer vector of lags

- alternating_values:

  Numeric vector (\>=2)

- block_lags:

  Integer length-2

- background:

  List controlling background

- bounds:

  Optional numeric length-2

- seed:

  Optional seed

- cumulative:

  Logical. If TRUE, returns cumulative (cumsum) profile

## Value

\- If n = 1: list(profile, meta) - If n \> 1: list(profiles, meta)
