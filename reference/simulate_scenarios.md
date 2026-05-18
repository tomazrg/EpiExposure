# Simulate epidemiological DLNM scenarios

Simulate epidemiological DLNM scenarios

## Usage

``` r
simulate_scenarios(
  fit,
  dat,
  wx_long,
  lag_windows,
  scenarios,
  lag_max,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns",
  ref_vals = NULL,
  pop_level = TRUE
)
```

## Arguments

- fit:

  Fitted model

- dat:

  Original design matrix

- wx_long:

  Long-format weather data

- lag_windows:

  Lag windows

- scenarios:

  Output from simulate_range()

- lag_max:

  Maximum lag

- df_var:

  Degrees of freedom (exposure)

- df_lag:

  Degrees of freedom (lag)

- fun_var:

  Basis function for exposure

- fun_lag:

  Basis function for lag

- ref_vals:

  Reference values

- pop_level:

  Logical
