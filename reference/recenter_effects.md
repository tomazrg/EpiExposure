# Recenter DLNM effects using a new reference value

Recomputes the DLNM effect surface using a different centering value,
without refitting the model.

## Usage

``` r
recenter_effects(
  fit,
  wx_long,
  var,
  lag_max,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns",
  ref = list(method = "median", value = NULL),
  probs = seq(0.05, 0.95, by = 0.01)
)
```

## Arguments

- fit:

  Fitted model from fit_epidlnm()

- wx_long:

  Long-format weather data

- var:

  Exposure variable (e.g. "tmax")

- lag_max:

  Maximum lag

- df_var:

  Degrees of freedom (exposure)

- df_lag:

  Degrees of freedom (lag)

- fun_var:

  Basis ("ns","bs","poly","lin")

- fun_lag:

  Basis ("ns","ps","lin")

- ref:

  New reference definition

- probs:

  Quantiles for exposure grid

## Value

crosspred object
