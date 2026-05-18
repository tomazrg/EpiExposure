# Predict full DLNM exposure-lag-response surface

Predict full DLNM exposure-lag-response surface

## Usage

``` r
predict_surface(
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

  ...

- wx_long:

  ...

- var:

  ...

- lag_max:

  ...

- df_var:

  Degrees of freedom (exposure dimension)

- df_lag:

  Degrees of freedom (lag dimension)

- fun_var:

  Basis function ("ns","bs","poly","lin")

- fun_lag:

  Basis function ("ns","ps","lin")

- ref:

  Reference exposure definition

- probs:

  Quantiles defining exposure grid

## Value

crosspred object
