# Check DLNM cross-basis identifiability

Evaluates whether a DLNM cross-basis matrix is full rank.

## Usage

``` r
check_identifiability(
  wx_long,
  var,
  lag_max,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns"
)
```

## Arguments

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

## Value

TRUE/FALSE
