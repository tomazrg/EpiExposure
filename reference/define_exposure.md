# Define DLNM exposure templates with flexible splines

Define DLNM exposure templates with flexible splines

## Usage

``` r
define_exposure(
  wx_long,
  vars,
  lag_max,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns"
)
```

## Arguments

- wx_long:

  Long-format weather data (epi_id, dpp, variables)

- vars:

  Character vector of exposure variables

- lag_max:

  Maximum lag

- df_var:

  Degrees of freedom for exposure dimension

- df_lag:

  Degrees of freedom for lag dimension

- fun_var:

  Basis function for exposure ("ns","bs","poly","lin")

- fun_lag:

  Basis function for lag ("ns","ps","lin")

## Value

Named list of crossbasis templates
