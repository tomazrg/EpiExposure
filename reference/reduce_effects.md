# Reduce DLNM effects to one dimension (article-consistent)

Reduce DLNM effects to one dimension (article-consistent)

## Usage

``` r
reduce_effects(
  fit,
  wx_long,
  var,
  lag_max,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns",
  type = c("overall", "lag", "var"),
  value = NULL,
  scale = c("percent", "response", "link")
)
```

## Arguments

- fit:

  Fitted model object

- wx_long:

  Long-format weather data

- var:

  Exposure variable (ex: "tmax")

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

- type:

  "overall", "lag", "var"

- value:

  Required for type = "lag" or "var"

- scale:

  "link", "response", "percent"

## Value

data.frame
