# Summarise DLNM effects (article-consistent)

Summarise DLNM effects (article-consistent)

## Usage

``` r
summarise_effects(
  fit,
  wx_long,
  var,
  lag_max,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns",
  scale = c("daily", "accumulated"),
  lag_windows = NULL,
  probs = seq(0.05, 0.95, by = 0.01),
  ref = list(method = "median", value = NULL),
  effect_measure = c("percent", "ratio", "linear")
)
```

## Arguments

- fit:

  Fitted model

- wx_long:

  Long-format weather data

- var:

  Exposure variable name (e.g. "tmax")

- lag_max:

  Maximum lag

- df_var:

  Degrees of freedom for exposure

- df_lag:

  Degrees of freedom for lag

- fun_var:

  Exposure basis function

- fun_lag:

  Lag basis function

- scale:

  "daily" or "accumulated"

- lag_windows:

  Lag windows (required for accumulated)

- probs:

  Quantiles

- ref:

  Reference definition list

- effect_measure:

  "percent","ratio","linear"

## Value

data.frame
