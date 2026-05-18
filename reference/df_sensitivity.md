# Compute sensitivity of DLNM effects (analytical or finite derivative)

Calculates the derivative of the exposure–response relationship,
optionally including elasticity and critical point detection.

## Usage

``` r
df_sensitivity(
  df,
  x = "value",
  y = "prediction",
  scenario_var = "scenario",
  k = 10,
  method = c("analytical", "finite"),
  elasticity = TRUE,
  critical = TRUE,
  eps = 1e-06
)
```

## Arguments

- df:

  Data frame with columns: value and prediction/effect, and optionally a
  scenario column

- x:

  Name of predictor variable (default = "value")

- y:

  Name of response variable (default = "prediction")

- scenario_var:

  Optional grouping variable (default = "scenario")

- k:

  Basis dimension for GAM smoothing (default = 10)

- method:

  Derivative method: "analytical" (default) or "finite"

- elasticity:

  Logical, compute elasticity (default = TRUE)

- critical:

  Logical, detect critical points (default = TRUE)

- eps:

  Small epsilon for numerical derivative (default = 1e-6)

## Value

data.frame with sensitivity, optional elasticity and critical points
