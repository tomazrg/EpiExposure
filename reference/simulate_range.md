# Generate DLNM simulation scenarios (scenario or profile mode)

Creates structured input for simulate_scenarios(), supporting both
discrete scenarios and continuous profiles (ranges).

## Usage

``` r
simulate_range(
  lag_windows,
  vary = list(),
  fixed = list(),
  scenario_type = c("grid", "paired"),
  mode = c("scenario", "profile"),
  scenario_names = NULL,
  scenario_var = NULL
)
```

## Arguments

- lag_windows:

  Output from define_lag_windows()

- vary:

  Named list of variables to vary (vectors of values)

- fixed:

  Named list of fixed values

- scenario_type:

  "grid" (all combinations) or "paired"

- mode:

  "scenario" (default) or "profile"

- scenario_names:

  Optional custom names

- scenario_var:

  Optional variable used to define scenario naming (recommended)

## Value

Named list of scenarios
