## -----------------------------------------------------------------------------
library(EpiExposure)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
data("epi_data")

cb <- define_exposures(
  data = epi_data,
  vars = c(
    "tmean",
    "rain",
    "wetness"
  ),
  max_lag = 85,
  df_var = 3,
  df_lag = 2
)

dat <- build_design(
  data = epi_data,
  cb_templates = cb,
  max_lag = 85
)

dat <- prepare_response(
  data = dat,
  response = "y",
  family = "beta"
)

fit <- fit_epidlnm(
  data = dat,
  model_engine = "glmmTMB",
  family = "beta",
  epiexposure_spec = attr(cb,"spec")
)

## -----------------------------------------------------------------------------
periods <- define_periods(
  max_lag = 85,
  cuts = c(
    20,
    40,
    60,
    85
  )
)

## -----------------------------------------------------------------------------
scenario_grid <- simulate_ranges(
  periods = periods,
  vary = list(
    tmean = c(20,34),
    rain = seq(0,15,1),
    wetness = c(2,10)
  ),
  scenario_type = "grid",
  scenario_var = "rain"
)

## -----------------------------------------------------------------------------
predictions <- simulate_scenarios(
  fit = fit,
  scenarios = scenario_grid,
  data = epi_data,
  uncertainty = TRUE,
  output = "summary",
  n_samples = 100
)

## -----------------------------------------------------------------------------

head(predictions)


## -----------------------------------------------------------------------------
slope = 49.3
intercept = 9691

## -----------------------------------------------------------------------------
losses <- simulate_losses(
  data = predictions,
  y = "prediction",
  slope = slope,
  intercept = intercept,
  attainable_yield = seq(
    4000,
    12000,
    by = 50
  ),
  price = 300
)

## -----------------------------------------------------------------------------

head(losses)


## -----------------------------------------------------------------------------
losses %>%
  ggplot(
    aes(
      rain,
      attainable_yield,
      fill = yield_loss
    )
  ) +
  geom_raster() +
  facet_wrap(
    ~ tmean + wetness,
    ncol = 4
  ) +
  theme_bw() +
  scale_fill_viridis_c() +
  labs(
    x = "Rainfall (mm)",
    y = "Attainable yield (kg ha⁻¹)",
    fill = "Yield loss\n(kg ha⁻¹)"
  )

## -----------------------------------------------------------------------------

price = 300


## -----------------------------------------------------------------------------
losses %>%
  ggplot(
    aes(
      rain,
      attainable_yield,
      fill = economic_loss
    )
  ) +
  geom_raster() +
  facet_wrap(
    ~ tmean + wetness,
    ncol = 4
  ) +
  theme_bw() +
  scale_fill_viridis_c() +
  labs(
    x = "Rainfall (mm)",
    y = "Attainable yield (kg ha⁻¹)",
    fill = "Economic loss\n(USD ha⁻¹)"
  )

## -----------------------------------------------------------------------------
low_pressure <- losses %>%
  filter(
    rain == min(rain),
    wetness == min(wetness)
  )

## -----------------------------------------------------------------------------
high_pressure <- losses %>%
  filter(
    rain == max(rain),
    wetness == max(wetness)
  )

## -----------------------------------------------------------------------------
bind_rows(
  mutate(low_pressure,
         scenario = "Low pressure"),
  mutate(high_pressure,
         scenario = "High pressure")
) %>%
  ggplot(
    aes(
      economic_loss,
      fill = scenario
    )
  ) +
  geom_histogram(
    alpha = 0.6,
    position = "identity"
  ) +
  theme_bw()

## -----------------------------------------------------------------------------
losses %>%
  ggplot(
    aes(
      attainable_yield,
      economic_loss
    )
  ) +
  geom_smooth() +
  theme_bw()

