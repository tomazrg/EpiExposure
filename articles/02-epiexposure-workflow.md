# An epidemiological workflow for DLNM interpretation

## Objective

`Explicar o pipeline conceitual` `Aqui NÃO precisa de muito código`

`Exposição no tempo → Estrutura lag → Modelo → Interpretação`

DLNM como framework epidemiológico, não só estatístico por que epidemias
têm memória temporal

``` r
title("DLNM basis for temperature")
```

``` r
library(EpiExposure)
library(dlnm)

sim_data$dpp = sim_data$dpp-1


t = define_exposure(
  wx_long = TargetSpot,
  vars = c("tmean","vpd","rain"),
  lag_max = 85,
  df_var = 4,
  df_lag = 4,
  fun_var = "ns",
  fun_lag = "ns"
)
```

``` r
t2 <- build_design_matrix(sim_data, t, lag_max = 85)
```

``` r
t3 <- prepare_response(sim_data, y_var = "y", family_choice = "beta")
```

``` r
y_df <- t3 |>
  group_by(epi_id) |>
  summarise(
    y_model = first(y_model),
    .groups = "drop"
  )
```

``` r
dat_model <- left_join(t2, y_df, by = "epi_id")
```

``` r
fit <- fit_epidlnm(
  dat = dat_model,
  family = "beta",
  model_engine = "glmmTMB"
)

summary(fit)
```
