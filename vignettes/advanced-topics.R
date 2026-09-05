## -----------------------------------------------------------------------------
fit_glmmTMB <- fit_epidlnm(
  data = dat,
  model_engine = "glmmTMB",
  family = "beta",
  epiexposure_spec = attr(cb,"spec")
)


## -----------------------------------------------------------------------------
fit_gam <- fit_epidlnm(
  data = dat,
  model_engine = "gam",
  family = "beta",
  epiexposure_spec = attr(cb,"spec")
)


## -----------------------------------------------------------------------------
fit_brms <- fit_epidlnm(
  data = dat,
  model_engine = "brms",
  family = "beta",
  epiexposure_spec = attr(cb, "spec"))


## -----------------------------------------------------------------------------
fit_inla <- fit_epidlnm(
  data = dat,
  model_engine = "INLA",
  family = "beta",
  epiexposure_spec = attr(cb, "spec"))


## -----------------------------------------------------------------------------
output = "summary"


## -----------------------------------------------------------------------------
output = "samples"


## -----------------------------------------------------------------------------
epi_l <- summarise_effects(
  fit = fit,
  data = epi_data,
  uncertainty = TRUE,
  output = "summary"
)


## -----------------------------------------------------------------------------
epi_l_samples <- summarise_effects(
  fit = fit,
  data = epi_data,
  uncertainty = TRUE,
  n_samples = 100,
  output = "samples"
)


## -----------------------------------------------------------------------------
pred <- predict_outcome(
  fit = fit,
  profiles = scenario,
  uncertainty = TRUE,
  n_samples = 100,
  output = "samples"
)


## -----------------------------------------------------------------------------
ggplot(
  pred,
  aes(prediction)
) +
  geom_histogram(
    bins = 20
  ) +
  theme_bw()


## -----------------------------------------------------------------------------
check_identifiability(
  data = epi_data,

  var = c(
    "tmean",
    "rain",
    "wetness"
  ),

  df_var = 3,
  df_lag = 2,

  max_lag = 85
)


## -----------------------------------------------------------------------------
df_var = 3
df_lag = 2


## -----------------------------------------------------------------------------
df_var = 4
df_lag = 4


## -----------------------------------------------------------------------------
df_var = 2
df_lag = 2


## -----------------------------------------------------------------------------
splines::ns()


## -----------------------------------------------------------------------------
set.seed(123) # verificar se nao seria o "seed" usado nas funcoes ao inves de "set.seed(123")

