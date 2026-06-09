# Getting started with EpiExposure

## What EpiExposure does?

`Introduzir o usuário ao workflow básico: da exposição ao modelo.`

``` r
remotes::install_github("tomazrg/EpiExposure")
library(EpiExposure)

# carregar dados
data("TargetSpot")

simula
```

``` r
# definir exposição (DLNM structure)
cb <- define_exposure(
  wx_long = TargetSpot,
  vars = c("tmean", "vpd"),
  lag_max = 85
)
```

``` r
# construir design matrix
dat <- build_design_matrix(
  wx_long = TargetSpot,
  cb_templates = cb,
  lag_max = 85
)
```

``` r
# preparar resposta
dat <- prepare_response(dat, y_var = "severity", family = "beta")
```

``` r
# ajustar modelo
fit <- fit_epidlnm(
  dat = dat,
  model_engine = "glmmTMB",
  family = "beta",
  random_effect = "epi_id",
  epiexposure_spec = attr(cb, "spec")
)
```

`DLNM não é só modelo — é estrutura epidemiológica`
`define_exposure() define como a exposição atua no tempo`
`fit_epidlnm() estimula a relação exposição–doença`

``` r
hist(dat$y_model, main="Distribution of disease response", col="steelblue")
```
