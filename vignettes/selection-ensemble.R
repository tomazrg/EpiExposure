## -----------------------------------------------------------------------------
library(EpiExposure)
library(dplyr)

data("epi_data")


## -----------------------------------------------------------------------------
library(future)

future::plan(
  future::multisession,
  workers = 4
)

