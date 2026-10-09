# Simulated negative-binomial epidemiological dataset

A simulated longitudinal dataset containing complete environmental
exposure profiles and an overdispersed count disease outcome. The
dataset is intended for demonstrating quadratic-variance
negative-binomial models in EpiExposure.

## Usage

``` r
nb2_data
```

## Format

A data frame with rows representing temporal observations within
epidemic exposure profiles and the following variables:

- epi_id:

  Unique numeric identifier for each epidemic.

- block:

  Experimental block identifier.

- time:

  Chronological time index within each epidemic profile.

- tmean:

  Daily mean temperature in degrees Celsius.

- wetness:

  Daily leaf wetness duration in hours.

- rain:

  Daily rainfall in millimeters.

- y:

  Non-negative count disease outcome repeated across the rows of each
  epidemic profile.

## Source

Simulated for the EpiExposure package.

## Details

The response is repeated across the temporal rows belonging to the same
epidemic because each complete exposure profile is associated with one
final epidemic-level disease outcome.

## Examples

``` r
data("nb2_data", package = "EpiExposure")
head(nb2_data)
#>   epi_id block time    tmean  wetness      rain y
#> 1      1   B01    0 28.51652 3.304398  9.234377 7
#> 2      1   B01    1 25.78990 5.600390 19.225955 7
#> 3      1   B01    2 26.86089 7.060567  0.000000 7
#> 4      1   B01    3 25.19746 5.064939  8.700396 7
#> 5      1   B01    4 25.80954 5.634410  0.000000 7
#> 6      1   B01    5 25.28982 6.813331  5.096005 7
```
