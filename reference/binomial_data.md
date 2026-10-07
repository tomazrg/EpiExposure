# Simulated binomial epidemiological dataset

A simulated longitudinal dataset containing complete environmental
exposure histories and a binary plant disease outcome. The dataset is
intended for demonstrating binomial models in EpiExposure.

## Usage

``` r
binomial_data
```

## Format

A data frame with rows representing temporal observations within
epidemic exposure histories and the following variables:

- epi_id:

  Unique numeric identifier for each epidemic.

- block:

  Experimental block identifier.

- time:

  Chronological time index within each epidemic history.

- tmean:

  Daily mean temperature in degrees Celsius.

- rain:

  Daily rainfall in millimeters.

- wetness:

  Daily leaf wetness duration in hours.

- y:

  Binary disease outcome coded as 0 or 1 and repeated across the rows of
  each epidemic history.

## Source

Simulated for the EpiExposure package.

## Details

The response is repeated across the temporal rows belonging to the same
epidemic because each complete exposure history is associated with one
final epidemic-level disease outcome.

## Examples

``` r
data("binomial_data", package = "EpiExposure")
head(binomial_data)
#>   epi_id block time    tmean rain  wetness y
#> 1      1   B01    0 35.00193    0 5.599322 0
#> 2      1   B01    1 34.64017    0 6.490799 0
#> 3      1   B01    2 33.77685    0 5.479461 0
#> 4      1   B01    3 33.42277    0 3.781691 0
#> 5      1   B01    4 33.27345    0 3.357283 0
#> 6      1   B01    5 33.25479    0 4.610523 0
```
