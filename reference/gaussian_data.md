# Simulated gaussian epidemiological dataset

A simulated longitudinal dataset containing complete environmental
exposure histories and a continuous plant disease outcome. The dataset
is intended for demonstrating gaussian models in EpiExposure.

## Usage

``` r
gaussian_data
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

- wetness:

  Daily leaf wetness duration in hours.

- rain:

  Daily rainfall in millimeters.

- y:

  Continuous disease outcome repeated across the rows of each epidemic
  history.

## Source

Simulated for the EpiExposure package.

## Details

The response is repeated across the temporal rows belonging to the same
epidemic because each complete exposure history is associated with one
final epidemic-level disease outcome.

## Examples

``` r
data("gaussian_data", package = "EpiExposure")
head(gaussian_data)
#>   epi_id block time    tmean  wetness     rain       y
#> 1      1   B01    0 29.16636 11.58879  0.00000 12.4501
#> 2      1   B01    1 27.30294 12.71218  0.00000 12.4501
#> 3      1   B01    2 26.86088 12.74422  0.00000 12.4501
#> 4      1   B01    3 26.29974 11.92657 10.12782 12.4501
#> 5      1   B01    4 26.67772 13.41525  0.00000 12.4501
#> 6      1   B01    5 26.63974 15.59525  0.00000 12.4501
```
