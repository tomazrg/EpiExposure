# Simulated poisson epidemiological dataset

A simulated longitudinal dataset containing complete environmental
exposure histories and a count plant disease outcome. The dataset is
intended for demonstrating poisson models in EpiExposure.

## Usage

``` r
poisson_data
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

  Non-negative count disease outcome repeated across the rows of each
  epidemic history.

## Source

Simulated for the EpiExposure package.

## Details

The response is repeated across the temporal rows belonging to the same
epidemic because each complete exposure history is associated with one
final epidemic-level disease outcome.

## Examples

``` r
data("poisson_data", package = "EpiExposure")
head(poisson_data)
#>   epi_id block time    tmean  wetness      rain y
#> 1      1   B01    0 25.32685 14.55902 11.629805 5
#> 2      1   B01    1 25.72697 13.35332  0.000000 5
#> 3      1   B01    2 24.88694 13.00895  8.036259 5
#> 4      1   B01    3 26.04031 13.22694  0.000000 5
#> 5      1   B01    4 25.94690 15.30767  6.709351 5
#> 6      1   B01    5 24.85839 17.41708  8.826920 5
```
