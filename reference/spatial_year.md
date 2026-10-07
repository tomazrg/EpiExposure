# Simulated epidemiological data with independent spatial fields by year

A simulated longitudinal dataset containing complete environmental
exposure histories for fixed field locations observed across multiple
years. The dataset is intended for demonstrating distributed
exposure-lag models with separate spatial fields by year.

## Usage

``` r
spatial_year
```

## Format

A data frame with rows representing temporal observations within
epidemic exposure histories and the following variables:

- epi_id:

  Unique identifier combining field location and year.

- location_id:

  Identifier for the field location.

- year:

  Year identifying the independent spatial field.

- block:

  Experimental block identifier within a year.

- block_id:

  Identifier combining year and experimental block.

- time:

  Chronological time index within each epidemic history.

- x_coord:

  First numeric spatial coordinate of the field location.

- y_coord:

  Second numeric spatial coordinate of the field location.

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

The response and location-level descriptors are repeated across the
temporal rows belonging to the same epidemic because each complete
exposure history is associated with one final epidemic-level disease
outcome.

## Examples

``` r
data("spatial_year", package = "EpiExposure")
head(spatial_year)
#>      epi_id location_id year block block_id time    x_coord  y_coord    tmean
#> 1 L001_2022        L001 2022   B01 2022_B01    0 -0.6698125 0.622697 24.38006
#> 2 L001_2022        L001 2022   B01 2022_B01    1 -0.6698125 0.622697 24.53832
#> 3 L001_2022        L001 2022   B01 2022_B01    2 -0.6698125 0.622697 24.90186
#> 4 L001_2022        L001 2022   B01 2022_B01    3 -0.6698125 0.622697 26.04159
#> 5 L001_2022        L001 2022   B01 2022_B01    4 -0.6698125 0.622697 23.60107
#> 6 L001_2022        L001 2022   B01 2022_B01    5 -0.6698125 0.622697 25.79864
#>     wetness      rain y
#> 1  4.903199 12.904436 3
#> 2  6.781023  7.651138 3
#> 3 11.195548  0.000000 3
#> 4 10.530333  0.000000 3
#> 5 10.322681  0.000000 3
#> 6 12.433350  0.000000 3
```
