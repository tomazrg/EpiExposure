# Simulated Gamma epidemiological dataset

A simulated longitudinal dataset containing complete environmental
exposure histories and a positive continuous plant disease outcome. The
dataset is intended for demonstrating Gamma models in EpiExposure.

## Usage

``` r
gamma_data
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

  Positive continuous disease outcome repeated across the rows of each
  epidemic history.

## Source

Simulated for the EpiExposure package.

## Details

The response is repeated across the temporal rows belonging to the same
epidemic because each complete exposure history is associated with one
final epidemic-level disease outcome.

## Examples

``` r
data("gamma_data", package = "EpiExposure")
head(gamma_data)
#>   epi_id block time    tmean   wetness rain        y
#> 1      1   B01    0 20.65737 12.402201    0 4.828342
#> 2      1   B01    1 21.75465  9.531243    0 4.828342
#> 3      1   B01    2 22.85553 12.114093    0 4.828342
#> 4      1   B01    3 22.68880  9.185416    0 4.828342
#> 5      1   B01    4 22.62336 11.514700    0 4.828342
#> 6      1   B01    5 23.76208 12.257759    0 4.828342
```
