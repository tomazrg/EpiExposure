# Simulated spatial Poisson epidemiological dataset

A simulated longitudinal dataset containing complete environmental
exposure histories and spatial coordinates for demonstrating distributed
exposure-lag models with Matérn spatial dependence and a Poisson disease
outcome.

## Usage

``` r
st_poisson
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

- x_coord:

  First numeric spatial coordinate of the epidemic location.

- y_coord:

  Second numeric spatial coordinate of the epidemic location.

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

The response, block, and spatial coordinates are repeated across the
temporal rows belonging to the same epidemic because each complete
exposure history is associated with one final epidemic-level disease
outcome.

## Examples

``` r
data("st_poisson", package = "EpiExposure")
head(st_poisson)
#>   epi_id block time    x_coord    y_coord    tmean  wetness     rain y
#> 1      1   B01    0 -0.4210796 -0.1649767 24.06189 10.81782 0.000000 4
#> 2      1   B01    1 -0.4210796 -0.1649767 24.82615 11.15336 0.000000 4
#> 3      1   B01    2 -0.4210796 -0.1649767 26.11794 12.20862 0.000000 4
#> 4      1   B01    3 -0.4210796 -0.1649767 26.17299 11.70761 6.440981 4
#> 5      1   B01    4 -0.4210796 -0.1649767 26.39594 13.08992 0.000000 4
#> 6      1   B01    5 -0.4210796 -0.1649767 26.23450 14.78350 0.000000 4
```
