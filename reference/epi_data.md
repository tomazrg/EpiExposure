# Simulated Beta epidemiological dataset

A simulated longitudinal dataset containing complete environmental
exposure histories and a continuous plant disease outcome expressed as a
proportion. The dataset is intended for demonstrating Beta regression
and the main modeling, interpretation, simulation, and prediction
workflows available in EpiExposure.

## Usage

``` r
epi_data
```

## Format

A data frame with rows representing temporal observations within
epidemic exposure histories and the following variables:

- epi_id:

  Unique numeric identifier for each epidemic.

- time:

  Chronological time index within each epidemic history.

- tmean:

  Daily mean temperature in degrees Celsius.

- rain:

  Daily rainfall in millimeters.

- wetness:

  Daily leaf wetness duration in hours.

- y:

  Continuous disease response expressed as a proportion strictly between
  0 and 1 and repeated across the rows of each epidemic history.

## Source

Simulated for the EpiExposure package.

## Details

The response is repeated across the temporal rows belonging to the same
epidemic because each complete exposure history is associated with one
final epidemic-level disease outcome.

## Examples

``` r
data("epi_data", package = "EpiExposure")
head(epi_data)
#> # A tibble: 6 × 6
#>   epi_id  time tmean  rain wetness     y
#>    <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
#> 1      1     0  22.0  0       8.83 0.330
#> 2      1     1  21.5  0      11.9  0.330
#> 3      1     2  21.5  5.38   12.5  0.330
#> 4      1     3  21.8 14.6    12.3  0.330
#> 5      1     4  20.8  0      11.4  0.330
#> 6      1     5  23.0  8.24    7.98 0.330
```
