# Simulated epidemiological data with independent spatial replicates

A simulated longitudinal dataset containing complete environmental
exposure profiles for fixed field locations observed across independent
epidemiological replicates. The dataset is intended for demonstrating
distributed exposure-lag models with separate spatial fields across
epidemic replicates.

## Usage

``` r
spatial_replicate
```

## Format

A data frame with rows representing temporal observations within
epidemic exposure profiles and the following variables:

- epi_id:

  Unique identifier combining field location and epidemiological
  replicate.

- location_id:

  Identifier for the field location.

- epidemic_replicate:

  Identifier for an independent epidemiological replicate.

- block:

  Experimental block identifier within a replicate.

- block_id:

  Identifier combining epidemiological replicate and experimental block.

- time:

  Chronological time index within each epidemic profile.

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
  epidemic profile.

## Source

Simulated for the EpiExposure package.

## Details

The response and location-level descriptors are repeated across the
temporal rows belonging to the same epidemic because each complete
exposure profile is associated with one final epidemic-level disease
outcome.

## Examples

``` r
data("spatial_replicate", package = "EpiExposure")
head(spatial_replicate)
#>             epi_id location_id epidemic_replicate block        block_id time
#> 1 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    0
#> 2 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    1
#> 3 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    2
#> 4 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    3
#> 5 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    4
#> 6 L001_Replicate_1        L001        Replicate_1   B01 Replicate_1_B01    5
#>     x_coord    y_coord    tmean   wetness rain y
#> 1 0.7586476 -0.4480224 23.89557  6.798144    0 5
#> 2 0.7586476 -0.4480224 24.45162  7.355172    0 5
#> 3 0.7586476 -0.4480224 24.49880  8.565241    0 5
#> 4 0.7586476 -0.4480224 24.68348 12.894026    0 5
#> 5 0.7586476 -0.4480224 25.65530 14.958210    0 5
#> 6 0.7586476 -0.4480224 23.98108 13.466749    0 5
```
