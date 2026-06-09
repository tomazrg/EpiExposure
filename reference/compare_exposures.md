# Compare exposure-lag profiles (pairwise or multiple)

Compare exposure-lag profiles (pairwise or multiple)

## Usage

``` r
compare_exposures(
  exposures = NULL,
  exposure1 = NULL,
  exposure2 = NULL,
  lags = NULL,
  q = 0.8,
  eps = 1e-12
)
```

## Arguments

- exposures:

  Named list of exposure profiles (preferred)

- exposure1:

  Legacy single exposure

- exposure2:

  Legacy second exposure

- lags:

  Optional integer vector (same length as exposures) exp - lags = c(7,
  14, 21, 28) is equal to value 1 = lag 7, valor 2 = lag 14...

- q:

  Quantile threshold for high-exposure overlap

- eps:

  Small constant for stability

## Value

data.frame
