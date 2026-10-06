# EpiExposure

`EpiExposure` is an R package for modeling and interpreting epidemiological exposure–lag–response relationships, simulating exposure histories, and predicting plant disease outcomes from environmental drivers.

---

## Core workflow

<div class="row">

<div class="col-md-4">

### Model

Construct and fit DLNM-based models using multiple statistical frameworks and response distributions.

</div>

<div class="col-md-4">

### Interpret

Summarise exposure–lag relationships across lags, epidemiological periods, and complete exposure histories.

</div>

<div class="col-md-4">

### Predict & Simulate

Evaluate observed or hypothetical exposure histories, compare scenarios, and predict expected disease outcomes.

</div>

</div>

---

## Why EpiExposure?

> In plant disease epidemiology, environmental effects are often non-linear, cumulative, and time-dependent. Traditional approaches commonly summarize weather conditions within predefined temporal windows, which may obscure delayed or continuously varying associations. `EpiExposure` provides a unified DLNM-based workflow for preserving exposure histories, modeling lagged associations, and translating fitted models into epidemiologically interpretable summaries, simulations, and predictions.

---

## What does it do?

`EpiExposure` enables you to:

- Model non-linear exposure–lag–response relationships
- Quantify lag- and period-specific environmental associations
- Define biologically meaningful epidemiological periods
- Summarise complete exposure histories using the Exposure Cumulative Impact (ECI)
- Decompose cumulative exposure impacts into exact lag-specific contributions
- Simulate and compare hypothetical exposure–lag scenarios
- Predict expected disease outcomes from complete exposure histories
- Compare and ensemble alternative DLNM specifications using cross-validation

---

## Learning path

1. Define exposure–lag structures with `define_exposures()`.
2. Build the DLNM design with `build_design()`.
3. Fit the epidemiological model with `fit_epidlnm()`.
4. Interpret exposure–lag relationships using `summarise_effects()`.
5. Explore hypothetical exposure histories with `simulate_exposures()`.
6. Predict expected outcomes using `predict_outcomes()`.
7. Evaluate alternative model structures with `find_bestfit()` and `ensemble_bestfit()`.

---

## Installation

`EpiExposure` is currently available from GitHub. Install the package with:

```r
install.packages("pak")
pak::pkg_install("tomazrg/EpiExposure")
```

Then load the package:

```r
library(EpiExposure)
```
