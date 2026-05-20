# EpiExposure

![](reference/figures/sticker.png)

`EpiExposure` is an R package designed to advance the epidemiological
interpretation of environmental drivers of plant disease epidemics.

------------------------------------------------------------------------

## Core workflow

### Analyse

Model DLNM effects using flexible and integrated statistical frameworks.

### Interpret

Summarise exposure-lag-response relationships in epidemiological terms.

### Simulate

Create exposure scenarios across lag windows for epidemiological
inference.

------------------------------------------------------------------------

## Why EpiExposure?

> `In plant disease epidemiology`, environmental effects are inherently
> dynamic, cumulative, and time-dependent. However, traditional
> analytical approaches often rely on arbitrarily defined temporal
> windows, which can fragment continuous processes and obscure critical
> epidemiological signals. `EpiExposure` addresses these limitations by
> providing tools that explicitly account for temporal continuity,
> lagged effects, and biologically interpretable summaries.

------------------------------------------------------------------------

## What does it do?

`EpiExposure` enables you to:

- Define biologically meaningful exposure windows  
- Quantify daily, cumulative, and phase-specific environmental effects  
- Identify critical time periods driving epidemics  
- Separate statistical modelling from epidemiological interpretation  
- Compare exposure–disease relationships across epidemics

------------------------------------------------------------------------

## Learning path

1.  Start with
    [`summarise_effects()`](https://tomazrg.github.io/EpiExposure/reference/summarise_effects.md)
2.  Explore
    [`predict_surface()`](https://tomazrg.github.io/EpiExposure/reference/predict_surface.md)
3.  Use
    [`simulate_scenarios()`](https://tomazrg.github.io/EpiExposure/reference/simulate_scenarios.md)
    for scenario analysis

## Installation

``` r
# development version
remotes::install_github("tomazrg/EpiExposure")
```
