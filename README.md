<img src="reference/figures/sticker.png"
     align="right"
     width="380px"
     style="margin-left: 20px; margin-right: -100px; margin-top: 10px;"/>

# EpiExposure

`EpiExposure` is an R package for modeling and interpreting epidemiological exposure–lag relationships between environmental drivers and plant disease outcomes.

---

## Core workflow

<div class="row">

<div class="col-md-4">

### Analyse

Model DLNM effects using flexible and integrated statistical frameworks.

</div>

<div class="col-md-4">

### Interpret

Summarise exposure-lag-response relationships in epidemiological terms.

</div>

<div class="col-md-4">

### Simulate

Create exposure scenarios across lag windows for epidemiological inference.

</div>

</div>

---

## Why EpiExposure?

> `In plant disease epidemiology`, environmental effects are inherently dynamic, cumulative, and time-dependent. However, traditional analytical approaches often rely on arbitrarily defined temporal windows, which can fragment continuous processes and obscure critical epidemiological signals. `EpiExposure` addresses these limitations by providing tools that explicitly account for temporal continuity, lagged effects, and biologically interpretable summaries.

---

## What does it do?

`EpiExposure` enables you to:

- Quantify daily, cumulative, and phase-specific environmental effects  
- Define biologically meaningful exposure periods 
- Identify critical time periods driving diseases  
- Separate statistical modeling from epidemiological interpretation  
- Compare exposure–disease relationships across epidemiological scenarios  

---

## Learning path

1. Start with `summarise_effects()` to interpret the exposure-lag-effects.
2. Explore `predict_surface()` to visualize the extent of the effects of 3D exposure delay.
3. Use `simulate_scenarios()` to simulate and analyze epidemiological issues.
4. Use `predict_outcome()` to predict outcomes.

## How to install the `EpiExposure`?

```r
# development version
remotes::install_github("tomazrg/EpiExposure")
```
