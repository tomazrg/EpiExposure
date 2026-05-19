<img src="reference/figures/sticker.png"
     align="right"
     width="380px"
     style="margin-left: 20px; margin-right: -100px; margin-top: 10px;"/>


## EpiExposure

`EpiExposure` is an R package designed to advance the epidemiological interpretation of environmental drivers of plant disease epidemics. It provides a structured and reproducible framework for analysing exposure–lag–response relationships, enabling researchers to move beyond descriptive associations toward biologically meaningful inference.

In plant disease epidemiology, environmental effects are inherently dynamic, cumulative, and time-dependent. However, traditional analytical approaches often rely on arbitrarily defined temporal windows, which can fragment continuous processes and obscure critical epidemiological signals.

EpiExposure addresses these limitations by providing tools that explicitly account for temporal continuity, lagged effects, and biologically interpretable summaries. By integrating epidemiological reasoning directly into the analytical workflow, the package allows users to identify critical exposure periods, quantify cumulative effects, and interpret environmental drivers within the context of disease development.

---

## Scope and Scientific Contribution

EpiExposure introduces an epidemiological layer for modelling environmental effects on plant disease dynamics. Rather than focusing solely on statistical estimation, the package structures the entire analytical pipeline around epidemiologically meaningful concepts.

Specifically, EpiExposure enables:

- Formal definition of exposure windows aligned with biological processes . 
- Quantification of cumulative and phase-specific effects across epidemics.
- Consistent separation between model fitting and epidemiological interpretation.  
- Identification of critical time periods driving disease development.  
- Reproducible workflows for comparative and multi-epidemic analyses.

By embedding these components into a unified framework, the package transforms how exposure–disease relationships are analysed, interpreted, and communicated.

---

## Conceptual Framework

Plant diseases arise from complex and dynamic interactions between host, pathogen, and environment. Environmental drivers do not act instantaneously but accumulate and interact over time, influencing infection cycles, latent periods, and epidemic progression.

EpiExposure operationalizes this perspective by treating environmental exposure as a continuous and structured process. Instead of relying on discrete summaries, it enables the estimation and interpretation of time-resolved exposure effects, capturing both immediate and lagged influences.

This approach allows researchers to:

- Detect non-linear and delayed environmental effects  
- Distinguish between short-term and cumulative drivers  
- Link climatic variability to specific epidemiological phases  
- Improve mechanistic understanding of disease–environment interactions  

---

## Key Features

- Epidemiology-driven analytical framework for exposure–lag modelling  
- Structured definition of biologically meaningful exposure windows  
- Integrated tools for daily and cumulative effect estimation  
- Scenario simulation for exploring environmental impacts on epidemics  
- Compatibility with multiple statistical engines (frequentist and Bayesian)  
- Reproducible workflows for research, teaching, and decision support  

---

## Why EpiExposure?

EpiExposure is not simply an implementation tool — it is a framework for translating complex statistical outputs into epidemiological knowledge.

In short:

> EpiExposure enables researchers to move from statistical estimation to epidemiological understanding of environmental drivers of disease.


## Core workflow

<div class="row">

<div class="col-md-4">

### Simulate

Create exposure scenarios across lag windows for ecological inference.

</div>

<div class="col-md-4">

### Analyse

Model DLNM effects using flexible statistical frameworks.

</div>

<div class="col-md-4">

### Interpret

Summarise exposure-lag-response relationships in epidemiological terms.

</div>

</div>

## Learning path

1. Start with `summarise_effects()`
2. Explore `predict_surface()`
3. Use `simulate_scenarios()` for scenario analysis

## Installation

```r
# development version
remotes::install_github("tomazrg/EpiExposure")
```
