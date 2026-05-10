# EpiExposure

## Overview

**EpiExposure** is an R package that provides an epidemiological framework for
structuring, summarising, and interpreting exposure–lag–response associations
fitted using **Distributed Lag Non-Linear Models (DLNMs)**.

The package is designed to work *on top of* the [`dlnm`](https://CRAN.R-project.org/package=dlnm) package,
adding domain-specific structure required in epidemiological studies,
particularly in plant disease epidemiology and climate–disease research.

While `dlnm` offers a flexible statistical infrastructure for modelling
exposure–lag–response relationships, **EpiExposure formalises epidemiological
decisions**, such as exposure windows and biologically meaningful summaries,
that are typically implemented in an ad hoc manner.

---

## Scope and Philosophy

EpiExposure **does not replace** `dlnm` and **does not redefine its statistical methodology**.
Instead, it provides:

- Explicit definition of epidemiological exposure windows
- Standardised workflows for cumulative and window-specific effects
- Clear separation between model fitting and epidemiological interpretation
- Reproducible summaries aligned with biological cycles (e.g. crop development)

In short:

> **`dlnm` answers *how to estimate*. 
> EpiExposure answers *how to interpret epidemiologically*.**

---

## Key Features

- Wrapper-based integration with `dlnm`
- Epidemiological exposure objects
- Window-based cumulative effect summaries
- Domain-aware plotting utilities
- Reproducible and comparable workflows across studies

---

## Installation

```r
# development version
remotes::install_github("yourusername/EpiExposure")

