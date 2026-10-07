# Package index

## Data & Exposure Setup

- [`epi_data`](https://tomazrg.github.io/EpiExposure/reference/epi_data.md)
  : Simulated Beta epidemiological dataset
- [`binomial_data`](https://tomazrg.github.io/EpiExposure/reference/binomial_data.md)
  : Simulated binomial epidemiological dataset
- [`gamma_data`](https://tomazrg.github.io/EpiExposure/reference/gamma_data.md)
  : Simulated Gamma epidemiological dataset
- [`gaussian_data`](https://tomazrg.github.io/EpiExposure/reference/gaussian_data.md)
  : Simulated Gaussian epidemiological dataset
- [`nb2_data`](https://tomazrg.github.io/EpiExposure/reference/nb2_data.md)
  : Simulated negative-binomial epidemiological dataset
- [`poisson_data`](https://tomazrg.github.io/EpiExposure/reference/poisson_data.md)
  : Simulated Poisson epidemiological dataset
- [`spatial_replicate`](https://tomazrg.github.io/EpiExposure/reference/spatial_replicate.md)
  : Simulated epidemiological data with independent spatial replicates
- [`spatial_year`](https://tomazrg.github.io/EpiExposure/reference/spatial_year.md)
  : Simulated epidemiological data with independent spatial fields by
  year
- [`st_poisson`](https://tomazrg.github.io/EpiExposure/reference/st_poisson.md)
  : Simulated spatial Poisson epidemiological dataset
- [`define_exposures()`](https://tomazrg.github.io/EpiExposure/reference/define_exposures.md)
  : Define DLNM exposure templates
- [`check_identifiability()`](https://tomazrg.github.io/EpiExposure/reference/check_identifiability.md)
  : Diagnose identifiability and numerical stability of DLNM cross-basis
  designs
- [`build_design()`](https://tomazrg.github.io/EpiExposure/reference/build_design.md)
  : Build an epidemic-level DLNM design matrix
- [`prepare_response()`](https://tomazrg.github.io/EpiExposure/reference/prepare_response.md)
  : Prepare and validate the epidemic-level response for modeling

## Model Fitting & Selection

- [`fit_epidlnm()`](https://tomazrg.github.io/EpiExposure/reference/fit_epidlnm.md)
  : Fit a harmonized DLNM inferential model
- [`find_bestfit()`](https://tomazrg.github.io/EpiExposure/reference/find_bestfit.md)
  : Find the best DLNM model structure using grouped validation
- [`ensemble_bestfit()`](https://tomazrg.github.io/EpiExposure/reference/ensemble_bestfit.md)
  : Create model and lag-contribution ensembles from best-fit DLNM
  models

## Predictions & Simulations

- [`predict_outcomes()`](https://tomazrg.github.io/EpiExposure/reference/predict_outcomes.md)
  : Predict outcomes from fitted EpiExposure DLNM models
- [`compare_exposures()`](https://tomazrg.github.io/EpiExposure/reference/compare_exposures.md)
  : Compare chronological exposure profiles
- [`compare_predictions()`](https://tomazrg.github.io/EpiExposure/reference/compare_predictions.md)
  : Compare predicted outcomes between exposure scenarios
- [`simulate_exposures()`](https://tomazrg.github.io/EpiExposure/reference/simulate_exposures.md)
  : Simulate chronological exposure profiles
- [`simulate_ranges()`](https://tomazrg.github.io/EpiExposure/reference/simulate_ranges.md)
  : Generate exposure-value scenarios across DLNM lag periods
- [`simulate_scenarios()`](https://tomazrg.github.io/EpiExposure/reference/simulate_scenarios.md)
  : Simulate epidemiological exposure-history scenarios
- [`simulate_losses()`](https://tomazrg.github.io/EpiExposure/reference/simulate_losses.md)
  : Simulate yield and economic losses

## Periods & Summaries

- [`define_periods()`](https://tomazrg.github.io/EpiExposure/reference/define_periods.md)
  : Define epidemiological lag periods
- [`summarise_effects()`](https://tomazrg.github.io/EpiExposure/reference/summarise_effects.md)
  : Summarise DLNM exposure-lag effects on link and response scales
- [`reduce_effects()`](https://tomazrg.github.io/EpiExposure/reference/reduce_effects.md)
  : Reduce a fitted DLNM to a one-dimensional association
- [`compare_periods()`](https://tomazrg.github.io/EpiExposure/reference/compare_periods.md)
  : Compare period-specific DLNM effects

## Effect Decomposition & Analysis

- [`compute_eci()`](https://tomazrg.github.io/EpiExposure/reference/compute_eci.md)
  : Compute Exposure Cumulative Impact
- [`compute_ecilag()`](https://tomazrg.github.io/EpiExposure/reference/compute_ecilag.md)
  : Compute lag-specific decomposition of Exposure Cumulative Impact
  (ECI)
- [`epi_sensitivity()`](https://tomazrg.github.io/EpiExposure/reference/epi_sensitivity.md)
  : Compute sensitivity of exposure-response curves

## Plotting

- [`plot_effects()`](https://tomazrg.github.io/EpiExposure/reference/plot_effects.md)
  : Plot lag-specific or period-specific effects
- [`plot_eci()`](https://tomazrg.github.io/EpiExposure/reference/plot_eci.md)
  : Plot ECI and lag-specific contribution results
- [`plot_scenarios()`](https://tomazrg.github.io/EpiExposure/reference/plot_scenarios.md)
  : Plot predictions across epidemiological scenarios
- [`plot_losses()`](https://tomazrg.github.io/EpiExposure/reference/plot_losses.md)
  : Plot simulated yield and economic losses
- [`plot_performance()`](https://tomazrg.github.io/EpiExposure/reference/plot_performance.md)
  : Plot observed versus predicted model performance
