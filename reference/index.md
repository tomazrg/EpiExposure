# Package index

## Data & Exposure Setup

- [`define_exposure()`](https://tomazrg.github.io/EpiExposure/reference/define_exposure.md)
  : Define DLNM exposure templates with flexible splines
- [`check_identifiability()`](https://tomazrg.github.io/EpiExposure/reference/check_identifiability.md)
  : Check DLNM cross-basis identifiability
- [`build_design_matrix()`](https://tomazrg.github.io/EpiExposure/reference/build_design_matrix.md)
  : Build epidemic-level DLNM design matrix
- [`prepare_response()`](https://tomazrg.github.io/EpiExposure/reference/prepare_response.md)
  : Prepare response variable for modeling

## Model Fitting

- [`fit_epidlnm()`](https://tomazrg.github.io/EpiExposure/reference/fit_epidlnm.md)
  : Fit DLNM inferential model

## Predictions & Simulations

- [`predict_surface()`](https://tomazrg.github.io/EpiExposure/reference/predict_surface.md)
  : Predict full DLNM exposure-lag-response surface Computes the full
  DLNM exposure–lag–response surface using the fitted model,#' without
  refitting. The surface is evaluated over a grid of exposure values and
  lags defined by the model specification or user input.
- [`predict_outcome()`](https://tomazrg.github.io/EpiExposure/reference/predict_outcome.md)
  : Predict outcome under user-defined exposure-lag profile(s)
- [`simulate_exposure()`](https://tomazrg.github.io/EpiExposure/reference/simulate_exposure.md)
  : Simulate exposure history across lags (full or patterned)
- [`simulate_range()`](https://tomazrg.github.io/EpiExposure/reference/simulate_range.md)
  : Generate DLNM simulation scenarios (scenario or profile mode)
- [`simulate_scenarios()`](https://tomazrg.github.io/EpiExposure/reference/simulate_scenarios.md)
  : Simulate epidemiological DLNM scenarios
- [`compare_exposures()`](https://tomazrg.github.io/EpiExposure/reference/compare_exposures.md)
  : Compare exposure-lag profiles (pairwise or multiple)
- [`compare_predictions()`](https://tomazrg.github.io/EpiExposure/reference/compare_predictions.md)
  : Compare predicted outcomes between multiple exposure scenarios
- [`compare_periods()`](https://tomazrg.github.io/EpiExposure/reference/compare_periods.md)
  : Compare accumulated DLNM effects between periods

## Periods & Summaries

- [`define_periods()`](https://tomazrg.github.io/EpiExposure/reference/define_periods.md)
  : Define epidemiological lag periods
- [`summarise_effects()`](https://tomazrg.github.io/EpiExposure/reference/summarise_effects.md)
  : Summarise DLNM effects
- [`recenter_effects()`](https://tomazrg.github.io/EpiExposure/reference/recenter_effects.md)
  : Recenter DLNM effects using a new reference value
- [`reduce_effects()`](https://tomazrg.github.io/EpiExposure/reference/reduce_effects.md)
  : Reduce DLNM effects to one dimension (article-consistent)

## Effect Decomposition & Analysis

- [`compute_eci()`](https://tomazrg.github.io/EpiExposure/reference/compute_eci.md)
  : Compute Exposure Cumulative Impact (ECI - the exposure profile with
  the fitted model coefficients)
- [`compute_ecilag()`](https://tomazrg.github.io/EpiExposure/reference/compute_ecilag.md)
  : Compute lag-specific decomposition of Exposure Cumulative Impact
  (ECI)
- [`identify_critical_lags()`](https://tomazrg.github.io/EpiExposure/reference/identify_critical_lags.md)
  : Identify critical lags based on daily DLNM effects
- [`lag_contribution()`](https://tomazrg.github.io/EpiExposure/reference/lag_contribution.md)
  : Quantify lag contributions to the cumulative DLNM effect
- [`df_sensitivity()`](https://tomazrg.github.io/EpiExposure/reference/df_sensitivity.md)
  : Compute sensitivity of DLNM effects (analytical or finite
  derivative)

## Datasets

- [`TargetSpot`](https://tomazrg.github.io/EpiExposure/reference/TargetSpot.md)
  : Example dataset for EpiExposure package
