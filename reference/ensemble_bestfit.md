# Create model and lag-contribution ensembles from best-fit DLNM models

Combines out-of-fold (OOF) predictions produced by the current
\`find_bestfit()\` into one or more model-level ensembles and,
optionally, combines lag-specific contributions from retained full-data
model refits.

## Usage

``` r
ensemble_bestfit(
  bestfit,
  predictions = NULL,
  lag_data = NULL,
  fit_list = NULL,
  data = NULL,
  group = "epi_id",
  var = NULL,
  compute_ecilag_args = list(),
  ensemble_scope = c("model", "lag", "both"),
  method = "unweighted",
  top_n = 3,
  model_id = NULL,
  weight_metric = NULL,
  threshold = NULL,
  weight_transform = c("softmax", "positive", "rank_inverse", "uniform"),
  stacking_model = c("ridge", "lm"),
  stack_objective = c("regularized", "min_error", "constrained"),
  lambda = 1e-06,
  stack_intercept = TRUE,
  model_col = "model_id",
  id_cols = c("group", "fold"),
  lag_group_cols = c("var", "lag"),
  lg_strategy = c("requested_available", "model", "strict"),
  rl_weights = TRUE,
  verbose = TRUE
)
```

## Arguments

- bestfit:

  Non-empty data frame returned by the current \`find_bestfit()\`. The
  function requires the standardized metadata stored by that function,
  including \`family\`, \`outcome_type\`, \`rank_metric\`,
  \`threshold\`, the prediction contract, the cross-validation method,
  \`fold_assignments\`, the common fitted \`max_lag\`, the expected
  profile length (\`max_lag + 1\`), and the EpiExposure exact-profile
  contract.

- predictions:

  Optional data frame of OOF predictions. If \`NULL\`, \`attr(bestfit,
  "predictions")\` is used. Under the current \`find_bestfit()\`
  contract, required columns include \`model_col\`, \`group\`, \`fold\`,
  \`observed\`, and \`predicted\`, in addition to any columns listed in
  \`id_cols\`.

  All selected models must contain predictions for the same OOF
  identifiers. The function deliberately does not silently reduce the
  ensemble to an intersection of different validation supports, because
  model weights and ensemble metrics would otherwise be evaluated on a
  different sample from the individual-model ranking.

- lag_data:

  Optional data frame of deterministic lag-specific ECI decompositions.
  Required columns are \`model_col\`, all \`lag_group_cols\`, and
  \`ECI_weighted\`. If \`NULL\` and a lag ensemble is requested,
  contributions are generated with \`compute_ecilag()\` from retained
  full-data fits.

- fit_list:

  Optional named list of fitted models used for automatic lag
  decomposition. If \`NULL\`, \`attr(bestfit, "fits")\` is used.

- data:

  Optional long-format exposure data used only for automatic lag
  decomposition, while ignored/not required for ensemble_scope =
  "model". It must contain \`group\`, a column named \`time\`, and the
  exposures requested from the retained fitted models. Rows are ordered
  internally by \`group\` and \`time\`, and time must satisfy the
  regular-series requirements used elsewhere in EpiExposure. Every group
  must contain exactly the common fitted \`max_lag + 1\` observations.
  Longer profiles are not truncated and shorter profiles are not padded.

- group:

  Character scalar naming the grouping column in \`data\` for automatic
  lag decomposition.

- var:

  Optional unique character vector of exposure variables requested for
  automatic lag decomposition.

- compute_ecilag_args:

  Named list of additional arguments passed to \`compute_ecilag()\`. It
  cannot override \`data\`, \`group\`, \`fit\`, \`var\`, or \`profile\`.
  \`uncertainty = TRUE\` is not supported by \`ensemble_bestfit()\` v1
  for lag ensembles because parameter draws from separately fitted
  models do not form a defined joint draw distribution.

- ensemble_scope:

  Character. One of \`"model"\`, \`"lag"\`, or \`"both"\`.

- method:

  One or more methods selected from \`"best"\`, \`"hard_voting"\`,
  \`"unweighted"\`, \`"weighted"\`, and \`"stacked"\`. Multiple methods
  may be evaluated in one call, for example \`c("unweighted",
  "weighted", "stacked")\`.

  \`"stacked"\` and \`"hard_voting"\` are model-level methods only.
  \`"hard_voting"\` additionally requires a binomial outcome.

- top_n:

  Positive integer number of top-ranked models selected when \`model_id
  = NULL\`.

- model_id:

  Optional vector of unique model identifiers. When supplied, \`top_n\`
  is ignored. Every requested model must have a finite value for
  \`weight_metric\`.

- weight_metric:

  Optional performance metric used to rank/select models. When
  \`model_id = NULL\`, the \`top_n\` models are selected according to
  this metric. For \`method = "weighted"\`, it is also used to construct
  model weights. For \`method = "stacked"\`, it determines which models
  are selected but is not transformed directly into stacking weights;
  those weights are estimated from the aligned OOF predictions by the
  level-2 model. If \`NULL\`, the function inherits \`attr(bestfit,
  "rank_metric")\`. Metric direction is determined by
  \`.metric_direction()\`.

- threshold:

  Optional probability threshold strictly between 0 and 1 for binomial
  classification metrics and \`hard_voting\`. If \`NULL\`,
  \`attr(bestfit, "threshold")\` is used.

  When \`weight_metric\` is threshold-dependent (\`Accuracy\`,
  \`Balanced_Accuracy\`, \`Sensitivity\`, \`Specificity\`, \`F1\`,
  \`MCC\`, or \`Precision\`), a user-supplied threshold must equal the
  threshold used by \`find_bestfit()\`. Otherwise the stored
  model-ranking metric and the ensemble threshold would refer to
  different classification rules, so the function stops explicitly.

- weight_transform:

  Character. Transformation used only by \`method = "weighted"\`:

  \`"softmax"\`

  : Stable softmax of the direction-adjusted metric.

  \`"positive"\`

  : For maximize metrics, negative values are truncated to zero and
    positive values are normalized. For minimize metrics, reciprocal
    loss is normalized. If no positive utility exists, the function
    stops instead of silently switching to uniform weights.

  \`"rank_inverse"\`

  : Weights proportional to inverse selected-model rank.

  \`"uniform"\`

  : Equal weights.

- stacking_model:

  Character. \`"ridge"\` or \`"lm"\`. It must be coherent with
  \`stack_objective\`: \`"regularized"\` requires \`"ridge"\` and
  \`"min_error"\` requires \`"lm"\`. With \`"constrained"\`, \`"ridge"\`
  adds the penalty \`lambda\` whereas \`"lm"\` uses no ridge penalty.

- stack_objective:

  Character. \`"regularized"\`, \`"min_error"\`, or \`"constrained"\`.

  \`"regularized"\` fits ridge stacking. \`"min_error"\` fits ordinary
  least-squares stacking. \`"constrained"\` minimizes squared error
  subject to non-negative model weights summing to one, using
  projected-gradient optimization on the simplex.

- lambda:

  Non-negative ridge penalty used by regularized stacking and by
  constrained stacking when \`stacking_model = "ridge"\`.

- stack_intercept:

  Logical. Include an intercept for unconstrained stacking.
  \`stack_objective = "constrained"\` requires \`stack_intercept =
  FALSE\`, because an unconstrained intercept would destroy the
  convex-combination interpretation.

- model_col:

  Character scalar naming the model identifier column.

- id_cols:

  Unique character vector identifying an OOF prediction row. Under the
  current \`find_bestfit()\` output the default \`c("group", "fold")\`
  is appropriate for both LOOCV and grouped k-fold.

- lag_group_cols:

  Unique character vector defining one lag-ensemble unit. It must
  include \`"lag"\`. Include the epidemic/group column as well when lag
  contributions are group-specific, for example \`c("epi_id", "var",
  "lag")\`.

- lg_strategy:

  Character. Controls automatic lag-variable availability:
  \`"requested_available"\` uses requested variables available in each
  model, \`"model"\` uses every exposure fitted in each model, and
  \`"strict"\` requires every requested variable in every selected
  model. For lag aggregation, \`"strict"\` additionally requires every
  lag unit to be available from all selected models.

- rl_weights:

  Logical. For lag ensembles, renormalize model weights among models
  actually available within each lag unit. If \`FALSE\`, missing models
  retain their absent weight mass and the lag contribution can therefore
  shrink toward zero.

- verbose:

  Logical. Print concise information.

## Value

A list containing:

- \`ensemble_summary\`:

  Model-level performance, one row per requested method.

- \`ensemble_predictions\`:

  OOF ensemble predictions with an explicit \`method\` column and the
  aligned individual-model predictions.

- \`ensemble_by_lag\`:

  Ensemble lag-specific ECI decomposition.

- \`lag_data\`:

  Lag data used for aggregation.

- \`model_weights\`:

  Method-specific model weights. For stacking, these are the
  \*\*final\*\* level-2 weights fitted to all available OOF rows after
  cross-fitted ensemble performance has been evaluated.

- \`selected_models\`:

  Selected rows from \`bestfit\`.

- \`stacking_coefficients_by_fold\`:

  For stacked methods, level-2 coefficients fitted without the current
  OOF fold and used to predict that fold.

- \`stacking_final_coefficients\`:

  Final stacking coefficients fitted to the complete OOF matrix for
  downstream use with full-data base-model refits.

- \`metric_warnings\`:

  Warnings captured while evaluating ensemble performance metrics.

The returned list also carries the inherited family, threshold,
population/expected-response prediction contract, CV metadata
(\`cv_method\`, \`cv_scheme\`, \`k\`, \`cv_stratified\`,
\`fold_assignments\`, and \`fold_balance\`), and temporal metadata
(\`max_lag\`, \`history_length\`, \`history_contract\`, and
\`time_step\`).

## Details

The function inherits the EpiExposure v1 prediction contract from
\`find_bestfit()\`: model-level ensemble predictions represent
population/fixed-component \*\*expected responses\*\*, with fitted
random effects set to zero. The OOF predictions used here are
deterministic predictions from the harmonized central parameter
estimate; coefficient/posterior uncertainty, residual noise,
future-observation noise, and group-specific random effects are not
added.

\## Model-level OOF ensemble contract

The current \`find_bestfit()\` produces one deterministic OOF
expected-response prediction per successfully validated group and
candidate model. This function combines those predictions; it does not
call the original model prediction methods again.

Selected models must share the same OOF support. This strict alignment
is intentional. Silently using only complete rows across models could
make an ensemble appear better simply because difficult groups were
removed.

For non-binary outcomes, ensemble performance uses CCC, Cb, Pearson
correlation, RMSE, and MAE. For binomial outcomes, probability-based
ensembles retain probabilities for ROC AUC, Brier score, and Log Loss
and apply \`threshold\` only for classification metrics.

\`hard_voting\` is different: each individual predicted probability is
thresholded first, and the final class is the majority vote. Exact ties
are resolved with the class produced by the highest-ranked selected
model. \`vote_fraction\` is reported for interpretation but is not
treated as a calibrated probability, so ROC AUC, Brier score, and Log
Loss are returned as \`NA\` for hard voting.

\## Stacking and the CV folds

Stacking has two separate outputs that should not be confused.

Ensemble performance is calculated from \*\*fold-cross-fitted level-2
predictions\*\*. For OOF fold \\f\\, stacking coefficients are fitted
using the stored OOF rows from folds other than \\f\\, and those
coefficients are then applied to fold \\f\\. With LOOCV this reduces to
the historical leave-one-row-out level-2 calculation. With grouped
k-fold, every epidemic in the same validation fold is held out from
level-2 fitting together.

In \`method = "stacked"\`, \`weight_metric\` is used to rank and select
the base models when \`model_id = NULL\`; it is not transformed directly
into stacking weights. Stacking weights are level-2 regression
coefficients estimated from the aligned OOF predictions of the selected
models. Therefore, a stacking weight represents the contribution of one
model conditional on the other selected models and does not necessarily
follow the individual-model ranking.

During cross-fitted performance evaluation, these coefficients are
re-estimated separately using the OOF rows outside each validation fold.
After this evaluation, a final stacking model is fitted to the complete
OOF prediction matrix. The resulting coefficients are reported in
\`model_weights\` and \`stacking_final_coefficients\` and are intended
for combining the selected base models after they have been refitted to
all available data.

For unconstrained ridge or least-squares stacking, these coefficients
may be negative and are not required to sum to one. Only
\`stack_objective = "constrained"\` enforces non-negative weights that
sum to one, with no stacking intercept.

This is a level-2 cross-fitting procedure based on the OOF prediction
matrix already produced by \`find_bestfit()\`. It is not a fully nested
re-fitting of every base learner inside an additional outer validation
loop; therefore ensemble performance should be interpreted as OOF
ensemble performance for model development rather than as an independent
external-validation estimate.

\## Exact common lag/profile contract

\`ensemble_bestfit()\` inherits the exact-profile contract established
by the current \`find_bestfit()\`: all candidate/fitted exposures use
one common non-negative integer \`max_lag\`, and every original exposure
profile contains exactly \`max_lag + 1\` equally spaced observations.

For automatic lag decomposition, the retained fitted models must report
the same \`max_lag\`, \`history_length\`, and \`history_contract\` as
\`bestfit\`, and every group supplied in \`data\` must contain exactly
that profile length. Profiles are never truncated, padded, or silently
realigned.

If precomputed \`lag_data\` is supplied, the original long-format
exposure profiles are no longer available to this function. In that
route, \`ensemble_bestfit()\` validates the lag values against the
inherited common \`max_lag\` but cannot reconstruct the original profile
row counts.

\## Lag ensembles

Lag ensembles combine lag-specific ECI contributions, not raw DLNM
coefficients:

\$\$ECI\_{ens,l} = \sum_m \alpha_m ECI\_{m,l}.\$\$

The model weights \\\alpha_m\\ come from OOF model performance, while
lag-specific ECI values are obtained from \`compute_ecilag()\` applied
to the selected models refitted on the complete dataset.

Percent contribution is based on the absolute contribution magnitude
within each non-lag grouping unit:

\$\$100 \|ECI\_{ens,l}\| / \sum_l \|ECI\_{ens,l}\|.\$\$

Ensemble uncertainty for lag contributions is deliberately not produced
in EpiExposure v1. Independently sampled coefficient/posterior draws
from separately fitted models do not automatically form a valid joint
draw distribution, so arbitrary draw-index pairing is not performed.
