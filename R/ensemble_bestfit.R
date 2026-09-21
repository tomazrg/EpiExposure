#' Create model and lag-contribution ensembles from best-fit DLNM models
#'
#' Combines out-of-fold (OOF) predictions produced by the current
#' `find_bestfit()` into one or more model-level ensembles and, optionally,
#' combines lag-specific contributions from retained full-data model refits.
#'
#' The function inherits the EpiExposure v1 prediction contract from
#' `find_bestfit()`: model-level ensemble predictions represent
#' population/fixed-component **expected responses**, with fitted random effects
#' set to zero. The OOF predictions used here are deterministic predictions from
#' the harmonized central parameter estimate; coefficient/posterior uncertainty,
#' residual noise, future-observation noise, and group-specific random effects
#' are not added.
#'
#' @param bestfit Non-empty data frame returned by the current
#'   `find_bestfit()`. The function requires the standardized metadata stored by
#'   that function, including `family`, `outcome_type`, `rank_metric`,
#'   `threshold`, the prediction contract, the cross-validation method,
#'   `fold_assignments`, the common fitted `max_lag`, the expected history
#'   length (`max_lag + 1`), and the EpiExposure exact-history contract.
#' @param predictions Optional data frame of OOF predictions. If `NULL`,
#'   `attr(bestfit, "predictions")` is used. Under the current
#'   `find_bestfit()` contract, required columns include `model_col`, `group`,
#'   `fold`, `observed`, and `predicted`, in addition to any columns listed in
#'   `id_cols`.
#'
#'   All selected models must contain predictions for the same OOF identifiers.
#'   The function deliberately does not silently reduce the ensemble to an
#'   intersection of different validation supports, because model weights and
#'   ensemble metrics would otherwise be evaluated on a different sample from
#'   the individual-model ranking.
#' @param lag_data Optional data frame of deterministic lag-specific
#'   ECI decompositions. Required columns are `model_col`, all `lag_group_cols`, and
#'   `ECI_weighted`. If `NULL` and a lag ensemble is requested, contributions
#'   are generated with `compute_ecilag()` from retained full-data fits.
#' @param fit_list Optional named list of fitted models used for automatic lag
#'   decomposition. If `NULL`, `attr(bestfit, "fits")` is used.
#' @param data Optional long-format exposure data used only for automatic lag
#'   decomposition, while ignored/not required for ensemble_scope = "model". It must contain `group`, a column named `time`, and the
#'   exposures requested from the retained fitted models. Rows are ordered
#'   internally by `group` and `time`, and time must satisfy the regular-series
#'   requirements used elsewhere in EpiExposure. Every group must contain
#'   exactly the common fitted `max_lag + 1` observations. Longer histories are
#'   not truncated and shorter histories are not padded.
#' @param group Character scalar naming the grouping column in `data` for
#'   automatic lag decomposition.
#' @param var Optional unique character vector of exposure variables requested
#'   for automatic lag decomposition.
#' @param compute_ecilag_args Named list of additional arguments passed to
#'   `compute_ecilag()`. It cannot override `data`, `group`, `fit`, `var`, or
#'   `profile`. `uncertainty = TRUE` is not supported by
#'   `ensemble_bestfit()` v1 for lag ensembles because parameter draws from
#'   separately fitted models do not form a defined joint draw distribution.
#' @param ensemble_scope Character. One of `"model"`, `"lag"`, or `"both"`.
#' @param method One or more methods selected from `"best"`, `"hard_voting"`,
#'   `"unweighted"`, `"weighted"`, and `"stacked"`. Multiple methods may be
#'   evaluated in one call, for example
#'   `c("unweighted", "weighted", "stacked")`.
#'
#'   `"stacked"` and `"hard_voting"` are model-level methods only.
#'   `"hard_voting"` additionally requires a binomial outcome.
#' @param top_n Positive integer number of top-ranked models selected when
#'   `model_id = NULL`.
#' @param model_id Optional vector of unique model identifiers. When supplied,
#'   `top_n` is ignored. Every requested model must have a finite value for
#'   `weight_metric`.
#' @param weight_metric Optional performance metric used to rank/select models.
#'   When `model_id = NULL`, the `top_n` models are selected according to this
#'   metric. For `method = "weighted"`, it is also used to construct model
#'   weights. For `method = "stacked"`, it determines which models are selected
#'   but is not transformed directly into stacking weights; those weights are
#'   estimated from the aligned OOF predictions by the level-2 model. If `NULL`,
#'   the function inherits `attr(bestfit, "rank_metric")`. Metric direction is
#'   determined by `.metric_direction()`.
#' @param threshold Optional probability threshold strictly between 0 and 1 for
#'   binomial classification metrics and `hard_voting`. If `NULL`,
#'   `attr(bestfit, "threshold")` is used.
#'
#'   When `weight_metric` is threshold-dependent (`Accuracy`,
#'   `Balanced_Accuracy`, `Sensitivity`, `Specificity`, `F1`, `MCC`, or
#'   `Precision`), a user-supplied threshold must equal the threshold used by
#'   `find_bestfit()`. Otherwise the stored model-ranking metric and the
#'   ensemble threshold would refer to different classification rules, so the
#'   function stops explicitly.
#' @param weight_transform Character. Transformation used only by
#'   `method = "weighted"`:
#'
#'   \describe{
#'     \item{`"softmax"`}{Stable softmax of the direction-adjusted metric.}
#'     \item{`"positive"`}{For maximize metrics, negative values are truncated
#'       to zero and positive values are normalized. For minimize metrics,
#'       reciprocal loss is normalized. If no positive utility exists, the
#'       function stops instead of silently switching to uniform weights.}
#'     \item{`"rank_inverse"`}{Weights proportional to inverse selected-model
#'       rank.}
#'     \item{`"uniform"`}{Equal weights.}
#'   }
#' @param stacking_model Character. `"ridge"` or `"lm"`. It must be coherent
#'   with `stack_objective`: `"regularized"` requires `"ridge"` and
#'   `"min_error"` requires `"lm"`. With `"constrained"`, `"ridge"` adds the
#'   penalty `lambda` whereas `"lm"` uses no ridge penalty.
#' @param stack_objective Character. `"regularized"`, `"min_error"`, or
#'   `"constrained"`.
#'
#'   `"regularized"` fits ridge stacking. `"min_error"` fits ordinary
#'   least-squares stacking. `"constrained"` minimizes squared error subject to
#'   non-negative model weights summing to one, using projected-gradient
#'   optimization on the simplex.
#' @param lambda Non-negative ridge penalty used by regularized stacking and by
#'   constrained stacking when `stacking_model = "ridge"`.
#' @param stack_intercept Logical. Include an intercept for unconstrained
#'   stacking. `stack_objective = "constrained"` requires
#'   `stack_intercept = FALSE`, because an unconstrained intercept would destroy
#'   the convex-combination interpretation.
#' @param model_col Character scalar naming the model identifier column.
#' @param id_cols Unique character vector identifying an OOF prediction row.
#'   Under the current `find_bestfit()` output the default
#'   `c("group", "fold")` is appropriate for both LOOCV and grouped k-fold.
#' @param lag_group_cols Unique character vector defining one lag-ensemble unit.
#'   It must include `"lag"`. Include the epidemic/group column as well when lag
#'   contributions are group-specific, for example
#'   `c("epi_id", "var", "lag")`.
#' @param lg_strategy Character. Controls automatic lag-variable availability:
#'   `"requested_available"` uses requested variables available in each model,
#'   `"model"` uses every exposure fitted in each model, and `"strict"` requires
#'   every requested variable in every selected model. For lag aggregation,
#'   `"strict"` additionally requires every lag unit to be available from all
#'   selected models.
#' @param rl_weights Logical. For lag ensembles, renormalize model weights among
#'   models actually available within each lag unit. If `FALSE`, missing models
#'   retain their absent weight mass and the lag contribution can therefore
#'   shrink toward zero.
#' @param verbose Logical. Print concise information.
#'
#' @return A list containing:
#'
#'   \describe{
#'     \item{`ensemble_summary`}{Model-level performance, one row per requested
#'       method.}
#'     \item{`ensemble_predictions`}{OOF ensemble predictions with an explicit
#'       `method` column and the aligned individual-model predictions.}
#'     \item{`ensemble_by_lag`}{Ensemble lag-specific ECI decomposition.}
#'     \item{`lag_data`}{Lag data used for aggregation.}
#'     \item{`model_weights`}{Method-specific model weights. For stacking,
#'       these are the **final** level-2 weights fitted to all available OOF
#'       rows after cross-fitted ensemble performance has been evaluated.}
#'     \item{`selected_models`}{Selected rows from `bestfit`.}
#'     \item{`stacking_coefficients_by_fold`}{For stacked methods, level-2
#'       coefficients fitted without the current OOF fold and used to predict
#'       that fold.}
#'     \item{`stacking_final_coefficients`}{Final stacking coefficients fitted
#'       to the complete OOF matrix for downstream use with full-data base-model
#'       refits.}
#'     \item{`metric_warnings`}{Warnings captured while evaluating ensemble
#'       performance metrics.}
#'   }
#'
#'   The returned list also carries the inherited family, threshold,
#'   population/expected-response prediction contract, CV metadata
#'   (`cv_method`, `cv_scheme`, `k`, `cv_stratified`, `fold_assignments`, and
#'   `fold_balance`), and temporal metadata (`max_lag`, `history_length`,
#'   `history_contract`, and `time_step`).
#'
#' @details
#' ## Model-level OOF ensemble contract
#'
#' The current `find_bestfit()` produces one deterministic OOF expected-response
#' prediction per successfully validated group and candidate model. This
#' function combines those predictions; it does not call the original model
#' prediction methods again.
#'
#' Selected models must share the same OOF support. This strict alignment is
#' intentional. Silently using only complete rows across models could make an
#' ensemble appear better simply because difficult groups were removed.
#'
#' For non-binary outcomes, ensemble performance uses CCC, Cb, Pearson
#' correlation, RMSE, and MAE. For binomial outcomes, probability-based
#' ensembles retain probabilities for ROC AUC, Brier score, and Log Loss and
#' apply `threshold` only for classification metrics.
#'
#' `hard_voting` is different: each individual predicted probability is
#' thresholded first, and the final class is the majority vote. Exact ties are
#' resolved with the class produced by the highest-ranked selected model.
#' `vote_fraction` is reported for interpretation but is not treated as a
#' calibrated probability, so ROC AUC, Brier score, and Log Loss are returned
#' as `NA` for hard voting.
#'
#' ## Stacking and the CV folds
#'
#' Stacking has two separate outputs that should not be confused.
#'
#' Ensemble performance is calculated from **fold-cross-fitted level-2
#' predictions**. For OOF fold \eqn{f}, stacking coefficients are fitted using
#' the stored OOF rows from folds other than \eqn{f}, and those coefficients are
#' then applied to fold \eqn{f}. With LOOCV this reduces to the historical
#' leave-one-row-out level-2 calculation. With grouped k-fold, every epidemic
#' in the same validation fold is held out from level-2 fitting together.
#'
#' In `method = "stacked"`, `weight_metric` is used to rank and select the base
#' models when `model_id = NULL`; it is not transformed directly into stacking
#' weights. Stacking weights are level-2 regression coefficients estimated from
#' the aligned OOF predictions of the selected models. Therefore, a stacking
#' weight represents the contribution of one model conditional on the other
#' selected models and does not necessarily follow the individual-model ranking.
#'
#' During cross-fitted performance evaluation, these coefficients are
#' re-estimated separately using the OOF rows outside each validation fold.
#' After this evaluation, a final stacking model is fitted to the complete OOF
#' prediction matrix. The resulting coefficients are reported in `model_weights`
#' and `stacking_final_coefficients` and are intended for combining the selected
#' base models after they have been refitted to all available data.
#'
#' For unconstrained ridge or least-squares stacking, these coefficients may be
#' negative and are not required to sum to one. Only
#' `stack_objective = "constrained"` enforces non-negative weights that sum to
#' one, with no stacking intercept.
#'
#' This is a level-2 cross-fitting procedure based on the OOF prediction matrix
#' already produced by `find_bestfit()`. It is not a fully nested re-fitting of
#' every base learner inside an additional outer validation loop; therefore
#' ensemble performance should be interpreted as OOF ensemble performance for
#' model development rather than as an independent external-validation estimate.
#'
#' ## Exact common lag/history contract
#'
#' `ensemble_bestfit()` inherits the exact-history contract established by the
#' current `find_bestfit()`: all candidate/fitted exposures use one common
#' non-negative integer `max_lag`, and every original exposure history contains
#' exactly `max_lag + 1` equally spaced observations.
#'
#' For automatic lag decomposition, the retained fitted models must report the
#' same `max_lag`, `history_length`, and `history_contract` as `bestfit`, and
#' every group supplied in `data` must contain exactly that history length.
#' Histories are never truncated, padded, or silently realigned.
#'
#' If precomputed `lag_data` is supplied, the original long-format exposure
#' histories are no longer available to this function. In that route,
#' `ensemble_bestfit()` validates the lag values against the inherited common
#' `max_lag` but cannot reconstruct the original history row counts.
#'
#' ## Lag ensembles
#'
#' Lag ensembles combine lag-specific ECI contributions, not raw DLNM
#' coefficients:
#'
#' \deqn{ECI_{ens,l} = \sum_m \alpha_m ECI_{m,l}.}
#'
#' The model weights \eqn{\alpha_m} come from OOF model performance, while
#' lag-specific ECI values are obtained from `compute_ecilag()` applied to the selected models
#' refitted on the complete dataset.
#'
#' Percent contribution is based on the absolute contribution magnitude within
#' each non-lag grouping unit:
#'
#' \deqn{100 |ECI_{ens,l}| / \sum_l |ECI_{ens,l}|.}
#'
#' Ensemble uncertainty for lag contributions is deliberately not produced in
#' EpiExposure v1. Independently sampled coefficient/posterior draws from
#' separately fitted models do not automatically form a valid joint draw
#' distribution, so arbitrary draw-index pairing is not performed.
#'
#' @export
ensemble_bestfit <- function(
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
    lambda = 1e-6,
    stack_intercept = TRUE,
    model_col = "model_id",
    id_cols = c("group", "fold"),
    lag_group_cols = c("var", "lag"),
    lg_strategy = c("requested_available", "model", "strict"),
    rl_weights = TRUE,
    verbose = TRUE
) {

  # ==========================================================================
  # ARGUMENT MATCHING AND SMALL VALIDATORS
  # ==========================================================================

  ensemble_scope <- match.arg(ensemble_scope)
  weight_transform <- match.arg(weight_transform)
  stacking_model <- match.arg(stacking_model)
  stack_objective <- match.arg(stack_objective)
  lg_strategy <- match.arg(lg_strategy)

  predictions_needed <-
    ensemble_scope %in% c("model", "both")
  lag_needed <-
    ensemble_scope %in% c("lag", "both")

  allowed_methods <- c(
    "unweighted",
    "weighted",
    "hard_voting",
    "best",
    "stacked"
  )

  if (!is.character(method) || !length(method) || anyNA(method) ||
      any(!nzchar(method))) {
    stop(
      "`method` must contain one or more non-empty ensemble method names.",
      call. = FALSE
    )
  }

  invalid_methods <- setdiff(method, allowed_methods)
  if (length(invalid_methods)) {
    stop(
      "Unsupported `method`: ",
      paste(invalid_methods, collapse = ", "),
      ". Available methods are: ",
      paste(allowed_methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  methods <- unique(method)

  response_col <- "observed"
  prediction_col <- "predicted"
  prediction_group_col <- "group"
  prediction_fold_col <- "fold"
  lag_col <- "lag"
  lag_contribution_col <- "ECI_weighted"
  lag_weight_col <- "weight"

  valid_names <- function(x) {
    is.character(x) &&
      length(x) >= 1L &&
      !anyNA(x) &&
      all(nzchar(x)) &&
      !anyDuplicated(x)
  }

  valid_scalar_name <- function(x) {
    valid_names(x) && length(x) == 1L
  }

  valid_flag <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  is_whole_scalar <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x) &&
      x == as.integer(x)
  }

  same_numeric <- function(a, b, tolerance = 1e-12) {
    is.numeric(a) &&
      is.numeric(b) &&
      length(a) == length(b) &&
      all(is.finite(a)) &&
      all(is.finite(b)) &&
      all(abs(a - b) <= tolerance * pmax(1, abs(a), abs(b)))
  }

  if (!is.data.frame(bestfit) || !nrow(bestfit)) {
    stop(
      "`bestfit` must be a non-empty data.frame returned by `find_bestfit()`.",
      call. = FALSE
    )
  }

  if (!valid_scalar_name(model_col)) {
    stop(
      "`model_col` must be one non-empty column name.",
      call. = FALSE
    )
  }

  if (!valid_names(id_cols)) {
    stop(
      "`id_cols` must contain unique non-empty column names.",
      call. = FALSE
    )
  }

  if (predictions_needed &&
      !all(c(prediction_group_col, prediction_fold_col) %in% id_cols)) {
    stop(
      "For model-level ensembles under the current EpiExposure CV contract, ",
      "`id_cols` must include both 'group' and 'fold'.",
      call. = FALSE
    )
  }

  if (!valid_names(lag_group_cols)) {
    stop(
      "`lag_group_cols` must contain unique non-empty column names.",
      call. = FALSE
    )
  }

  if (lag_needed &&
      !lag_col %in% lag_group_cols) {
    stop(
      "`lag_group_cols` must include 'lag' so each lag contribution is a ",
      "separate ensemble unit.",
      call. = FALSE
    )
  }

  if (!valid_scalar_name(group)) {
    stop(
      "`group` must be one non-empty column name.",
      call. = FALSE
    )
  }

  if (!is.null(var) &&
      (!is.character(var) || !length(var) || anyNA(var) ||
       any(!nzchar(var)) || anyDuplicated(var))) {
    stop(
      "`var` must be NULL or a unique non-empty character vector.",
      call. = FALSE
    )
  }

  if (!model_col %in% names(bestfit)) {
    stop(
      "`model_col` ('", model_col, "') was not found in `bestfit`.",
      call. = FALSE
    )
  }

  if (anyNA(bestfit[[model_col]]) ||
      anyDuplicated(as.character(bestfit[[model_col]]))) {
    stop(
      "`model_col` must contain unique non-missing model identifiers in ",
      "`bestfit`.",
      call. = FALSE
    )
  }

  if (!is_whole_scalar(top_n) || top_n < 1L) {
    stop(
      "`top_n` must be a positive integer.",
      call. = FALSE
    )
  }
  top_n <- as.integer(top_n)

  if (!is.numeric(lambda) || length(lambda) != 1L ||
      is.na(lambda) || !is.finite(lambda) || lambda < 0) {
    stop(
      "`lambda` must be one non-negative finite numeric value.",
      call. = FALSE
    )
  }

  if (!valid_flag(stack_intercept)) {
    stop(
      "`stack_intercept` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!valid_flag(rl_weights)) {
    stop(
      "`rl_weights` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!valid_flag(verbose)) {
    stop(
      "`verbose` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.list(compute_ecilag_args)) {
    stop(
      "`compute_ecilag_args` must be a named list.",
      call. = FALSE
    )
  }

  if (length(compute_ecilag_args)) {
    if (is.null(names(compute_ecilag_args)) ||
        anyNA(names(compute_ecilag_args)) ||
        any(!nzchar(names(compute_ecilag_args))) ||
        anyDuplicated(names(compute_ecilag_args))) {
      stop(
        "`compute_ecilag_args` must have unique non-empty names.",
        call. = FALSE
      )
    }

    reserved_args <- c(
      "data", "group", "fit", "var", "profile"
    )
    forbidden_args <- intersect(
      names(compute_ecilag_args),
      reserved_args
    )

    if (length(forbidden_args)) {
      stop(
        "`compute_ecilag_args` cannot override: ",
        paste(forbidden_args, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (lag_needed &&
        is.null(lag_data) &&
        "uncertainty" %in% names(compute_ecilag_args) &&
        isTRUE(compute_ecilag_args$uncertainty)) {
      stop(
        "Lag-ensemble uncertainty is not supported in EpiExposure v1. ",
        "Set `compute_ecilag_args$uncertainty = FALSE` or omit that argument.",
        call. = FALSE
      )
    }
  }

  if (!is.null(model_id)) {
    if (!length(model_id) || anyNA(model_id) ||
        anyDuplicated(as.character(model_id))) {
      stop(
        "`model_id` must be NULL or contain unique non-missing model ",
        "identifiers.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # STRICT BESTFIT METADATA CONTRACT
  # ==========================================================================

  required_attributes <- c(
    "family",
    "outcome_type",
    "rank_metric",
    "threshold",
    "prediction_level",
    "prediction_estimand",
    "prediction_contract",
    "cv_method",
    "cv_scheme",
    "fold_assignments",
    "basis_training_only",
    "max_lag",
    "history_length",
    "history_contract"
  )

  missing_attributes <- required_attributes[
    vapply(
      required_attributes,
      function(nm) is.null(attr(bestfit, nm, exact = TRUE)),
      logical(1)
    )
  ]

  if (length(missing_attributes)) {
    stop(
      "`bestfit` is missing required metadata from the current ",
      "`find_bestfit()`: ",
      paste(missing_attributes, collapse = ", "),
      ". Rerun model selection with the current function rather than inferring ",
      "an ensemble contract from an older object.",
      call. = FALSE
    )
  }

  family_name <- .resolve_family_name(
    attr(bestfit, "family", exact = TRUE)
  )

  supported_families <- c(
    "beta",
    "binomial",
    "poisson",
    "gamma",
    "gaussian",
    "negative_binomial"
  )

  if (identical(family_name, "ordinal")) {
    stop(
      "Ordinal outcomes are not supported in EpiExposure v1.",
      call. = FALSE
    )
  }

  if (!family_name %in% supported_families) {
    stop(
      "Unsupported EpiExposure v1 family in `bestfit`: '",
      family_name,
      "'.",
      call. = FALSE
    )
  }

  outcome_type <- .resolve_outcome_type(family_name)
  stored_outcome_type <- as.character(
    attr(bestfit, "outcome_type", exact = TRUE)
  )

  if (length(stored_outcome_type) != 1L ||
      is.na(stored_outcome_type) ||
      !identical(stored_outcome_type, outcome_type)) {
    stop(
      "`bestfit` contains inconsistent family/outcome metadata.",
      call. = FALSE
    )
  }

  prediction_level <- attr(
    bestfit,
    "prediction_level",
    exact = TRUE
  )
  prediction_estimand <- attr(
    bestfit,
    "prediction_estimand",
    exact = TRUE
  )
  prediction_contract <- attr(
    bestfit,
    "prediction_contract",
    exact = TRUE
  )

  if (!identical(prediction_level, "population")) {
    stop(
      "`ensemble_bestfit()` requires population-level OOF predictions.",
      call. = FALSE
    )
  }

  if (!identical(prediction_estimand, "expected_response")) {
    stop(
      "`ensemble_bestfit()` requires the expected-response prediction estimand.",
      call. = FALSE
    )
  }

  if (!identical(
    prediction_contract,
    "central_expected_response"
  )) {
    stop(
      "`ensemble_bestfit()` requires deterministic central expected-response ",
      "OOF predictions from the current `find_bestfit()`.",
      call. = FALSE
    )
  }

  if (!isTRUE(
    attr(bestfit, "basis_training_only", exact = TRUE)
  )) {
    stop(
      "`bestfit` does not certify training-only basis construction.",
      call. = FALSE
    )
  }

  cv_method <- as.character(
    attr(bestfit, "cv_method", exact = TRUE)
  )
  cv_scheme <- as.character(
    attr(bestfit, "cv_scheme", exact = TRUE)
  )

  if (length(cv_method) != 1L ||
      is.na(cv_method) ||
      !cv_method %in% c("LOOCV", "k-fold")) {
    stop(
      "`bestfit` contains invalid `cv_method` metadata.",
      call. = FALSE
    )
  }

  valid_cv_schemes <- c(
    LOOCV = "leave_one_group_out",
    `k-fold` = if (identical(outcome_type, "binary")) {
      "stratified_grouped_k_fold"
    } else {
      "grouped_k_fold"
    }
  )

  if (length(cv_scheme) != 1L ||
      is.na(cv_scheme) ||
      !identical(
        cv_scheme,
        unname(valid_cv_schemes[[cv_method]])
      )) {
    stop(
      "`bestfit` contains inconsistent `cv_method`/`cv_scheme` metadata.",
      call. = FALSE
    )
  }

  effective_k <- attr(bestfit, "k", exact = TRUE)

  if (is.null(effective_k) ||
      !is_whole_scalar(effective_k) ||
      effective_k < 2L) {
    stop(
      "`bestfit` contains invalid cross-validation `k` metadata.",
      call. = FALSE
    )
  }
  effective_k <- as.integer(effective_k)

  cv_stratified <- attr(
    bestfit,
    "cv_stratified",
    exact = TRUE
  )

  if (is.null(cv_stratified) ||
      !valid_flag(cv_stratified)) {
    stop(
      "`bestfit` contains invalid `cv_stratified` metadata.",
      call. = FALSE
    )
  }

  expected_stratified <-
    identical(cv_method, "k-fold") &&
    identical(outcome_type, "binary")

  if (!identical(cv_stratified, expected_stratified)) {
    stop(
      "`bestfit` contains inconsistent binomial CV stratification metadata.",
      call. = FALSE
    )
  }

  fold_assignments <- attr(
    bestfit,
    "fold_assignments",
    exact = TRUE
  )

  if (!is.data.frame(fold_assignments) ||
      !nrow(fold_assignments) ||
      !all(
        c("group", "fold", "observed") %in%
        names(fold_assignments)
      )) {
    stop(
      "`bestfit` contains invalid `fold_assignments` metadata.",
      call. = FALSE
    )
  }

  if (anyNA(fold_assignments$group) ||
      anyDuplicated(as.character(fold_assignments$group)) ||
      !is.numeric(fold_assignments$fold) ||
      anyNA(fold_assignments$fold) ||
      any(fold_assignments$fold !=
          as.integer(fold_assignments$fold)) ||
      any(fold_assignments$fold < 1L |
          fold_assignments$fold > effective_k) ||
      !is.numeric(fold_assignments$observed) ||
      anyNA(fold_assignments$observed) ||
      any(!is.finite(fold_assignments$observed))) {
    stop(
      "`bestfit$fold_assignments` metadata are internally invalid.",
      call. = FALSE
    )
  }

  fold_assignments$group <- as.character(
    fold_assignments$group
  )
  fold_assignments$fold <- as.integer(
    fold_assignments$fold
  )

  fold_balance <- attr(
    bestfit,
    "fold_balance",
    exact = TRUE
  )

  if (!is.null(fold_balance) &&
      !is.data.frame(fold_balance)) {
    stop(
      "`bestfit` contains invalid `fold_balance` metadata.",
      call. = FALSE
    )
  }

  cv_seed <- attr(bestfit, "cv_seed", exact = TRUE)
  max_lag_metadata <- attr(bestfit, "max_lag", exact = TRUE)
  history_length_metadata <- attr(
    bestfit,
    "history_length",
    exact = TRUE
  )
  history_contract_metadata <- attr(
    bestfit,
    "history_contract",
    exact = TRUE
  )
  time_step_metadata <- attr(bestfit, "time_step", exact = TRUE)

  if (!is_whole_scalar(max_lag_metadata) ||
      max_lag_metadata < 0L) {
    stop(
      "`bestfit` contains invalid common `max_lag` metadata.",
      call. = FALSE
    )
  }
  max_lag_metadata <- as.integer(max_lag_metadata)

  if (!is_whole_scalar(history_length_metadata) ||
      history_length_metadata < 1L) {
    stop(
      "`bestfit` contains invalid `history_length` metadata.",
      call. = FALSE
    )
  }
  history_length_metadata <- as.integer(history_length_metadata)

  if (!identical(
    history_length_metadata,
    max_lag_metadata + 1L
  )) {
    stop(
      "`bestfit$history_length` metadata must equal `max_lag + 1`. ",
      "Expected ", max_lag_metadata + 1L,
      " but found ", history_length_metadata, ".",
      call. = FALSE
    )
  }

  if (!identical(
    history_contract_metadata,
    "all_fitted_exposures_same_exact_max_lag_plus_one"
  )) {
    stop(
      "`bestfit` does not satisfy the current EpiExposure exact-history ",
      "contract.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # METRIC AND THRESHOLD CONTRACT
  # ==========================================================================

  available_metrics <- .available_metrics(family_name)

  inherited_metric <- attr(
    bestfit,
    "rank_metric",
    exact = TRUE
  )

  if (is.null(weight_metric)) {
    if (!is.character(inherited_metric) ||
        length(inherited_metric) != 1L ||
        is.na(inherited_metric) ||
        !nzchar(inherited_metric)) {
      stop(
        "`bestfit` does not contain one valid `rank_metric`.",
        call. = FALSE
      )
    }

    weight_metric <- inherited_metric
  }

  if (!is.character(weight_metric) ||
      length(weight_metric) != 1L ||
      is.na(weight_metric) ||
      !nzchar(weight_metric)) {
    stop(
      "`weight_metric` must be NULL or one non-empty metric name.",
      call. = FALSE
    )
  }

  if (!weight_metric %in% available_metrics) {
    stop(
      "`weight_metric = ",
      sQuote(weight_metric),
      "` is not available for family '",
      family_name,
      "'. Available metrics are: ",
      paste(available_metrics, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!weight_metric %in% names(bestfit)) {
    stop(
      "`weight_metric` ('",
      weight_metric,
      "') was not found in `bestfit`.",
      call. = FALSE
    )
  }

  if (!is.numeric(bestfit[[weight_metric]])) {
    stop(
      "`weight_metric` must identify a numeric column in `bestfit`.",
      call. = FALSE
    )
  }

  metric_direction <- .metric_direction(weight_metric)

  stored_threshold <- attr(
    bestfit,
    "threshold",
    exact = TRUE
  )

  if (!is.numeric(stored_threshold) ||
      length(stored_threshold) != 1L ||
      is.na(stored_threshold) ||
      !is.finite(stored_threshold) ||
      stored_threshold <= 0 ||
      stored_threshold >= 1) {
    stop(
      "`bestfit` contains invalid threshold metadata.",
      call. = FALSE
    )
  }

  threshold_was_supplied <- !is.null(threshold)

  if (is.null(threshold)) {
    threshold <- as.numeric(stored_threshold)
  }

  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      is.na(threshold) ||
      !is.finite(threshold) ||
      threshold <= 0 ||
      threshold >= 1) {
    stop(
      "`threshold` must be one finite numeric probability strictly between ",
      "0 and 1.",
      call. = FALSE
    )
  }

  classification_metrics <- c(
    "Accuracy",
    "Balanced_Accuracy",
    "Sensitivity",
    "Specificity",
    "F1",
    "MCC",
    "Precision"
  )

  if (identical(outcome_type, "binary") &&
      weight_metric %in% classification_metrics &&
      threshold_was_supplied &&
      !same_numeric(threshold, stored_threshold)) {
    stop(
      "`weight_metric = '",
      weight_metric,
      "' was calculated by `find_bestfit()` using threshold = ",
      format(stored_threshold),
      ", but `ensemble_bestfit()` received threshold = ",
      format(threshold),
      ". Use the stored threshold, choose a threshold-independent probability ",
      "metric, or rerun `find_bestfit()` with the desired threshold.",
      call. = FALSE
    )
  }

  if (!identical(outcome_type, "binary")) {
    if (threshold_was_supplied && verbose) {
      message(
        "`threshold` is ignored for non-binomial ensemble performance."
      )
    }

    # Keep a valid internal value for the common performance helper, but report
    # threshold as NA in non-binary ensemble summaries.
    threshold <- as.numeric(stored_threshold)
  }

  # ==========================================================================
  # METHOD-SPECIFIC CONTRACTS
  # ==========================================================================

  if (ensemble_scope %in% c("lag", "both") &&
      "stacked" %in% methods) {
    stop(
      "`method = 'stacked'` is supported only for ",
      "`ensemble_scope = 'model'`.",
      call. = FALSE
    )
  }

  if (ensemble_scope %in% c("lag", "both") &&
      "hard_voting" %in% methods) {
    stop(
      "`method = 'hard_voting'` is supported only for ",
      "`ensemble_scope = 'model'`.",
      call. = FALSE
    )
  }

  if ("hard_voting" %in% methods &&
      !identical(outcome_type, "binary")) {
    stop(
      "`method = 'hard_voting'` is available only for binomial outcomes.",
      call. = FALSE
    )
  }

  if ("stacked" %in% methods) {
    if (identical(stack_objective, "regularized") &&
        !identical(stacking_model, "ridge")) {
      stop(
        "`stack_objective = 'regularized'` requires ",
        "`stacking_model = 'ridge'`.",
        call. = FALSE
      )
    }

    if (identical(stack_objective, "min_error") &&
        !identical(stacking_model, "lm")) {
      stop(
        "`stack_objective = 'min_error'` requires ",
        "`stacking_model = 'lm'`.",
        call. = FALSE
      )
    }

    if (identical(stack_objective, "constrained") &&
        isTRUE(stack_intercept)) {
      stop(
        "`stack_objective = 'constrained'` requires ",
        "`stack_intercept = FALSE` so the ensemble remains a convex ",
        "combination of model predictions.",
        call. = FALSE
      )
    }
  }

  if (verbose) {
    message(
      "Ensemble contract: family = ",
      family_name,
      ", outcome_type = ",
      outcome_type,
      ", CV = ",
      cv_method,
      if (identical(cv_method, "k-fold")) {
        paste0(" (k=", effective_k, ")")
      } else {
        ""
      },
      ", weight_metric = ",
      weight_metric,
      " (",
      metric_direction,
      ")."
    )

    if (identical(outcome_type, "binary") &&
        ensemble_scope %in% c("model", "both")) {
      message(
        "Classification threshold = ",
        format(threshold),
        ". Probability metrics retain probabilities; thresholding is used ",
        "only for classification metrics",
        if ("hard_voting" %in% methods) " and model votes" else "",
        "."
      )
    }
  }

  # ==========================================================================
  # MODEL SELECTION
  # ==========================================================================

  rank_models <- function(model_table) {
    if (!nrow(model_table)) {
      return(model_table)
    }

    metric <- model_table[[weight_metric]]
    finite <- is.finite(metric)

    if (!any(finite)) {
      return(model_table[FALSE, , drop = FALSE])
    }

    model_table <- model_table[finite, , drop = FALSE]
    metric <- model_table[[weight_metric]]

    primary <- if (identical(
      metric_direction,
      "maximize"
    )) {
      -metric
    } else {
      metric
    }

    tie_rank <- if (
      "rank" %in% names(model_table) &&
      is.numeric(model_table$rank)
    ) {
      model_table$rank
    } else {
      seq_len(nrow(model_table))
    }

    ord <- order(
      primary,
      tie_rank,
      as.character(model_table[[model_col]]),
      na.last = TRUE
    )

    model_table[ord, , drop = FALSE]
  }

  bestfit_ids_chr <- as.character(
    bestfit[[model_col]]
  )

  if (!is.null(model_id)) {
    requested_ids_chr <- as.character(model_id)

    missing_models <- setdiff(
      requested_ids_chr,
      bestfit_ids_chr
    )

    if (length(missing_models)) {
      stop(
        "The following `model_id` were not found in `bestfit`: ",
        paste(missing_models, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    requested_index <- match(
      requested_ids_chr,
      bestfit_ids_chr
    )

    selected_models <- bestfit[
      requested_index,
      ,
      drop = FALSE
    ]

    nonfinite_selected <- !is.finite(
      selected_models[[weight_metric]]
    )

    if (any(nonfinite_selected)) {
      stop(
        "Explicitly selected model(s) have non-finite `",
        weight_metric,
        "`: ",
        paste(
          selected_models[[model_col]][
            nonfinite_selected
          ],
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    selected_models <- rank_models(
      selected_models
    )
  } else {
    ranked_finite <- rank_models(bestfit)

    if (!nrow(ranked_finite)) {
      stop(
        "No candidate model has a finite `",
        weight_metric,
        "` value.",
        call. = FALSE
      )
    }

    selected_models <- utils::head(
      ranked_finite,
      top_n
    )

    if (verbose &&
        nrow(ranked_finite) < top_n) {
      message(
        "Only ",
        nrow(ranked_finite),
        " model(s) have finite ",
        weight_metric,
        "; all available models will be used."
      )
    }
  }

  if (!nrow(selected_models)) {
    stop(
      "No models were selected for the ensemble.",
      call. = FALSE
    )
  }

  selected_ids <- selected_models[[model_col]]
  selected_ids_chr <- as.character(selected_ids)
  n_selected <- length(selected_ids_chr)

  if (verbose) {
    message(
      "Selected ",
      n_selected,
      " model(s): ",
      paste(selected_ids_chr, collapse = ", "),
      ". Methods: ",
      paste(methods, collapse = ", "),
      "."
    )
  }

  # ==========================================================================
  # MODEL WEIGHT TRANSFORMS
  # ==========================================================================

  compute_model_weights <- function(
    model_table,
    ensemble_method
  ) {
    n_models <- nrow(model_table)

    if (n_models < 1L) {
      stop(
        "No models are available for ensemble weighting.",
        call. = FALSE
      )
    }

    ids <- as.character(
      model_table[[model_col]]
    )

    if (identical(ensemble_method, "unweighted")) {
      weights <- rep(
        1 / n_models,
        n_models
      )

    } else if (identical(ensemble_method, "best")) {
      weights <- c(
        1,
        rep(0, n_models - 1L)
      )

    } else if (identical(ensemble_method, "weighted")) {
      metric <- model_table[[weight_metric]]

      if (any(!is.finite(metric))) {
        stop(
          "Weighted ensembles require finite `",
          weight_metric,
          "` values for every selected model.",
          call. = FALSE
        )
      }

      if (identical(weight_transform, "uniform")) {
        weights <- rep(
          1 / n_models,
          n_models
        )

      } else if (identical(
        weight_transform,
        "rank_inverse"
      )) {
        utility <- 1 / seq_len(n_models)
        weights <- utility / sum(utility)

      } else if (identical(
        weight_transform,
        "positive"
      )) {
        if (identical(
          metric_direction,
          "maximize"
        )) {
          utility <- pmax(metric, 0)

          if (!any(utility > 0)) {
            stop(
              "`weight_transform = 'positive'` produced no positive utility ",
              "because every selected value of `",
              weight_metric,
              "` is less than or equal to zero. Use `softmax`, ",
              "`rank_inverse`, or `uniform` instead.",
              call. = FALSE
            )
          }

        } else {
          if (any(metric < 0)) {
            stop(
              "`weight_transform = 'positive'` expects non-negative loss/error ",
              "values for minimize metrics, but `",
              weight_metric,
              "` contains a negative value.",
              call. = FALSE
            )
          }

          utility <- 1 / pmax(
            metric,
            .Machine$double.eps
          )
        }

        if (any(!is.finite(utility)) ||
            sum(utility) <= 0) {
          stop(
            "Could not construct finite positive model utilities.",
            call. = FALSE
          )
        }

        weights <- utility / sum(utility)

      } else {
        score <- if (identical(
          metric_direction,
          "maximize"
        )) {
          metric
        } else {
          -metric
        }

        shifted <- score - max(score)
        exp_score <- exp(shifted)

        if (any(!is.finite(exp_score)) ||
            sum(exp_score) <= 0) {
          stop(
            "Softmax model weighting failed to produce finite weights.",
            call. = FALSE
          )
        }

        weights <- exp_score / sum(exp_score)
      }

    } else {
      stop(
        "Internal model weighting is not defined for method '",
        ensemble_method,
        "'.",
        call. = FALSE
      )
    }

    if (length(weights) != n_models ||
        any(!is.finite(weights)) ||
        any(weights < 0) ||
        abs(sum(weights) - 1) >
        1e-10 * max(1, n_models)) {
      stop(
        "Internal ensemble weights are invalid.",
        call. = FALSE
      )
    }

    names(weights) <- ids
    weights
  }

  # ==========================================================================
  # OOF PREDICTION ALIGNMENT
  # ==========================================================================

  make_key <- function(d, columns) {
    if (!nrow(d)) {
      return(character(0))
    }

    encoded <- lapply(
      columns,
      function(column) {
        if (!column %in% names(d)) {
          stop(
            "Key column '",
            column,
            "' is absent.",
            call. = FALSE
          )
        }

        value <- d[[column]]

        if (anyNA(value)) {
          stop(
            "Key column '",
            column,
            "' cannot contain missing values.",
            call. = FALSE
          )
        }

        as.character(value)
      }
    )

    do.call(
      paste,
      c(
        encoded,
        list(sep = "\u001f")
      )
    )
  }

  validate_predictions_against_cv <- function(
    prediction_data
  ) {
    required_cv <- c(
      prediction_group_col,
      prediction_fold_col,
      response_col
    )

    missing_cv <- setdiff(
      required_cv,
      names(prediction_data)
    )

    if (length(missing_cv)) {
      stop(
        "OOF predictions are missing current CV-contract column(s): ",
        paste(missing_cv, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    cv_rows <- unique(
      prediction_data[
        ,
        required_cv,
        drop = FALSE
      ]
    )

    group_split <- split(
      cv_rows,
      as.character(
        cv_rows[[prediction_group_col]]
      )
    )

    bad_group <- vapply(
      group_split,
      function(x) nrow(x) != 1L,
      logical(1)
    )

    if (any(bad_group)) {
      stop(
        "At least one OOF group is associated with multiple folds or observed ",
        "values in `predictions`.",
        call. = FALSE
      )
    }

    match_index <- match(
      as.character(cv_rows[[prediction_group_col]]),
      fold_assignments$group
    )

    if (anyNA(match_index)) {
      stop(
        "`predictions` contains group(s) absent from ",
        "`attr(bestfit, 'fold_assignments')`.",
        call. = FALSE
      )
    }

    expected_fold <- fold_assignments$fold[
      match_index
    ]
    expected_observed <- fold_assignments$observed[
      match_index
    ]

    if (any(
      as.integer(cv_rows[[prediction_fold_col]]) !=
      expected_fold
    )) {
      stop(
        "OOF fold identifiers in `predictions` do not match ",
        "`attr(bestfit, 'fold_assignments')`.",
        call. = FALSE
      )
    }

    observed_value <- as.numeric(
      cv_rows[[response_col]]
    )

    if (any(!is.finite(observed_value)) ||
        !same_numeric(
          observed_value,
          as.numeric(expected_observed)
        )) {
      stop(
        "Observed outcomes in `predictions` do not match ",
        "`attr(bestfit, 'fold_assignments')`.",
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  build_prediction_matrix <- function(
    prediction_data
  ) {
    required <- unique(
      c(
        model_col,
        id_cols,
        prediction_group_col,
        prediction_fold_col,
        response_col,
        prediction_col
      )
    )

    missing <- setdiff(
      required,
      names(prediction_data)
    )

    if (length(missing)) {
      stop(
        "Missing OOF prediction column(s): ",
        paste(missing, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    prediction_data <- prediction_data[
      as.character(
        prediction_data[[model_col]]
      ) %in% selected_ids_chr,
      ,
      drop = FALSE
    ]

    if (!nrow(prediction_data)) {
      stop(
        "No OOF predictions were found for the selected models.",
        call. = FALSE
      )
    }

    present_models <- unique(
      as.character(
        prediction_data[[model_col]]
      )
    )

    missing_selected_predictions <- setdiff(
      selected_ids_chr,
      present_models
    )

    if (length(missing_selected_predictions)) {
      stop(
        "OOF predictions are missing for selected model(s): ",
        paste(
          missing_selected_predictions,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    if (!is.numeric(prediction_data[[response_col]]) ||
        anyNA(prediction_data[[response_col]]) ||
        any(!is.finite(
          prediction_data[[response_col]]
        ))) {
      stop(
        "`observed` in OOF predictions must contain only finite numeric values.",
        call. = FALSE
      )
    }

    if (!is.numeric(prediction_data[[prediction_col]]) ||
        anyNA(prediction_data[[prediction_col]]) ||
        any(!is.finite(
          prediction_data[[prediction_col]]
        ))) {
      stop(
        "`predicted` in OOF predictions must contain only finite numeric values.",
        call. = FALSE
      )
    }

    key_cols <- c(
      model_col,
      id_cols
    )

    if (anyDuplicated(
      prediction_data[key_cols]
    )) {
      stop(
        "`predictions` contains duplicate rows for the same model and OOF ",
        "identifier.",
        call. = FALSE
      )
    }

    validate_predictions_against_cv(
      prediction_data
    )

    prediction_data$.obs_key <- make_key(
      prediction_data,
      id_cols
    )

    # Strictly require identical OOF support across selected models.
    support_by_model <- lapply(
      selected_ids_chr,
      function(model_id) {
        sort(
          unique(
            prediction_data$.obs_key[
              as.character(
                prediction_data[[model_col]]
              ) == model_id
            ]
          )
        )
      }
    )
    names(support_by_model) <- selected_ids_chr

    reference_support <- support_by_model[[1L]]

    for (model_id in selected_ids_chr[-1L]) {
      current_support <- support_by_model[[model_id]]

      if (!identical(
        current_support,
        reference_support
      )) {
        missing_current <- setdiff(
          reference_support,
          current_support
        )
        extra_current <- setdiff(
          current_support,
          reference_support
        )

        stop(
          "Selected models do not share identical OOF validation support. ",
          "Model ",
          model_id,
          " differs from model ",
          selected_ids_chr[1L],
          " (missing identifier(s): ",
          length(missing_current),
          "; extra identifier(s): ",
          length(extra_current),
          "). Choose models with the same successful OOF groups or rerun ",
          "`find_bestfit()` so selected candidates complete the same CV support.",
          call. = FALSE
        )
      }
    }

    reference_model <- selected_ids_chr[1L]

    obs_df <- prediction_data[
      as.character(
        prediction_data[[model_col]]
      ) == reference_model,
      c(
        ".obs_key",
        id_cols,
        response_col
      ),
      drop = FALSE
    ]

    reference_order <- order(
      match(
        as.character(
          obs_df[[prediction_group_col]]
        ),
        fold_assignments$group
      ),
      obs_df[[prediction_fold_col]]
    )

    obs_df <- obs_df[
      reference_order,
      ,
      drop = FALSE
    ]
    rownames(obs_df) <- NULL

    obs_keys <- obs_df$.obs_key
    n_obs <- length(obs_keys)

    if (n_obs < 2L) {
      stop(
        "At least two aligned OOF groups are required for ensemble evaluation.",
        call. = FALSE
      )
    }

    pred_mat <- matrix(
      NA_real_,
      nrow = n_obs,
      ncol = n_selected,
      dimnames = list(
        obs_keys,
        selected_ids_chr
      )
    )

    for (j in seq_along(selected_ids_chr)) {
      model_id <- selected_ids_chr[j]

      temporary <- prediction_data[
        as.character(
          prediction_data[[model_col]]
        ) == model_id,
        ,
        drop = FALSE
      ]

      idx <- match(
        obs_keys,
        temporary$.obs_key
      )

      if (anyNA(idx)) {
        stop(
          "Internal OOF support alignment failed for model ",
          model_id,
          ".",
          call. = FALSE
        )
      }

      model_observed <- temporary[[response_col]][
        idx
      ]

      if (!same_numeric(
        as.numeric(model_observed),
        as.numeric(obs_df[[response_col]])
      )) {
        stop(
          "Observed outcomes differ across selected models for the same OOF ",
          "identifier.",
          call. = FALSE
        )
      }

      pred_mat[, j] <- temporary[[prediction_col]][
        idx
      ]
    }

    if (any(!is.finite(pred_mat))) {
      stop(
        "Aligned OOF prediction matrix contains non-finite values.",
        call. = FALSE
      )
    }

    observed <- as.numeric(
      obs_df[[response_col]]
    )
    fold_vector <- as.integer(
      obs_df[[prediction_fold_col]]
    )

    if (length(unique(fold_vector)) < 2L) {
      stop(
        "At least two distinct CV folds are required for ensemble evaluation.",
        call. = FALSE
      )
    }

    if (identical(outcome_type, "binary")) {
      if (!all(observed %in% c(0, 1))) {
        stop(
          "Binomial ensemble evaluation requires observed outcomes coded 0/1.",
          call. = FALSE
        )
      }

      if (any(pred_mat < 0 | pred_mat > 1)) {
        stop(
          "Binomial OOF model predictions must be probabilities in [0, 1].",
          call. = FALSE
        )
      }

      if (all(pred_mat %in% c(0, 1))) {
        warning(
          "Selected binomial-model OOF predictions contain only 0/1 values. ",
          "Probability-based metrics will have no probability-resolution ",
          "information.",
          call. = FALSE
        )
      }
    }

    colnames(pred_mat) <- paste0(
      "model_",
      selected_ids_chr
    )

    list(
      pred_mat = pred_mat,
      obs_df = obs_df,
      obs = observed,
      folds = fold_vector,
      n_oof = n_obs
    )
  }

  prediction_matrix <- NULL

  if (predictions_needed) {
    if (is.null(predictions)) {
      predictions <- attr(
        bestfit,
        "predictions",
        exact = TRUE
      )
    }

    if (!is.data.frame(predictions) ||
        !nrow(predictions)) {
      stop(
        "`predictions` must be supplied or stored as non-empty ",
        "`attr(bestfit, 'predictions')`.",
        call. = FALSE
      )
    }

    prediction_matrix <- build_prediction_matrix(
      predictions
    )
  }

  # ==========================================================================
  # STACKING
  # ==========================================================================

  project_simplex <- function(v) {
    v <- as.numeric(v)

    if (!length(v) ||
        any(!is.finite(v))) {
      stop(
        "Simplex projection requires finite numeric values.",
        call. = FALSE
      )
    }

    u <- sort(
      v,
      decreasing = TRUE
    )
    cssv <- cumsum(u)

    candidates <- which(
      u + (1 - cssv) /
        seq_along(u) > 0
    )

    if (!length(candidates)) {
      stop(
        "Simplex projection failed.",
        call. = FALSE
      )
    }

    rho <- max(candidates)
    theta <- (
      cssv[rho] - 1
    ) / rho

    w <- pmax(
      v - theta,
      0
    )

    total <- sum(w)

    if (!is.finite(total) ||
        total <= 0) {
      stop(
        "Simplex projection produced invalid weights.",
        call. = FALSE
      )
    }

    w / total
  }

  fit_constrained_stack <- function(
    X,
    y
  ) {
    X <- as.matrix(X)
    y <- as.numeric(y)

    if (nrow(X) != length(y) ||
        !nrow(X) ||
        !ncol(X) ||
        any(!is.finite(X)) ||
        any(!is.finite(y))) {
      stop(
        "Invalid data supplied to constrained stacking.",
        call. = FALSE
      )
    }

    p <- ncol(X)

    if (p == 1L) {
      return(
        list(
          intercept = 0,
          weights = 1,
          converged = TRUE,
          iterations = 0L
        )
      )
    }

    penalty_lambda <- if (
      identical(stacking_model, "ridge")
    ) {
      lambda
    } else {
      0
    }

    n <- nrow(X)
    H <- crossprod(X) / n
    if (penalty_lambda > 0) {
      H <- H + diag(
        penalty_lambda,
        p
      )
    }

    b <- as.numeric(
      crossprod(X, y) / n
    )

    eigen_values <- tryCatch(
      eigen(
        H,
        symmetric = TRUE,
        only.values = TRUE
      )$values,
      error = function(e) NULL
    )

    if (is.null(eigen_values) ||
        any(!is.finite(eigen_values))) {
      stop(
        "Could not determine a stable constrained-stacking step size.",
        call. = FALSE
      )
    }

    lipschitz <- 2 * max(
      eigen_values,
      0
    )

    if (!is.finite(lipschitz) ||
        lipschitz <= 0) {
      lipschitz <- 1
    }

    step_size <- 1 / lipschitz
    w <- rep(
      1 / p,
      p
    )

    tolerance <- 1e-10
    max_iterations <- 10000L
    converged <- FALSE
    iteration <- 0L

    for (iteration in seq_len(max_iterations)) {
      gradient <- 2 * (
        as.numeric(H %*% w) - b
      )

      w_new <- project_simplex(
        w - step_size * gradient
      )

      if (max(abs(w_new - w)) <= tolerance) {
        w <- w_new
        converged <- TRUE
        break
      }

      w <- w_new
    }

    if (!converged) {
      stop(
        "Constrained stacking did not converge within ",
        max_iterations,
        " projected-gradient iterations.",
        call. = FALSE
      )
    }

    list(
      intercept = 0,
      weights = w,
      converged = TRUE,
      iterations = as.integer(iteration)
    )
  }

  fit_stack_model <- function(
    X,
    y
  ) {
    X <- as.matrix(X)
    y <- as.numeric(y)

    if (nrow(X) != length(y) ||
        !nrow(X) ||
        !ncol(X) ||
        any(!is.finite(X)) ||
        any(!is.finite(y))) {
      stop(
        "Invalid data supplied to stacking.",
        call. = FALSE
      )
    }

    if (identical(
      stack_objective,
      "constrained"
    )) {
      fit <- fit_constrained_stack(
        X,
        y
      )

      names(fit$weights) <- selected_ids_chr
      return(fit)
    }

    design <- if (stack_intercept) {
      cbind(
        "(Intercept)" = 1,
        X
      )
    } else {
      X
    }

    if (identical(
      stack_objective,
      "min_error"
    )) {
      lm_fit <- stats::lm.fit(
        x = design,
        y = y
      )

      if (lm_fit$rank < ncol(design)) {
        stop(
          "Least-squares stacking is rank deficient. Use fewer selected ",
          "models or `stacking_model = 'ridge', ",
          "stack_objective = 'regularized'`.",
          call. = FALSE
        )
      }

      coefficients <- as.numeric(
        lm_fit$coefficients
      )

      if (length(coefficients) !=
          ncol(design) ||
          any(!is.finite(coefficients))) {
        stop(
          "Least-squares stacking produced invalid coefficients.",
          call. = FALSE
        )
      }

    } else {
      penalty <- diag(
        lambda,
        ncol(design)
      )

      if (stack_intercept) {
        penalty[1L, 1L] <- 0
      }

      A <- crossprod(design) + penalty
      b <- crossprod(design, y)

      coefficients <- tryCatch(
        as.numeric(
          solve(A, b)
        ),
        error = function(e) {
          stop(
            "Ridge stacking failed to solve the penalized normal equations: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )

      if (length(coefficients) !=
          ncol(design) ||
          any(!is.finite(coefficients))) {
        stop(
          "Ridge stacking produced invalid coefficients.",
          call. = FALSE
        )
      }
    }

    intercept <- if (stack_intercept) {
      coefficients[1L]
    } else {
      0
    }

    weights <- if (stack_intercept) {
      coefficients[-1L]
    } else {
      coefficients
    }

    if (length(weights) !=
        ncol(X) ||
        any(!is.finite(weights)) ||
        !is.finite(intercept)) {
      stop(
        "Stacking coefficient extraction failed.",
        call. = FALSE
      )
    }

    names(weights) <- selected_ids_chr

    list(
      intercept = as.numeric(intercept),
      weights = as.numeric(weights),
      converged = TRUE,
      iterations = NA_integer_
    )
  }

  predict_stack_model <- function(
    stack_fit,
    X
  ) {
    X <- as.matrix(X)

    if (ncol(X) != length(
      stack_fit$weights
    )) {
      stop(
        "Stacking prediction dimension mismatch.",
        call. = FALSE
      )
    }

    prediction <- as.numeric(
      stack_fit$intercept +
        X %*% as.numeric(
          stack_fit$weights
        )
    )

    if (any(!is.finite(prediction))) {
      stop(
        "Stacking produced non-finite predictions.",
        call. = FALSE
      )
    }

    prediction
  }

  validate_binomial_stack_prediction <- function(
    prediction,
    stage
  ) {
    if (!identical(outcome_type, "binary")) {
      return(prediction)
    }

    tolerance <- 1e-12

    if (any(
      prediction < -tolerance |
      prediction > 1 + tolerance
    )) {
      stop(
        "Stacking produced values outside [0, 1] for a binomial outcome ",
        "during ",
        stage,
        ". These values cannot be interpreted as probabilities. Use ",
        "`stack_objective = 'constrained'` with ",
        "`stack_intercept = FALSE`, reduce model collinearity, or choose ",
        "another ensemble method.",
        call. = FALSE
      )
    }

    pmin(
      pmax(prediction, 0),
      1
    )
  }

  compute_stacking <- function(
    pred_mat,
    observed,
    folds
  ) {
    pred_mat <- as.matrix(pred_mat)
    observed <- as.numeric(observed)
    folds <- as.integer(folds)

    if (nrow(pred_mat) !=
        length(observed) ||
        length(observed) !=
        length(folds)) {
      stop(
        "Internal stacking dimensions are inconsistent.",
        call. = FALSE
      )
    }

    fold_levels <- sort(
      unique(folds)
    )

    if (length(fold_levels) < 2L) {
      stop(
        "Stacking requires at least two CV folds.",
        call. = FALSE
      )
    }

    cross_fitted <- rep(
      NA_real_,
      nrow(pred_mat)
    )

    coefficient_rows <- list()
    coefficient_index <- 1L

    for (fold_value in fold_levels) {
      test_index <- which(
        folds == fold_value
      )
      train_index <- which(
        folds != fold_value
      )

      if (!length(test_index) ||
          length(train_index) < 2L) {
        stop(
          "Insufficient level-2 training rows while cross-fitting stacker for ",
          "fold ",
          fold_value,
          ".",
          call. = FALSE
        )
      }

      fold_fit <- fit_stack_model(
        pred_mat[
          train_index,
          ,
          drop = FALSE
        ],
        observed[train_index]
      )

      fold_prediction <- predict_stack_model(
        fold_fit,
        pred_mat[
          test_index,
          ,
          drop = FALSE
        ]
      )

      fold_prediction <-
        validate_binomial_stack_prediction(
          fold_prediction,
          paste0(
            "fold-cross-fitted evaluation of fold ",
            fold_value
          )
        )

      cross_fitted[test_index] <-
        fold_prediction

      coefficient_row <- data.frame(
        fold = rep(
          fold_value,
          n_selected
        ),
        selected_ids,
        weight = as.numeric(
          fold_fit$weights
        ),
        intercept = rep(
          fold_fit$intercept,
          n_selected
        ),
        n_level2_train = rep(
          length(train_index),
          n_selected
        ),
        n_level2_test = rep(
          length(test_index),
          n_selected
        ),
        iterations = rep(
          fold_fit$iterations,
          n_selected
        ),
        stringsAsFactors = FALSE
      )
      names(coefficient_row)[2L] <- model_col

      coefficient_rows[[
        coefficient_index
      ]] <- coefficient_row

      coefficient_index <-
        coefficient_index + 1L
    }

    if (any(!is.finite(cross_fitted))) {
      stop(
        "Stacking failed to produce a finite cross-fitted prediction for every ",
        "OOF row.",
        call. = FALSE
      )
    }

    final_fit <- fit_stack_model(
      pred_mat,
      observed
    )

    final_prediction <- predict_stack_model(
      final_fit,
      pred_mat
    )

    final_prediction <-
      validate_binomial_stack_prediction(
        final_prediction,
        "final level-2 fit"
      )

    final_coefficients <- data.frame(
      selected_ids,
      weight = as.numeric(
        final_fit$weights
      ),
      intercept = rep(
        final_fit$intercept,
        n_selected
      ),
      n_level2_train = rep(
        nrow(pred_mat),
        n_selected
      ),
      iterations = rep(
        final_fit$iterations,
        n_selected
      ),
      stringsAsFactors = FALSE
    )
    names(final_coefficients)[1L] <- model_col

    by_fold <- do.call(
      rbind,
      coefficient_rows
    )
    rownames(by_fold) <- NULL
    rownames(final_coefficients) <- NULL

    list(
      cross_fitted_prediction = cross_fitted,
      final_training_prediction = final_prediction,
      final_weights = stats::setNames(
        as.numeric(final_fit$weights),
        selected_ids_chr
      ),
      final_intercept = final_fit$intercept,
      coefficients_by_fold = by_fold,
      final_coefficients = final_coefficients
    )
  }

  # ==========================================================================
  # PERFORMANCE METRIC COLLECTION
  # ==========================================================================

  metric_warning_rows <- list()
  metric_warning_index <- 1L

  collect_performance_metrics <- function(
    observed,
    predicted,
    method_name,
    hard_voting = FALSE
  ) {
    captured <- character(0)

    metrics <- withCallingHandlers(
      .compute_performance_metrics(
        observed = observed,
        predicted = predicted,
        family = family_name,
        threshold = threshold
      ),
      warning = function(w) {
        message_text <- conditionMessage(w)

        if (!(hard_voting &&
              grepl(
                "threshold.*predicted probabilities.*0/1",
                message_text,
                ignore.case = TRUE
              ))) {
          captured <<- c(
            captured,
            message_text
          )
        }

        invokeRestart("muffleWarning")
      }
    )

    if (hard_voting) {
      metrics$ROC_AUC <- NA_real_
      metrics$Brier <- NA_real_
      metrics$LogLoss <- NA_real_
    }

    if (length(captured)) {
      for (warning_text in unique(captured)) {
        metric_warning_rows[[
          metric_warning_index
        ]] <<- data.frame(
          method = method_name,
          warning = warning_text,
          stringsAsFactors = FALSE
        )

        metric_warning_index <<-
          metric_warning_index + 1L
      }
    }

    metrics
  }

  # ==========================================================================
  # MODEL-LEVEL ENSEMBLES
  # ==========================================================================

  ensemble_summary <- NULL
  ensemble_predictions <- NULL

  prediction_weights_by_method <- list()
  stack_intercepts <- stats::setNames(
    rep(
      NA_real_,
      length(methods)
    ),
    methods
  )

  stacking_coefficients_by_fold <- data.frame()
  stacking_final_coefficients <- data.frame()

  if (predictions_needed) {
    summary_list <- vector(
      "list",
      length(methods)
    )
    prediction_list <- vector(
      "list",
      length(methods)
    )

    stack_by_fold_list <- list()
    stack_final_list <- list()
    stack_index <- 1L

    for (method_index in seq_along(methods)) {
      method_name <- methods[[method_index]]
      prediction_weights <- NULL
      stack_intercept_value <- NA_real_
      vote_fraction <- rep(
        NA_real_,
        prediction_matrix$n_oof
      )

      if (identical(
        method_name,
        "stacked"
      )) {
        stacked <- compute_stacking(
          pred_mat = prediction_matrix$pred_mat,
          observed = prediction_matrix$obs,
          folds = prediction_matrix$folds
        )

        ensemble_prediction <-
          stacked$cross_fitted_prediction

        prediction_weights <-
          stacked$final_weights

        stack_intercept_value <-
          stacked$final_intercept

        stack_intercepts[[method_name]] <-
          stack_intercept_value

        current_by_fold <-
          stacked$coefficients_by_fold
        current_by_fold$method <-
          method_name

        current_final <-
          stacked$final_coefficients
        current_final$method <-
          method_name

        stack_by_fold_list[[stack_index]] <-
          current_by_fold
        stack_final_list[[stack_index]] <-
          current_final
        stack_index <- stack_index + 1L

      } else if (identical(
        method_name,
        "hard_voting"
      )) {
        class_matrix <-
          prediction_matrix$pred_mat >=
          threshold

        vote_fraction <- rowMeans(
          class_matrix
        )

        highest_ranked_class <- as.numeric(
          class_matrix[, 1L]
        )

        ensemble_prediction <- ifelse(
          vote_fraction > 0.5,
          1,
          ifelse(
            vote_fraction < 0.5,
            0,
            highest_ranked_class
          )
        )

        prediction_weights <- stats::setNames(
          rep(
            NA_real_,
            n_selected
          ),
          selected_ids_chr
        )

      } else {
        prediction_weights <-
          compute_model_weights(
            selected_models,
            method_name
          )

        ensemble_prediction <- as.numeric(
          prediction_matrix$pred_mat %*%
            as.numeric(
              prediction_weights
            )
        )
      }

      if (any(!is.finite(
        ensemble_prediction
      ))) {
        stop(
          "Method '",
          method_name,
          "' produced non-finite ensemble predictions.",
          call. = FALSE
        )
      }

      if (identical(outcome_type, "binary") &&
          !identical(
            method_name,
            "hard_voting"
          )) {
        if (any(
          ensemble_prediction < 0 |
          ensemble_prediction > 1
        )) {
          stop(
            "Method '",
            method_name,
            "' produced values outside [0, 1] for a binomial outcome.",
            call. = FALSE
          )
        }
      }

      prediction_weights_by_method[[
        method_name
      ]] <- prediction_weights

      current_predictions <-
        prediction_matrix$obs_df[
          ,
          id_cols,
          drop = FALSE
        ]

      current_predictions$method <-
        method_name
      current_predictions$observed <-
        prediction_matrix$obs
      current_predictions$predicted_ensemble <-
        ensemble_prediction

      if (identical(outcome_type, "binary")) {
        if (identical(
          method_name,
          "hard_voting"
        )) {
          current_predictions$predicted_probability <-
            NA_real_
          current_predictions$predicted_class <-
            as.integer(
              ensemble_prediction
            )
          current_predictions$vote_fraction <-
            vote_fraction

        } else {
          current_predictions$predicted_probability <-
            ensemble_prediction
          current_predictions$predicted_class <-
            as.integer(
              ensemble_prediction >= threshold
            )
          current_predictions$vote_fraction <-
            NA_real_
        }
      }

      base_prediction_frame <- as.data.frame(
        prediction_matrix$pred_mat,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      current_predictions <- cbind(
        current_predictions,
        base_prediction_frame
      )

      prediction_list[[method_index]] <-
        current_predictions

      metrics <- collect_performance_metrics(
        observed = prediction_matrix$obs,
        predicted = ensemble_prediction,
        method_name = method_name,
        hard_voting = identical(
          method_name,
          "hard_voting"
        )
      )

      summary_list[[method_index]] <- data.frame(
        method = method_name,
        n_models = n_selected,
        n_oof = prediction_matrix$n_oof,
        family = family_name,
        outcome_type = outcome_type,
        cv_method = cv_method,
        cv_scheme = cv_scheme,
        k = effective_k,
        weight_metric = weight_metric,
        metric_direction = metric_direction,
        threshold = if (
          identical(outcome_type, "binary")
        ) {
          threshold
        } else {
          NA_real_
        },
        stacking_model = if (
          identical(method_name, "stacked")
        ) {
          stacking_model
        } else {
          NA_character_
        },
        stack_objective = if (
          identical(method_name, "stacked")
        ) {
          stack_objective
        } else {
          NA_character_
        },
        stack_intercept = if (
          identical(method_name, "stacked")
        ) {
          stack_intercept_value
        } else {
          NA_real_
        },
        weight_transform = if (
          identical(method_name, "weighted")
        ) {
          weight_transform
        } else {
          NA_character_
        },
        metrics,
        stringsAsFactors = FALSE
      )
    }

    ensemble_summary <- do.call(
      rbind,
      summary_list
    )
    rownames(ensemble_summary) <- NULL

    ensemble_predictions <- do.call(
      rbind,
      prediction_list
    )
    rownames(ensemble_predictions) <- NULL

    if (length(stack_by_fold_list)) {
      stacking_coefficients_by_fold <-
        do.call(
          rbind,
          stack_by_fold_list
        )
      rownames(
        stacking_coefficients_by_fold
      ) <- NULL

      stacking_final_coefficients <-
        do.call(
          rbind,
          stack_final_list
        )
      rownames(
        stacking_final_coefficients
      ) <- NULL
    }
  }

  # ==========================================================================
  # AUTOMATIC LAG DECOMPOSITION
  # ==========================================================================

  ensemble_by_lag <- NULL
  lag_weights_by_method <- list()

  if (lag_needed &&
      is.null(lag_data)) {

    if (is.null(fit_list)) {
      fit_list <- attr(
        bestfit,
        "fits",
        exact = TRUE
      )
    }

    if (!is.list(fit_list) ||
        !length(fit_list)) {
      stop(
        "No retained full-data fitted models were found. Provide `fit_list`, ",
        "call `find_bestfit(..., keep_fits = TRUE)`, or supply deterministic ",
        "`lag_data`.",
        call. = FALSE
      )
    }

    if (is.null(names(fit_list)) ||
        anyNA(names(fit_list)) ||
        any(!nzchar(names(fit_list))) ||
        anyDuplicated(names(fit_list))) {
      stop(
        "`fit_list` must have unique non-empty model-ID names.",
        call. = FALSE
      )
    }

    missing_fits <- setdiff(
      selected_ids_chr,
      names(fit_list)
    )

    if (length(missing_fits)) {
      stop(
        "No retained full-data fit was found for selected model(s): ",
        paste(missing_fits, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.data.frame(data) ||
        !nrow(data)) {
      stop(
        "`data` is required for automatic lag decomposition.",
        call. = FALSE
      )
    }

    required_data <- c(
      group,
      "time"
    )

    missing_data <- setdiff(
      required_data,
      names(data)
    )

    if (length(missing_data)) {
      stop(
        "Automatic lag-decomposition `data` is missing column(s): ",
        paste(missing_data, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (anyNA(data[[group]])) {
      stop(
        "The lag-decomposition grouping column cannot contain missing values.",
        call. = FALSE
      )
    }

    if (!is.numeric(data$time) ||
        anyNA(data$time) ||
        any(!is.finite(data$time))) {
      stop(
        "`data$time` must contain only finite numeric values.",
        call. = FALSE
      )
    }

    .epix_validate_regular_time(
      data = data,
      group_index = as.character(
        data[[group]]
      ),
      time_col = "time"
    )

    lag_group_counts <- table(
      as.character(
        data[[group]]
      )
    )

    invalid_lag_groups <- names(lag_group_counts)[
      lag_group_counts != history_length_metadata
    ]

    if (length(invalid_lag_groups)) {
      examples <- paste0(
        invalid_lag_groups,
        "=",
        as.integer(
          lag_group_counts[
            invalid_lag_groups
          ]
        )
      )

      stop(
        "Every group in automatic lag-decomposition `data` must contain ",
        "exactly ", history_length_metadata,
        " observations for max_lag = ", max_lag_metadata,
        ". Non-matching group(s) include: ",
        paste(
          utils::head(
            examples,
            5L
          ),
          collapse = ", "
        ),
        if (length(examples) > 5L) "; ..." else ".",
        " Histories are not truncated, padded, or silently realigned.",
        call. = FALSE
      )
    }

    data_lag <- data[
      order(
        data[[group]],
        data$time
      ),
      ,
      drop = FALSE
    ]
    rownames(data_lag) <- NULL

    lag_list <- list()
    lag_index <- 1L

    for (model_id in selected_ids_chr) {
      fitted <- fit_list[[model_id]]

      fit_meta <- .get_epiexposure_metadata(
        fitted
      )

      if (!identical(
        fit_meta$family_name,
        family_name
      )) {
        stop(
          "Retained fit ",
          model_id,
          " has family '",
          fit_meta$family_name,
          "', inconsistent with `bestfit` family '",
          family_name,
          "'.",
          call. = FALSE
        )
      }

      if (!identical(
        fit_meta$prediction_level,
        "population"
      ) ||
      !identical(
        fit_meta$prediction_estimand,
        "expected_response"
      ) ||
      !identical(
        fit_meta$point_prediction_contract,
        "central_expected_response"
      )) {
        stop(
          "Retained fit ",
          model_id,
          " does not satisfy the population central expected-response contract.",
          call. = FALSE
        )
      }

      if (!identical(
        fit_meta$max_lag,
        max_lag_metadata
      ) ||
      !identical(
        fit_meta$history_length,
        history_length_metadata
      ) ||
      !identical(
        fit_meta$history_contract,
        history_contract_metadata
      )) {
        stop(
          "Retained fit ",
          model_id,
          " does not match the `bestfit` exact-history metadata. All selected ",
          "models must use the same max_lag and exactly max_lag + 1 time points.",
          call. = FALSE
        )
      }

      fitted_vars <- fit_meta$vars

      if (!length(fitted_vars)) {
        stop(
          "Retained fit ",
          model_id,
          " contains no fitted exposure metadata.",
          call. = FALSE
        )
      }

      if (identical(
        lg_strategy,
        "model"
      )) {
        selected_vars <- fitted_vars

      } else if (identical(
        lg_strategy,
        "requested_available"
      )) {
        selected_vars <- if (is.null(var)) {
          fitted_vars
        } else {
          intersect(
            var,
            fitted_vars
          )
        }

        if (!length(selected_vars)) {
          if (verbose) {
            warning(
              "Skipping lag decomposition for model ",
              model_id,
              ": none of the requested exposures are fitted in that model.",
              call. = FALSE
            )
          }
          next
        }

      } else {
        selected_vars <- if (is.null(var)) {
          fitted_vars
        } else {
          var
        }

        missing_vars <- setdiff(
          selected_vars,
          fitted_vars
        )

        if (length(missing_vars)) {
          stop(
            "Model ",
            model_id,
            " lacks requested exposure(s): ",
            paste(missing_vars, collapse = ", "),
            ".",
            call. = FALSE
          )
        }
      }

      missing_exposure_data <- setdiff(
        selected_vars,
        names(data_lag)
      )

      if (length(missing_exposure_data)) {
        stop(
          "`data` is missing exposure(s) required for lag decomposition of ",
          "model ",
          model_id,
          ": ",
          paste(missing_exposure_data, collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      for (variable in selected_vars) {
        if (!is.numeric(data_lag[[variable]]) ||
            anyNA(data_lag[[variable]]) ||
            any(!is.finite(data_lag[[variable]]))) {
          stop(
            "Exposure variable '", variable,
            "' must contain only finite numeric values for automatic lag ",
            "decomposition. All selected exposures must share the same validated ",
            "time rows.",
            call. = FALSE
          )
        }
      }

      call_args <- c(
        list(
          data = data_lag,
          group = group,
          fit = fitted,
          var = selected_vars
        ),
        compute_ecilag_args
      )

      decomposition <- tryCatch(
        do.call(
          compute_ecilag,
          call_args
        ),
        error = function(e) {
          stop(
            "Lag decomposition failed for model ",
            model_id,
            ": ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )

      if (!is.list(decomposition) ||
          !is.data.frame(
            decomposition$by_lag
          ) ||
          !nrow(
            decomposition$by_lag
          )) {
        stop(
          "`compute_ecilag()` returned no valid non-empty `by_lag` for model ",
          model_id,
          ".",
          call. = FALSE
        )
      }

      temporary <- decomposition$by_lag

      selected_position <- match(
        model_id,
        selected_ids_chr
      )

      temporary[[model_col]] <- rep(
        selected_ids[
          selected_position
        ],
        nrow(temporary)
      )

      lag_list[[lag_index]] <-
        temporary
      lag_index <- lag_index + 1L
    }

    if (!length(lag_list)) {
      stop(
        "No lag decompositions were generated for the selected models.",
        call. = FALSE
      )
    }

    lag_data <- do.call(
      rbind,
      lag_list
    )
    rownames(lag_data) <- NULL

    if (verbose) {
      message(
        "Deterministic lag contributions generated from retained full-data fits."
      )
    }
  }

  # ==========================================================================
  # LAG ENSEMBLES
  # ==========================================================================

  make_group_key <- function(
    d,
    columns
  ) {
    make_key(
      d,
      columns
    )
  }

  if (lag_needed) {
    if (!is.data.frame(lag_data) ||
        !nrow(lag_data)) {
      stop(
        "`lag_data` must be a non-empty data.frame.",
        call. = FALSE
      )
    }

    required_lag <- unique(
      c(
        model_col,
        lag_group_cols,
        lag_contribution_col
      )
    )

    missing_lag <- setdiff(
      required_lag,
      names(lag_data)
    )

    if (length(missing_lag)) {
      stop(
        "Missing lag-data column(s): ",
        paste(missing_lag, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.numeric(
      lag_data[[lag_contribution_col]]
    ) ||
    anyNA(
      lag_data[[lag_contribution_col]]
    ) ||
    any(!is.finite(
      lag_data[[lag_contribution_col]]
    ))) {
      stop(
        "Lag contributions must contain only finite numeric values.",
        call. = FALSE
      )
    }

    if (!is.numeric(lag_data[[lag_col]]) ||
        anyNA(lag_data[[lag_col]]) ||
        any(!is.finite(lag_data[[lag_col]])) ||
        any(lag_data[[lag_col]] != as.integer(lag_data[[lag_col]])) ||
        any(lag_data[[lag_col]] < 0L |
            lag_data[[lag_col]] > max_lag_metadata)) {
      stop(
        "`lag_data$lag` must contain integer lags from 0 through the inherited ",
        "common max_lag = ", max_lag_metadata, ".",
        call. = FALSE
      )
    }

    lag_use <- lag_data[
      as.character(
        lag_data[[model_col]]
      ) %in% selected_ids_chr,
      ,
      drop = FALSE
    ]

    if (!nrow(lag_use)) {
      stop(
        "No lag contributions were found for the selected models.",
        call. = FALSE
      )
    }

    lag_models_present <- unique(
      as.character(
        lag_use[[model_col]]
      )
    )

    if (identical(
      lg_strategy,
      "strict"
    )) {
      missing_lag_models <- setdiff(
        selected_ids_chr,
        lag_models_present
      )

      if (length(missing_lag_models)) {
        stop(
          "`lg_strategy = 'strict'` requires lag contributions from every ",
          "selected model. Missing model(s): ",
          paste(missing_lag_models, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    }

    lag_key_cols <- c(
      model_col,
      lag_group_cols
    )

    if (anyDuplicated(
      lag_use[lag_key_cols]
    )) {
      stop(
        "`lag_data` has multiple rows for the same model and lag unit. ",
        "Add the relevant grouping column to `lag_group_cols` or aggregate ",
        "before calling `ensemble_bestfit()`.",
        call. = FALSE
      )
    }

    lag_unit_key <- make_group_key(
      lag_use,
      lag_group_cols
    )

    if (identical(
      lg_strategy,
      "strict"
    )) {
      unit_model_counts <- vapply(
        split(
          as.character(
            lag_use[[model_col]]
          ),
          lag_unit_key
        ),
        function(x) length(unique(x)),
        integer(1)
      )

      if (any(
        unit_model_counts != n_selected
      )) {
        stop(
          "`lg_strategy = 'strict'` requires every lag unit to contain all ",
          n_selected,
          " selected models.",
          call. = FALSE
        )
      }
    }

    lag_method_rows <- vector(
      "list",
      length(methods)
    )

    for (method_index in seq_along(methods)) {
      method_name <- methods[[method_index]]

      lag_weights <- compute_model_weights(
        selected_models,
        method_name
      )

      lag_weights_by_method[[
        method_name
      ]] <- lag_weights

      current_lag_use <- lag_use
      current_lag_use$.model_weight <-
        as.numeric(
          lag_weights[
            match(
              as.character(
                current_lag_use[[model_col]]
              ),
              names(lag_weights)
            )
          ]
        )

      if (any(!is.finite(
        current_lag_use$.model_weight
      ))) {
        stop(
          "Could not match model weights to `lag_data`.",
          call. = FALSE
        )
      }

      current_group_key <- make_group_key(
        current_lag_use,
        lag_group_cols
      )

      ensemble_rows <- lapply(
        unique(current_group_key),
        function(current_key) {
          temporary <- current_lag_use[
            current_group_key == current_key,
            ,
            drop = FALSE
          ]

          raw_weights <- temporary$.model_weight
          raw_weight_sum <- sum(raw_weights)

          if (!is.finite(raw_weight_sum) ||
              raw_weight_sum < 0) {
            stop(
              "Invalid model-weight sum within a lag unit.",
              call. = FALSE
            )
          }

          weights <- raw_weights

          if (rl_weights) {
            if (raw_weight_sum <= 0) {
              stop(
                "Lag-unit model weights sum to zero and cannot be ",
                "renormalized.",
                call. = FALSE
              )
            }

            weights <- weights /
              raw_weight_sum
          }

          output <- temporary[
            1L,
            lag_group_cols,
            drop = FALSE
          ]

          output$method <- method_name
          output$ECI_weighted_ens <-
            sum(
              temporary[[
                lag_contribution_col
              ]] * weights
            )

          if (lag_weight_col %in%
              names(temporary)) {
            if (!is.numeric(
              temporary[[lag_weight_col]]
            ) ||
            anyNA(
              temporary[[lag_weight_col]]
            ) ||
            any(!is.finite(
              temporary[[lag_weight_col]]
            ))) {
              stop(
                "Optional lag `weight` values must be finite numeric values.",
                call. = FALSE
              )
            }

            output$ensemble_weight <- sum(
              temporary[[
                lag_weight_col
              ]] * weights
            )
          }

          output$n_models <-
            length(
              unique(
                temporary[[model_col]]
              )
            )
          output$n_selected_models <-
            n_selected
          output$raw_model_weight_sum <-
            raw_weight_sum
          output$weights_renormalized <-
            rl_weights
          output$models_used <- paste(
            sort(
              unique(
                as.character(
                  temporary[[model_col]]
                )
              )
            ),
            collapse = ", "
          )

          output
        }
      )

      method_lag <- do.call(
        rbind,
        ensemble_rows
      )
      rownames(method_lag) <- NULL

      percent_groups <- setdiff(
        lag_group_cols,
        lag_col
      )

      method_lag$ECI_percent_ens <-
        NA_real_

      if (!length(percent_groups)) {
        denominator <- sum(
          abs(
            method_lag$ECI_weighted_ens
          )
        )

        if (denominator > 0) {
          method_lag$ECI_percent_ens <-
            100 *
            abs(
              method_lag$ECI_weighted_ens
            ) /
            denominator
        }

      } else {
        percent_key <- make_group_key(
          method_lag,
          percent_groups
        )

        for (current_key in unique(
          percent_key
        )) {
          rows <- which(
            percent_key == current_key
          )

          denominator <- sum(
            abs(
              method_lag$ECI_weighted_ens[
                rows
              ]
            )
          )

          if (denominator > 0) {
            method_lag$ECI_percent_ens[
              rows
            ] <-
              100 *
              abs(
                method_lag$ECI_weighted_ens[
                  rows
                ]
              ) /
              denominator
          }
        }
      }

      lag_method_rows[[method_index]] <-
        method_lag
    }

    ensemble_by_lag <- do.call(
      rbind,
      lag_method_rows
    )
    rownames(ensemble_by_lag) <- NULL
  }

  # ==========================================================================
  # MODEL-WEIGHT TABLE
  # ==========================================================================

  model_weight_rows <- vector(
    "list",
    length(methods)
  )

  for (method_index in seq_along(methods)) {
    method_name <- methods[[method_index]]

    prediction_weights <-
      prediction_weights_by_method[[
        method_name
      ]]

    lag_weights <-
      lag_weights_by_method[[
        method_name
      ]]

    current <- data.frame(
      method = rep(
        method_name,
        n_selected
      ),
      selected_ids,
      metric_value = as.numeric(
        selected_models[[
          weight_metric
        ]]
      ),
      prediction_weight = if (
        !is.null(prediction_weights)
      ) {
        as.numeric(
          prediction_weights[
            match(
              selected_ids_chr,
              names(
                prediction_weights
              )
            )
          ]
        )
      } else {
        rep(
          NA_real_,
          n_selected
        )
      },
      lag_weight = if (
        !is.null(lag_weights)
      ) {
        as.numeric(
          lag_weights[
            match(
              selected_ids_chr,
              names(
                lag_weights
              )
            )
          ]
        )
      } else {
        rep(
          NA_real_,
          n_selected
        )
      },
      stringsAsFactors = FALSE
    )

    names(current)[2L] <- model_col

    model_weight_rows[[method_index]] <-
      current
  }

  model_weights <- do.call(
    rbind,
    model_weight_rows
  )
  rownames(model_weights) <- NULL

  if (length(methods) == 1L &&
      identical(methods, "stacked")) {
    attr(
      model_weights,
      "stack_intercept"
    ) <- stack_intercepts[["stacked"]]
  }

  # ==========================================================================
  # METRIC WARNINGS
  # ==========================================================================

  metric_warnings <- if (
    length(metric_warning_rows)
  ) {
    output <- do.call(
      rbind,
      metric_warning_rows
    )
    rownames(output) <- NULL
    unique(output)
  } else {
    data.frame(
      method = character(0),
      warning = character(0),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(metric_warnings)) {
    warning(
      "One or more ensemble performance metrics were undefined or required ",
      "special handling. See `result$metric_warnings`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # RETURN
  # ==========================================================================

  out <- list(
    ensemble_summary = ensemble_summary,
    ensemble_predictions = ensemble_predictions,
    ensemble_by_lag = ensemble_by_lag,
    lag_data = lag_data,
    model_weights = model_weights,
    selected_models = selected_models,
    stacking_coefficients_by_fold =
      stacking_coefficients_by_fold,
    stacking_final_coefficients =
      stacking_final_coefficients,
    metric_warnings = metric_warnings,
    method = methods,
    ensemble_scope = ensemble_scope,
    lg_strategy = lg_strategy,
    rl_weights = rl_weights,
    family = family_name,
    outcome_type = outcome_type,
    weight_metric = weight_metric,
    metric_direction = metric_direction,
    threshold = if (
      identical(outcome_type, "binary")
    ) {
      threshold
    } else {
      NA_real_
    },
    prediction_level = prediction_level,
    prediction_estimand = prediction_estimand,
    prediction_contract = prediction_contract,
    cv_method = cv_method,
    cv_scheme = cv_scheme,
    k = effective_k,
    cv_seed = cv_seed,
    cv_stratified = cv_stratified,
    fold_assignments = fold_assignments,
    fold_balance = fold_balance,
    basis_training_only = TRUE,
    max_lag = max_lag_metadata,
    history_length = history_length_metadata,
    history_contract = history_contract_metadata,
    time_step = time_step_metadata
  )

  out
}
