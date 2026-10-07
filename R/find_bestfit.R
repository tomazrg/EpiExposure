#' Find the best DLNM model structure using grouped validation
#'
#' Evaluates combinations of exposure variables and cross-basis dimensions using
#' grouped validation, then ranks candidate models with family-aware predictive
#' performance metrics. Three grouped validation strategies are available:
#'
#' - `"LOOCV"`: grouped leave-one-out cross-validation;
#' - `"k-fold"`: grouped k-fold cross-validation;
#' - `"holdout"`: one grouped proportional holdout validation split.
#'
#' The validation unit is always the complete column named by `group`. The
#' function creates the internal canonical alias `epi_id` from that column, but
#' never splits a group's exposure-history rows across training and evaluation.
#' With the default `group = "epi_id"`, complete epidemics are therefore kept
#' intact in every strategy.
#'
#' The statistical target used for every validation prediction is the same target
#' used by the EpiExposure v1 prediction layer: the population/fixed-component
#' expected response, with fitted random effects set to zero. Model ranking uses
#' deterministic predictions from the harmonized central parameter estimate and
#' does not propagate coefficient, posterior, residual, or future-observation
#' uncertainty.
#'
#' @param data Long-format data frame containing the response, grouping column,
#'   chronological time column, and candidate exposure variables.
#' @param response Character scalar naming the response column in `data`. The
#'   response may be repeated over exposure-history rows but must be constant
#'   within each complete validation group because one outcome is predicted per
#'   group.
#' @param group Character scalar naming the independent grouped-validation unit.
#'   Every row belonging to one group remains together. Internally,
#'   `data_long$epi_id` is a canonical alias of `data_long[[group]]`; it does not
#'   redefine the user's grouping structure.
#' @param time Character scalar naming the chronological time column. Time must
#'   be finite and unique within groups. Within every group observations must be
#'   complete and equally spaced, and the time step must be the same across
#'   groups. Every group must contain exactly `max_lag + 1` time points. Because
#'   all candidate exposure variables are finite columns of the same validated
#'   long-format rows, all variables necessarily use the same temporal support.
#'   Rows are ordered internally by `group` and `time`.
#' @param vars Unique character vector of candidate exposure variables. Exposure
#'   variables must be distinct from `response`, `group`, and `time`.
#' @param max_lag Non-negative integer maximum lag, expressed in observation
#'   steps. A two-element form such as `c(0, 85)` is accepted and normalized to
#'   its maximum. Because all groups must use the same time step, a lag index has
#'   the same temporal meaning across groups.
#' @param df_var_grid Positive integer candidate dimensions for the
#'   exposure-response basis.
#' @param df_lag_grid Positive integer candidate dimensions for the lag-response
#'   basis.
#' @param min_vars Positive integer minimum number of exposure variables.
#' @param max_vars Positive integer maximum number of exposure variables, or
#'   `NULL` to use `length(vars)`.
#' @param var_sets Optional list of exact variable combinations. When supplied,
#'   `min_vars` and `max_vars` do not determine the combinations, although their
#'   basic argument validation is still applied.
#' @param fun_var Character exposure-basis function passed to
#'   `define_exposures()`.
#' @param fun_lag Character lag-basis function passed to `define_exposures()`.
#' @param model_engine Character model engine supported by `fit_epidlnm()`:
#'   `"glm"`, `"glmmTMB"`, `"gam"`, `"gamm"`, `"gls"`, `"spamm"`,
#'   `"brms"`, `"inla"`, or `"bdlnm"`.
#' @param family Character or family object passed to `prepare_response()` and
#'   `fit_epidlnm()`. Supported EpiExposure v1 canonical families are `"beta"`,
#'   `"binomial"`, `"poisson"`, `"gamma"`, `"gaussian"`, and
#'   `"negative_binomial"` (NB2). Ordinal outcomes and NB1 are not supported in
#'   EpiExposure v1. Family-object links are checked against the selected engine
#'   before validation and passed unchanged to `fit_epidlnm()`. To specify a
#'   non-default link, supply a supported family object, such as
#'   `family = stats::Gamma(link = "inverse")`; character family names use the
#'   EpiExposure default link.
#' @param random_effect Optional character scalar naming one grouping column used
#'   as a conventional random intercept by `fit_epidlnm()`. It must be constant
#'   within each validation group. Random slopes, nested or crossed
#'   random-effect specifications, and arbitrary engine-specific random-effect
#'   expressions are not supported.
#'
#'   For `brms`, the random intercept is fitted as `(1 | group)`. For the
#'   INLA-backed `inla` and `bdlnm` engines, it is fitted as
#'   `INLA::f(group, model = "iid")`. EpiExposure v1 does not expose structured
#'   INLA latent effects such as `"rw1"`, `"rw2"`, `"ar1"`, `"besag"`, `"bym"`,
#'   `"bym2"`, or SPDE effects.
#'
#'   Random effects may contribute to model fitting, but validation predictions
#'   are always population-level and therefore set fitted random effects to
#'   zero. `glm` and `gls` do not accept `random_effect`; `gamm` requires it
#'   under the current EpiExposure v1 fitting contract.
#' @param random_effect_prior Optional named list defining the hyperprior for the
#'   precision of the IID random intercept when `model_engine = "inla"` or
#'   `"bdlnm"`. The list is passed unchanged to the `hyper` argument of the
#'   internally generated `INLA::f(..., model = "iid")` term. For example,
#'   `list(prec = list(prior = "pc.prec", param = c(1, 0.01)))`.
#'
#'   `NULL` retains the engine default. A non-`NULL` value requires
#'   `random_effect` and is currently supported only for INLA-backed engines.
#'   The argument customizes the hyperprior of the IID structure; it does not
#'   select another latent model. The same prior is used in every validation
#'   training fit and in retained full-data refits.
#' @param spatial_effect `NULL` (default), or two distinct numeric coordinate
#'   column names for a Matérn term fitted by `model_engine = "spamm"` only.
#'   Coordinates are constant within each epidemic but may be identical across
#'   different epidemics; they never redefine the exposure-history unit.
#' @param spatial_structure Only `"matern"` is supported (case-insensitive)
#'   when spatial coordinates are supplied; otherwise this setting is inert.
#' @param spatial_group `NULL` for one shared field or a factor, character, or
#'   integer grouping column such as `"year"` for independent spatial fields.
#'   Requires spatial coordinates and must be constant within each epidemic.
#'   Spatial effects and `random_effect` are distinct and may coexist.
#' @param min_success Positive integer of at least 2 giving the minimum number of
#'   **evaluation groups with a finite out-of-sample validation prediction**
#'   required for a candidate to be eligible for metric calculation and, when
#'   requested, full-data refitting. Under LOOCV and k-fold, all groups are
#'   evaluation groups. Under holdout, only the fixed test groups count toward
#'   `n_success`, `n_failed`, and `min_success`; holdout training groups are never
#'   counted as failed predictions. Candidate diagnostics retain
#'   `n_success_folds` and `n_failed_folds` as complete-validation-iteration
#'   counters: for holdout there is one iteration, successful only when every
#'   designated test group receives a finite prediction.
#' @param rank_metric Optional character scalar selecting the primary ranking
#'   metric. If `NULL`, the default is `"CCC"` for non-binary outcomes and
#'   `"ROC_AUC"` for binary outcomes. Non-binary metrics are `"CCC"`, `"Cb"`,
#'   `"rho"`, `"RMSE"`, and `"MAE"`. Binary metrics are `"ROC_AUC"`, `"Brier"`,
#'   `"LogLoss"`, `"Accuracy"`, `"Balanced_Accuracy"`, `"Sensitivity"`,
#'   `"Specificity"`, `"F1"`, `"MCC"`, and `"Precision"`. Metric direction is
#'   resolved automatically.
#' @param threshold Numeric probability strictly between 0 and 1 used only for
#'   binomial classification metrics. The default is `0.5`. Predicted
#'   probabilities are retained; thresholding is used only to create the
#'   auxiliary 0/1 class required by classification metrics. ROC AUC, Brier
#'   score, and Log Loss do not depend on the threshold.
#' @param top_n Positive integer or `Inf`; number of ranked candidates returned.
#' @param keep_fits Logical. Validation metrics are always calculated first and
#'   exclusively from out-of-sample predictions produced by the selected
#'   validation strategy. If `TRUE`, every candidate that reaches `min_success`
#'   is **then** refitted on the complete data using a full-data basis. These
#'   full-data refits do not generate, replace, or alter validation predictions
#'   and do not participate in metric calculation or ranking. After ranking and
#'   `top_n` filtering, only retained candidate fits are stored in
#'   `attr(result, "fits")`. They are intended for downstream visualization,
#'   epidemiological interpretation, exposure-lag-response surfaces, effect
#'   summaries, lag decomposition, scenario simulation, and later prediction.
#' @param verbose Logical. Print the validation strategy, candidate grid,
#'   execution mode, concise warning summaries, and completion information.
#' @param validation_method Character grouped-validation strategy. `"LOOCV"`
#'   leaves one complete group out per iteration. `"k-fold"` assigns complete
#'   groups to `k` approximately balanced folds. `"holdout"` creates one fixed
#'   grouped proportional training/test split and uses the test groups as the
#'   validation set for every candidate. All strategies operate on complete
#'   groups and never split temporal rows from one group across partitions.
#' @param k Positive integer number of folds used only when
#'   `validation_method = "k-fold"`. It must satisfy
#'   `2 <= k < number of groups`. Fold sizes differ by at most one group. For
#'   binomial outcomes, allocation is additionally stratified by the group-level
#'   0/1 response and each response class must contain at least `k` groups.
#' @param test_prop Numeric scalar strictly between 0 and 1 giving the requested
#'   approximate proportion of complete groups reserved for the test partition
#'   when `validation_method = "holdout"`. The default `0.30` requests
#'   approximately 70 percent of groups for training and 30 percent for testing.
#'   It is validated for every call but ignored by LOOCV and k-fold. For binary
#'   outcomes, the effective holdout size may be adjusted to the closest feasible
#'   value that keeps both outcomes 0 and 1 in both training and test partitions.
#' @param seed Optional non-negative integer controlling grouped allocation for
#'   `validation_method = "k-fold"` and the grouped train/test split for
#'   `validation_method = "holdout"`. The same data, group order, validation
#'   method, `test_prop`, and seed reproduce the same allocation. The caller's
#'   existing global `.Random.seed` is restored after allocation. `seed` does not
#'   alter LOOCV membership and is not passed to `fit_epidlnm()`.
#' @param ... Named additional arguments passed to `fit_epidlnm()`. Core
#'   arguments managed by `find_bestfit()` (`data`, `model_engine`, `family`,
#'   `random_effect`, `random_effect_prior`, `spatial_effect`,
#'   `spatial_structure`, `spatial_group`, `epiexposure_spec`, and
#'   `basis_objects`) cannot be supplied again through `...`. For INLA-backed
#'   engines, likelihood availability, `control.family$control.link$model`
#'   consistency, and `control.compute$config = TRUE` compatibility are checked
#'   before validation. `fit_epidlnm()` supplies required INLA controls during
#'   each fit.
#'
#' @return A data frame ranked by `rank_metric`. Non-binary families report
#'   `CCC`, `Cb`, `rho`, `RMSE`, and `MAE`; binomial models report `ROC_AUC`,
#'   `Brier`, `LogLoss`, `Accuracy`, `Balanced_Accuracy`, `Sensitivity`,
#'   `Specificity`, `F1`, `MCC`, and `Precision`.
#'
#'   Candidate diagnostics include `n_groups_total`, `n_predictions`,
#'   `n_success`, `n_failed`, `n_folds`, `n_success_folds`, and
#'   `n_failed_folds`. `n_predictions` is the number of groups assigned to
#'   out-of-sample evaluation: all groups for LOOCV/k-fold, but only test groups
#'   for holdout. The `*_folds` fields are retained as validation-iteration
#'   counters; holdout therefore has one iteration.
#'
#'   Diagnostic attributes include `"predictions"`, `"predictions_by_model"`,
#'   `"failures"`, `"warnings"`, `"warning_summary"`, `"timing"`, and
#'   `"candidate_times"`. With `keep_fits = TRUE`, `"fits"` contains retained
#'   full-data refits and `"retained_fit_scope"` is
#'   `"full_data_refit_after_validation"`; otherwise `"retained_fit_scope"` is
#'   `NA_character_`.
#'
#'   Grouped-validation metadata include:
#'
#'   - `"validation_method"`: `"LOOCV"`, `"k-fold"`, or `"holdout"`;
#'   - `"validation_scheme"`: normalized descriptive scheme;
#'   - `"k"`: number of LOOCV iterations for LOOCV, supplied `k` for k-fold,
#'     and `NA_integer_` for holdout;
#'   - `"validation_seed"`: supplied k-fold/holdout seed or `NA_integer_`;
#'   - `"validation_stratified"`: whether binomial stratification was used;
#'   - `"test_prop"`: requested holdout proportion, otherwise `NA_real_`;
#'   - `"n_groups_total"` and `"n_evaluation_groups"`;
#'   - `"training_groups"`, `"test_groups"`, and `"evaluation_groups"`;
#'   - `"validation_assignments"`: data frame with `group`, `partition`, `fold`,
#'     and `observed`; holdout uses fixed `train`/`test` partitions, whereas
#'     LOOCV/k-fold use `partition = "evaluation"` and the evaluation fold;
#'   - `"validation_balance"`: fold-level counts for LOOCV/k-fold or train/test
#'     counts for holdout, including 0/1 class counts for binomial outcomes.
#'
#'   Standard downstream metadata are retained in `"family"`, `"outcome_type"`,
#'   `"rank_metric"`, `"threshold"`, `"prediction_level"`,
#'   `"prediction_estimand"`, `"prediction_contract"`, `"random_effect"`,
#'   `"random_effect_model"`, `"random_effect_prior"`, `"spatial_effect"`,
#'   `"spatial_structure"`, `"spatial_group"`, `"spatial_term"`,
#'   `"has_spatial_effect"`, `"basis_training_only"`, `"max_lag"`,
#'   `"history_length"`, `"history_contract"`, and `"time_step"`.
#'
#' @details
#' ## Grouped validation and cross-validation
#'
#' Validation membership is constructed once before candidate evaluation and is
#' reused unchanged for every candidate. This guarantees fair comparison on the
#' same out-of-sample groups. LOOCV produces one evaluation iteration per group.
#' Grouped k-fold randomizes complete groups into approximately balanced folds;
#' for binomial outcomes, allocation is stratified so both classes occur in every
#' fold whenever the documented class-count requirement is satisfied.
#'
#' ## Grouped proportional holdout
#'
#' With `validation_method = "holdout"`, the function creates exactly one grouped
#' training/test split before candidate evaluation. `test_prop` determines the
#' requested approximate fraction of complete groups assigned to testing. For
#' non-binary outcomes, test groups are sampled directly. For binary outcomes,
#' sampling is stratified at the **group-level response**, not at longitudinal
#' row level, and both outcomes 0 and 1 are required in both partitions. At least
#' two complete groups in each class are therefore required.
#'
#' Every candidate uses the same fixed holdout partition. `define_exposures()` is
#' called only on the holdout training groups, so exposure-basis knots, boundary
#' knots, effective dimensions, `argvar`, `arglag`, templates, and related
#' cross-basis attributes are learned exclusively from training data. The stored
#' training template is then transported unchanged to the holdout test histories
#' by `.build_design_from_templates()`. Test values never redefine the basis.
#'
#' Holdout metrics are calculated exclusively from predictions for the fixed test
#' groups. Because that same test subset is used to compare and rank candidate
#' structures, it is a **validation** holdout, not an untouched final external
#' test set. A final independent external evaluation requires another data set or
#' partition that was never used for candidate selection. If `keep_fits = TRUE`,
#' eligible candidates are refitted on all available data only after their
#' validation metrics have already been calculated; these refits are downstream
#' models and never replace holdout predictions.
#'
#' ## Training-only basis construction
#'
#' For every validation iteration, `define_exposures()` receives only the
#' corresponding training groups. The effective `argvar`, `arglag`, lag range,
#' spline knots, boundary knots, and other returned basis attributes are reused
#' to transform the evaluation histories. For highly skewed predictors or many
#' repeated values, `splines::ns()` may warn when an interior knot coincides with
#' a boundary knot; the effective training basis returned after that adjustment
#' is the template transported to evaluation data.
#'
#' The function validates basis transport explicitly. Requested `max_lag`, the
#' lag stored in training cross-bases, and the lag stored in the exposure
#' specification must agree. Reconstructed training and evaluation cross-bases
#' must have matching dimensions and native column names/order. Canonical
#' `cb_<variable>_<index>` names are assigned deterministically. For
#' `model_engine = "bdlnm"`, an epidemic-level matrix-form cross-basis is built
#' from the same training parameters, numerically checked against the canonical
#' design, and exposure-prefixed so internal INLA names remain globally unique.
#'
#' ## Temporal requirements
#'
#' Every complete group must contain exactly `max_lag + 1` ordered, equally
#' spaced observations. Histories are never truncated, padded, or silently
#' realigned. All candidate exposures are finite columns of the same validated
#' long-format rows and therefore share identical temporal support.
#'
#' ## Random-effect scope
#'
#' Candidate models may include one conventional random intercept through
#' `random_effect`. The same random-intercept specification is used in every
#' validation training fit and, when `keep_fits = TRUE`, in the full-data
#' refits retained after validation.
#'
#' For Bayesian engines, EpiExposure v1 supports `(1 | group)` through `brms`
#' and `f(group, model = "iid")` through `inla` and `bdlnm`. The candidate grid
#' does not search over alternative latent-effect structures, and
#' `find_bestfit()` does not support INLA models such as `"rw1"`, `"rw2"`,
#' `"ar1"`, `"besag"`, `"bym"`, `"bym2"`, or SPDE effects.
#'
#' This restriction keeps the validation target and population-level
#' prediction contract harmonized across candidates and engines. All fitted
#' random-effect contributions are set to zero when out-of-sample validation
#' predictions are generated.
#'
#' ## Optional spatial covariance
#'
#' Only `spaMM` accepts `spatial_effect = c("x_coord", "y_coord")` through the
#' EpiExposure v1 interface. This adds a Matérn term during fitting;
#' `spatial_group` may define independent fields. This spaMM-specific interface
#' is distinct from INLA/SPDE or areal latent-effect models, which are not
#' currently exposed by `find_bestfit()`.
#'
#' Conventional random effects remain distinct. All fitted random and spatial
#' effects are set to zero for population-level validation predictions.
#'
#' ## Validation prediction target, warnings, and ranking
#'
#' Every successful evaluation group receives one deterministic out-of-sample
#' expected-response prediction from the harmonized central fixed/population
#' parameter estimate. For LOOCV and k-fold these are out-of-fold predictions;
#' for holdout they are out-of-sample holdout predictions. Coefficient/posterior
#' uncertainty and residual/future-observation noise are excluded from ranking.
#'
#' Warnings raised during fitting/prediction are captured and classified as
#' `"knot"`, `"convergence"`, `"hessian"`, or `"other"`; warnings do not by
#' themselves make a validation iteration fail. Candidate ranking uses
#' `rank_metric` first, with the remaining family-appropriate metrics as
#' deterministic tie-breakers in canonical order.
#'
#' Candidate-level parallelism is preserved. When the caller configures a
#' `future` plan with more than one worker, one future owns each candidate and
#' that candidate's validation iterations remain sequential. The function does
#' not change the caller's global `future::plan()`.
#'
#' @examples
#' \dontrun{
#' # Existing grouped-validation examples:
#' find_bestfit(dat, vars = "tmean", max_lag = 10, model_engine = "spamm",
#'              family = "poisson", random_effect = "epi_id")
#' find_bestfit(dat, vars = "tmean", max_lag = 10, model_engine = "spamm",
#'              family = "poisson", spatial_effect = c("x_coord", "y_coord"))
#' find_bestfit(dat, vars = "tmean", max_lag = 10, model_engine = "spamm",
#'              family = "poisson", random_effect = "block_id",
#'              spatial_effect = c("x_coord", "y_coord"))
#' find_bestfit(dat, vars = "tmean", max_lag = 10, model_engine = "spamm",
#'              family = "poisson", random_effect = "block_id",
#'              spatial_effect = c("x_coord", "y_coord"), spatial_group = "year")
#'
#' # Grouped proportional holdout for a non-binary response:
#' best_holdout <- find_bestfit(
#'   data = dat,
#'   response = "y",
#'   group = "epi_id",
#'   time = "time",
#'   vars = c("tmean", "rain", "wetness"),
#'   max_lag = 10,
#'   model_engine = "glm",
#'   family = "poisson",
#'   validation_method = "holdout",
#'   test_prop = 0.30,
#'   seed = 123,
#'   top_n = 1,
#'   keep_fits = TRUE
#' )
#'
#' # Grouped stratified holdout for a binary response:
#' best_binary_holdout <- find_bestfit(
#'   data = dat_binomial,
#'   response = "y",
#'   group = "epi_id",
#'   time = "time",
#'   vars = c("tmean", "rain"),
#'   max_lag = 10,
#'   model_engine = "glm",
#'   family = "binomial",
#'   validation_method = "holdout",
#'   test_prop = 0.30,
#'   seed = 123,
#'   top_n = 1,
#'   keep_fits = TRUE
#' )
#' }
#' @export
find_bestfit <- function(
    data,
    response = "y",
    group = "epi_id",
    time = "time",
    vars,
    max_lag,
    df_var_grid = c(3, 4, 5),
    df_lag_grid = c(3, 4, 5),
    min_vars = 1,
    max_vars = NULL,
    var_sets = NULL,
    fun_var = "ns",
    fun_lag = "ns",
    model_engine = "glmmTMB",
    family = "beta",
    random_effect = NULL,
    random_effect_prior = NULL,
    spatial_effect = NULL,
    spatial_structure = "matern",
    spatial_group = NULL,
    min_success = 2,
    rank_metric = NULL,
    threshold = 0.5,
    top_n = Inf,
    keep_fits = FALSE,
    verbose = TRUE,
    validation_method = c("LOOCV", "k-fold", "holdout"),
    k = 5,
    test_prop = 0.30,
    seed = NULL,
    ...
) {

  model_engine <- match.arg(
    model_engine,
    choices = c(
      "glm", "glmmTMB", "gam", "gamm", "gls",
      "spamm", "brms", "inla", "bdlnm"
    )
  )

  validation_method <- match.arg(
    validation_method,
    choices = c("LOOCV", "k-fold", "holdout")
  )

  function_start_time <- proc.time()[["elapsed"]]

  # Force `...` in the calling session before futures are created.
  fit_dots <- list(...)

  # Keep preflight family and link validation consistent with fit_epidlnm().
  # Native family constructors and engine validation remain in fit_epidlnm().
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  stopf <- function(...) stop(..., call. = FALSE)

  is_scalar_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }

  normalize_token <- function(x) {
    x <- tolower(trimws(as.character(x)[1]))
    x <- gsub("[[:space:]-]+", "_", x)
    gsub("[^a-z0-9_]", "", x)
  }

  extract_family_raw <- function(family_input) {
    if (is_scalar_string(family_input)) return(family_input)

    if (is.list(family_input) && !is.null(family_input$family) &&
        length(family_input$family) >= 1L) {
      return(as.character(family_input$family[[1]]))
    }

    stopf(
      "Unsupported `family` specification. Supply one supported canonical ",
      "family name or a family object containing a `family` field."
    )
  }

  resolve_family_info <- function(family_input) {
    raw <- extract_family_raw(family_input)
    z <- normalize_token(raw)

    if (z %in% c("ordinal", "cumulative")) {
      return(list(name = "ordinal", variant = "ordinal", raw = raw))
    }

    if (z %in% c(
      "nb1", "nbinom1", "negative_binomial_1", "negativebinomial1",
      "negative_binomial_type_1", "negativebinomialtype1"
    )) {
      return(list(name = "negative_binomial", variant = "NB1", raw = raw))
    }

    if (z %in% c(
      "nb2", "negbin", "nbinom", "nbinom2", "negative_binomial",
      "negative_binomial_2", "negativebinomial", "negativebinomial2",
      "negative_binomial_type_2", "negativebinomialtype2"
    )) {
      return(list(name = "negative_binomial", variant = "NB2", raw = raw))
    }

    if (grepl("negative_?binomial", z)) {
      stopf(
        "Unknown negative-binomial parameterization: '", raw, "'. ",
        "EpiExposure v1 accepts only the explicitly recognized NB2 variants."
      )
    }

    if (z %in% c("beta", "beta_family", "beta_proportion", "beta_regression", "betar", "beta_resp")) {
      return(list(name = "beta", variant = "mean_precision", raw = raw))
    }

    if (z %in% c("binomial", "bernoulli")) {
      return(list(name = "binomial", variant = "bernoulli", raw = raw))
    }

    if (z == "poisson") {
      return(list(name = "poisson", variant = "poisson", raw = raw))
    }

    if (z == "gamma") {
      return(list(name = "gamma", variant = "gamma", raw = raw))
    }

    if (z %in% c("gaussian", "normal")) {
      return(list(name = "gaussian", variant = "gaussian", raw = raw))
    }

    stopf(
      "Unsupported family: '", raw, "'. EpiExposure v1 supports: beta, ",
      "binomial, poisson, gamma, gaussian, and negative_binomial (NB2)."
    )
  }

  extract_input_link <- function(family_input) {
    if (is.list(family_input) && !is.null(family_input$link) &&
        length(family_input$link) >= 1L) {
      out <- tolower(trimws(as.character(family_input$link[[1]])))
      if (!is.na(out) && nzchar(out)) return(out)
    }
    NULL
  }

  default_family_link <- function(family_name) {
    switch(
      family_name,
      beta = "logit",
      binomial = "logit",
      poisson = "log",
      gamma = "log",
      gaussian = "identity",
      negative_binomial = "log",
      stopf("Could not determine the default link for family '", family_name, "'.")
    )
  }

  validate_link_name <- function(link_name) {
    if (!is_scalar_string(link_name)) stopf("Could not determine a valid model link.")
    link_name <- tolower(link_name)
    link_name
  }

  validate_family_link <- function(
    family_name,
    link_name,
    model_engine
  ) {
    engine_links <- list(
      glm = list(
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity")
      ),
      glmmTMB = list(
        beta = c("logit", "probit", "cloglog", "identity", "inverse", "sqrt"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      gam = list(
        beta = c("logit", "probit", "cloglog", "cauchit"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      gamm = list(
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity")
      ),
      gls = list(
        gaussian = "identity"
      ),
      spamm = list(
        beta = c("logit", "probit", "cloglog", "cauchit"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      brms = list(
        beta = c("logit", "probit", "cloglog", "cauchit"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      inla = list(
        beta = "logit",
        gaussian = "identity",
        binomial = "logit",
        poisson = "log",
        gamma = "log",
        negative_binomial = "log"
      ),
      bdlnm = list(
        beta = "logit",
        gaussian = "identity",
        binomial = "logit",
        poisson = "log",
        gamma = "log",
        negative_binomial = "log"
      )
    )

    registered_families <- engine_links[[model_engine]]

    if (is.null(registered_families)) {
      stopf("Unsupported `model_engine`: '", model_engine, "'.")
    }

    allowed <- registered_families[[family_name]]

    if (is.null(allowed)) {
      stopf(
        "Family '", family_name, "' is not supported for ",
        "model_engine = '", model_engine, "'. ",
        "Supported families: ",
        paste(names(registered_families), collapse = ", "), "."
      )
    }

    # GLS is intentionally restricted by the EpiExposure contract.
    if (!link_name %in% allowed) {
      stopf(
        "Link '", link_name,
        "' is not supported for family '", family_name,
        "' with `model_engine = '", model_engine, "'. ",
        "Supported link(s): ",
        paste(allowed, collapse = ", "),
        "."
      )
    }

    invisible(TRUE)
  }

  valid_scalar_name <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }

  valid_flag <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  is_whole_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) &&
      is.finite(x) && x == as.integer(x)
  }

  if (identical(validation_method, "k-fold")) {
    if (!is_whole_scalar(k) || k < 2L) {
      stop(
        "`k` must be an integer greater than or equal to 2 when ",
        "`validation_method = 'k-fold'`.",
        call. = FALSE
      )
    }
    k <- as.integer(k)
  }

  if (!is.numeric(test_prop) || length(test_prop) != 1L ||
      is.na(test_prop) || !is.finite(test_prop) ||
      test_prop <= 0 || test_prop >= 1) {
    stop(
      "`test_prop` must be one finite numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }
  test_prop <- as.numeric(test_prop)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
        seed != floor(seed)) {
      stop(
        "`seed` must be NULL or one non-negative integer representable by R.",
        call. = FALSE
      )
    }
    seed <- as.integer(seed)
  }

  .with_local_seed <- function(seed_value, code) {
    if (is.null(seed_value)) {
      return(force(code))
    }

    seed_exists <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )

    if (seed_exists) {
      old_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }

    on.exit(
      {
        if (seed_exists) {
          assign(
            ".Random.seed",
            old_seed,
            envir = .GlobalEnv
          )
        } else if (exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )) {
          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )

    set.seed(seed_value)
    force(code)
  }

  # Basic columns and user inputs

  if (!is.data.frame(data) || !nrow(data)) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }

  if (!valid_scalar_name(response)) {
    stop("`response` must be one non-empty column name.", call. = FALSE)
  }
  if (!valid_scalar_name(group)) {
    stop("`group` must be one non-empty column name.", call. = FALSE)
  }
  if (!valid_scalar_name(time)) {
    stop("`time` must be one non-empty column name.", call. = FALSE)
  }

  if (anyDuplicated(c(response, group, time))) {
    stop(
      "`response`, `group`, and `time` must name three distinct columns.",
      call. = FALSE
    )
  }

  missing_basic <- setdiff(c(response, group, time), names(data))
  if (length(missing_basic)) {
    stop(
      "Required columns missing from `data`: ",
      paste(missing_basic, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.character(vars) || !length(vars) || anyNA(vars) ||
      any(!nzchar(vars)) || anyDuplicated(vars)) {
    stop(
      "`vars` must contain unique non-empty variable names.",
      call. = FALSE
    )
  }

  role_overlap <- intersect(vars, c(response, group, time))
  if (length(role_overlap)) {
    stop(
      "Candidate exposure variables cannot also be `response`, `group`, or ",
      "`time`: ", paste(role_overlap, collapse = ", "), ".",
      call. = FALSE
    )
  }

  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop(
      "Variables missing from `data`: ",
      paste(missing_vars, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (anyNA(data[[group]])) {
    stop("The grouping column cannot contain missing values.", call. = FALSE)
  }

  if (!is.numeric(data[[time]]) || anyNA(data[[time]]) ||
      any(!is.finite(data[[time]]))) {
    stop(
      "The chronological time column must contain only finite numeric values.",
      call. = FALSE
    )
  }

  if (!is.numeric(data[[response]]) || anyNA(data[[response]]) ||
      any(!is.finite(data[[response]]))) {
    stop(
      "The response column must contain only finite numeric values.",
      call. = FALSE
    )
  }

  for (variable in vars) {
    if (!is.numeric(data[[variable]]) || anyNA(data[[variable]]) ||
        any(!is.finite(data[[variable]]))) {
      stop(
        "Exposure variable '", variable,
        "' must contain only finite numeric values.",
        call. = FALSE
      )
    }
  }

  # Canonical names are required internally by define_exposures/build design.
  # Exposures are deliberately not allowed to occupy these reserved roles.
  if (any(vars %in% c("epi_id", "time", "y"))) {
    stop(
      "Candidate exposure names cannot use the internal reserved names ",
      "'epi_id', 'time', or 'y'.",
      call. = FALSE
    )
  }

  # Coordinate names are validated after the main data roles are known, and
  # before any fold or candidate is constructed.
  spatial_spec <- .epix_validate_spatial_spec(
    data = data,
    model_engine = model_engine,
    spatial_effect = spatial_effect,
    spatial_structure = spatial_structure,
    spatial_group = spatial_group,
    forbidden = c(response, time, vars, "y_model", "y",
                  if (!identical(group, "epi_id")) "epi_id")
  )
  spatial_effect <- spatial_spec$effect
  spatial_structure <- spatial_spec$structure
  spatial_group <- spatial_spec$group
  spatial_term <- spatial_spec$term

  if (!is.numeric(max_lag) || !length(max_lag) || anyNA(max_lag) ||
      any(!is.finite(max_lag)) || any(max_lag < 0) ||
      any(max_lag != as.integer(max_lag))) {
    stop(
      "`max_lag` must contain non-negative finite integer lag values.",
      call. = FALSE
    )
  }
  max_lag <- as.integer(max(max_lag))
  required_history_length <- max_lag + 1L
  history_contract <- "all_fitted_exposures_same_exact_max_lag_plus_one"

  validate_df_grid <- function(x, argument) {
    if (!is.numeric(x) || !length(x) || anyNA(x) ||
        any(!is.finite(x)) || any(x <= 0) ||
        any(x != as.integer(x))) {
      stop(
        "`", argument, "` must contain positive finite integers.",
        call. = FALSE
      )
    }
    sort(unique(as.integer(x)))
  }

  df_var_grid <- validate_df_grid(df_var_grid, "df_var_grid")
  df_lag_grid <- validate_df_grid(df_lag_grid, "df_lag_grid")

  if (!is_whole_scalar(min_vars) || min_vars < 1) {
    stop("`min_vars` must be a positive integer.", call. = FALSE)
  }
  min_vars <- as.integer(min_vars)

  if (is.null(max_vars)) {
    max_vars <- length(vars)
  }

  if (!is_whole_scalar(max_vars) ||
      max_vars < min_vars || max_vars > length(vars)) {
    stop(
      "`max_vars` must be an integer between `min_vars` and `length(vars)`.",
      call. = FALSE
    )
  }
  max_vars <- as.integer(max_vars)

  if (!is.null(random_effect)) {
    if (!valid_scalar_name(random_effect)) {
      stop(
        "`random_effect` must be NULL or one non-empty column name.",
        call. = FALSE
      )
    }

    if (!random_effect %in% names(data)) {
      stop(
        "`random_effect` ('", random_effect, "') was not found in `data`.",
        call. = FALSE
      )
    }

    if (random_effect %in% c(response, time) ||
        random_effect %in% vars) {
      stop(
        "`random_effect` cannot also be the response, time column, or a ",
        "candidate exposure variable.",
        call. = FALSE
      )
    }

    if (startsWith(random_effect, "cb_")) {
      stop(
        "`random_effect` cannot start with the reserved cross-basis prefix ",
        "'cb_'.",
        call. = FALSE
      )
    }

    if (anyNA(data[[random_effect]])) {
      stop("`random_effect` cannot contain missing values.", call. = FALSE)
    }

    if (random_effect == "epi_id" && group != "epi_id") {
      stop(
        "`random_effect = 'epi_id'` is ambiguous when `group` is a different ",
        "column because `epi_id` is an internal canonical alias. Use the actual ",
        "grouping-column name intended for the random intercept.",
        call. = FALSE
      )
    }
  }

  if (!is.null(random_effect) && model_engine %in% c("glm", "gls")) {
    stop(
      "`model_engine = '", model_engine,
      "'` does not support `random_effect` in EpiExposure v1.",
      call. = FALSE
    )
  }

  if (is.null(random_effect) && identical(model_engine, "gamm")) {
    stop(
      "`model_engine = 'gamm'` requires one random-intercept column in the ",
      "current EpiExposure v1 interface.",
      call. = FALSE
    )
  }

  # Random-effect prior contract

  if (!is.null(random_effect_prior)) {

    if (is.null(random_effect)) {
      stop(
        "`random_effect_prior` requires a non-NULL `random_effect`.",
        call. = FALSE
      )
    }

    if (!model_engine %in% c("inla", "bdlnm")) {
      stop(
        "`random_effect_prior` is currently supported only for ",
        "`model_engine = 'inla'` or `model_engine = 'bdlnm'`.",
        call. = FALSE
      )
    }

    if (!is.list(random_effect_prior) ||
        !length(random_effect_prior) ||
        is.null(names(random_effect_prior)) ||
        anyNA(names(random_effect_prior)) ||
        any(!nzchar(names(random_effect_prior))) ||
        anyDuplicated(names(random_effect_prior))) {
      stop(
        "`random_effect_prior` must be NULL or a non-empty named list ",
        "accepted by `INLA::f(..., hyper = ...)`.",
        call. = FALSE
      )
    }
  }

  if (!is_whole_scalar(min_success) || min_success < 2) {
    stop(
      "`min_success` must be an integer greater than or equal to 2.",
      call. = FALSE
    )
  }
  min_success <- as.integer(min_success)

  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) ||
      top_n <= 0 ||
      (is.finite(top_n) && top_n != as.integer(top_n))) {
    stop(
      "`top_n` must be a positive integer or `Inf`.",
      call. = FALSE
    )
  }
  if (is.finite(top_n)) top_n <- as.integer(top_n)

  if (!valid_flag(keep_fits)) {
    stop("`keep_fits` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!valid_flag(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  # All arguments passed to fit_epidlnm() through ... must be named and must not
  # duplicate the arguments managed by find_bestfit().
  if (length(fit_dots)) {
    dot_names <- names(
      fit_dots
    )

    if (is.null(dot_names) ||
        anyNA(dot_names) ||
        any(!nzchar(dot_names))) {
      stop(
        "All arguments supplied through `...` must be explicitly named.",
        call. = FALSE
      )
    }

    if (anyDuplicated(dot_names)) {
      duplicated_dot_names <- unique(
        dot_names[
          duplicated(dot_names)
        ]
      )

      stop(
        "Arguments supplied through `...` must have unique names. ",
        "Duplicated argument(s): ",
        paste(
          duplicated_dot_names,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    reserved_fit_args <- c(
      "data",
      "model_engine",
      "family",
      "random_effect",
      "random_effect_prior",
      "spatial_effect",
      "spatial_structure",
      "spatial_group",
      "epiexposure_spec",
      "basis_objects"
    )

    conflicting_fit_args <- intersect(
      dot_names,
      reserved_fit_args
    )

    if (length(conflicting_fit_args)) {
      stop(
        "Argument(s) managed by `find_bestfit()` cannot be supplied again ",
        "through `...`: ",
        paste(
          conflicting_fit_args,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }
  }

  # Canonical family and early response/engine validation

  family_info <- resolve_family_info(family)
  family_name <- family_info$name

  if (identical(family_name, "ordinal")) {
    stopf(
      "Ordinal outcomes are not supported in EpiExposure v1. This deliberate ",
      "restriction keeps fitting, prediction, uncertainty, performance metrics, ",
      "and ensemble behavior harmonized across supported families."
    )
  }

  if (identical(family_info$variant, "NB1")) {
    stopf(
      "Negative-binomial NB1 was requested, but EpiExposure v1 standardizes ",
      "`negative_binomial` to the quadratic-variance NB2 parameterization. ",
      "Use an NB2-compatible family specification."
    )
  }

  input_link <- extract_input_link(family)

  link_name <- validate_link_name(
    input_link %||% default_family_link(family_name)
  )

  # Reject unsupported engine/family/link combinations BEFORE starting any
  # validation candidate. fit_epidlnm() repeats this validation before native fitting.
  validate_family_link(
    family_name = family_name,
    link_name = link_name,
    model_engine = model_engine
  )

  # The INLA availability check must cover both the direct and bdlnm routes.
  # Fit-time validation remains active as the final backend-specific barrier.
  if (model_engine %in% c("inla", "bdlnm")) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stopf(
        "Package 'INLA' is required to check likelihood availability for ",
        "model_engine = '", model_engine, "'."
      )
    }

    engine_family <- switch(
      family_name,
      beta = "beta",
      gaussian = "gaussian",
      poisson = "poisson",
      gamma = "gamma",
      binomial = "binomial",
      negative_binomial = "nbinomial"
    )

    available_likelihoods <- names(INLA::inla.models()$likelihood)

    if (!engine_family %in% available_likelihoods) {
      stopf(
        "Likelihood '", engine_family, "' for EpiExposure family '",
        family_name, "' and model_engine = '", model_engine, "' is unavailable. ",
        "The current INLA installation does not provide this likelihood."
      )
    }

    # Validate conflicting INLA controls once, before creating validation partitions. Do not
    # change fit_dots: fit_epidlnm() injects the required settings for every fit.
    control_family <- fit_dots$control.family %||% list()
    if (!is.list(control_family)) {
      stopf("INLA `control.family` supplied through `...` must be a list.")
    }

    control_link <- control_family$control.link %||% list()
    if (!is.list(control_link)) {
      stopf("INLA `control.family$control.link` must be a list.")
    }

    if (!is.null(control_link$model)) {
      supplied_link <- tolower(as.character(control_link$model)[1])
      if (!identical(supplied_link, link_name)) {
        stopf(
          "Conflicting INLA link specifications: `family` implies '", link_name,
          "' but `control.family$control.link$model` is '", supplied_link, "'."
        )
      }
    }

    control_compute <- fit_dots$control.compute %||% list()
    if (!is.list(control_compute)) {
      stopf(
        "INLA-backed engines require `control.compute` supplied through `...` ",
        "to be a list."
      )
    }

    if (!is.null(control_compute$config) && !isTRUE(control_compute$config)) {
      stopf(
        "EpiExposure requires `control.compute$config = TRUE` for INLA-backed ",
        "engines so posterior coefficient draws can be generated for ",
        "downstream uncertainty."
      )
    }
  }

  outcome_type <- .resolve_outcome_type(family_name)

  # Canonical working data

  data_long <- data
  data_long$epi_id <- data_long[[group]]
  data_long$time <- data_long[[time]]
  data_long$y <- data_long[[response]]
  data_long <- data_long[
    order(data_long$epi_id, data_long$time),
    ,
    drop = FALSE
  ]
  rownames(data_long) <- NULL

  response_check <- unique(
    data_long[, c("epi_id", "y"), drop = FALSE]
  )

  if (anyDuplicated(response_check$epi_id)) {
    stop(
      "Multiple distinct response values were found within at least one ",
      "validation group. One group-level outcome is required.",
      call. = FALSE
    )
  }

  y_group <- response_check$y

  if (identical(family_name, "beta") &&
      any(y_group <= 0 | y_group >= 1)) {
    stop(
      "`family = 'beta'` requires every group-level response to lie strictly ",
      "inside (0, 1).",
      call. = FALSE
    )
  }

  if (identical(family_name, "binomial") &&
      !all(y_group %in% c(0, 1))) {
    stop(
      "`family = 'binomial'` requires the group-level response to be coded ",
      "0/1 so probability and classification metrics are comparable.",
      call. = FALSE
    )
  }

  if (family_name %in% c("poisson", "negative_binomial") &&
      (any(y_group < 0) || any(y_group != floor(y_group)))) {
    stop(
      "`family = '", family_name,
      "'` requires non-negative integer group-level responses.",
      call. = FALSE
    )
  }

  if (identical(family_name, "gamma") && any(y_group <= 0)) {
    stop(
      "`family = 'gamma'` requires strictly positive group-level responses.",
      call. = FALSE
    )
  }

  # Temporal completeness and common lag unit

  .epix_validate_regular_time(
    data = data_long,
    group_index = as.character(data_long$epi_id),
    time_col = "time"
  )

  group_counts <- table(as.character(data_long$epi_id))
  invalid_history_groups <- names(group_counts)[
    group_counts != required_history_length
  ]

  if (length(invalid_history_groups)) {
    examples <- paste0(
      invalid_history_groups,
      "=",
      as.integer(group_counts[invalid_history_groups])
    )

    stop(
      "Every validation group must contain exactly ",
      required_history_length,
      " equally spaced observations for max_lag = ",
      max_lag,
      ". Non-matching group(s) include: ",
      paste(utils::head(examples, 5L), collapse = ", "),
      if (length(examples) > 5L) "; ..." else ".",
      " Histories are not truncated, padded, or silently realigned.",
      call. = FALSE
    )
  }

  # Record the common temporal spacing for downstream diagnostics.
  time_step <- NA_real_
  for (id in unique(as.character(data_long$epi_id))) {
    tt <- sort(data_long$time[as.character(data_long$epi_id) == id])
    if (length(tt) > 1L) {
      time_step <- as.numeric(diff(tt)[1L])
      break
    }
  }

  # A random-intercept column used in one epidemic-level model row must itself
  # be unique/constant inside each validation group.
  if (!is.null(random_effect) && random_effect != group) {
    re_check <- unique(
      data_long[, c("epi_id", random_effect), drop = FALSE]
    )
    if (anyDuplicated(re_check$epi_id)) {
      stop(
        "`random_effect` must be constant within each validation group because ",
        "the fitted design contains one row per group.",
        call. = FALSE
      )
    }
  }

  # Check the entire original data once, before creating validation partitions/candidates.
  # Keep epidemic histories separate even when coordinate pairs are repeated.
  .epix_validate_spatial_constancy(
    data = data_long,
    spatial_effect = spatial_effect,
    spatial_group = spatial_group
  )

  # Ranking metric contract

  available_metrics <- .available_metrics(family_name)

  if (is.null(rank_metric)) {
    rank_metric <- if (identical(outcome_type, "binary")) {
      "ROC_AUC"
    } else {
      "CCC"
    }
  } else if (!is.character(rank_metric) || length(rank_metric) != 1L ||
             is.na(rank_metric) || !nzchar(rank_metric)) {
    stop(
      "`rank_metric` must be NULL or one non-empty metric name.",
      call. = FALSE
    )
  }

  if (!rank_metric %in% available_metrics) {
    stop(
      "`rank_metric = ", sQuote(rank_metric),
      "` is not available for family '", family_name,
      "'. Available metrics are: ",
      paste(available_metrics, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop(
      "`threshold` must be one finite numeric probability strictly between ",
      "0 and 1.",
      call. = FALSE
    )
  }

  rank_direction <- .metric_direction(rank_metric)

  if (verbose) {
    message(
      "Performance ranking: family = ", family_name,
      ", outcome_type = ", outcome_type,
      ", rank_metric = ", rank_metric,
      " (", rank_direction, ")."
    )

    if (identical(outcome_type, "binary")) {
      message(
        "Classification threshold = ", format(threshold),
        ". It is applied only to out-of-sample predicted probabilities when ",
        "classification metrics are computed; ROC_AUC, Brier, and LogLoss ",
        "remain probability-scale metrics."
      )
    }
  }

  # Strict training-template transport

  .validate_training_templates <- function(templates, variables) {

    if (!requireNamespace("dlnm", quietly = TRUE)) {
      stop(
        "Package 'dlnm' is required for DLNM cross-basis construction.",
        call. = FALSE
      )
    }

    if (!is.list(templates) || !length(templates) ||
        is.null(names(templates)) || anyNA(names(templates)) ||
        any(!nzchar(names(templates))) || anyDuplicated(names(templates))) {
      stop(
        "`define_exposures()` must return a non-empty named list of unique ",
        "cross-basis templates.",
        call. = FALSE
      )
    }

    if (!setequal(names(templates), variables)) {
      missing_templates <- setdiff(variables, names(templates))
      extra_templates <- setdiff(names(templates), variables)
      parts <- c(
        if (length(missing_templates)) {
          paste0("missing: ", paste(missing_templates, collapse = ", "))
        },
        if (length(extra_templates)) {
          paste0("unexpected: ", paste(extra_templates, collapse = ", "))
        }
      )
      stop(
        "Training basis-template names must match the candidate variables ",
        "exactly (", paste(parts, collapse = "; "), ").",
        call. = FALSE
      )
    }

    # Read list-level specification metadata before reordering because ordinary
    # list subsetting may drop custom attributes.
    stored_spec <- attr(templates, "spec", exact = TRUE)
    templates <- templates[variables]

    if (!is.list(stored_spec) || is.null(names(stored_spec)) ||
        anyNA(names(stored_spec)) || any(!nzchar(names(stored_spec))) ||
        anyDuplicated(names(stored_spec)) ||
        !setequal(names(stored_spec), variables)) {
      stop(
        "`define_exposures()` returned missing or inconsistent `spec` metadata.",
        call. = FALSE
      )
    }

    stored_spec <- stored_spec[variables]
    effective_spec <- stored_spec
    cb_columns <- character(0)

    template_details <- vector("list", length(variables))
    names(template_details) <- variables

    for (variable in variables) {
      template <- templates[[variable]]

      if (!inherits(template, "crossbasis")) {
        stop(
          "Training template for variable '", variable,
          "' must inherit from `crossbasis`.",
          call. = FALSE
        )
      }

      template_argvar <- attr(template, "argvar", exact = TRUE)
      template_arglag <- attr(template, "arglag", exact = TRUE)
      template_lag <- attr(template, "lag", exact = TRUE)

      if (!is.list(template_argvar) || !is.list(template_arglag)) {
        stop(
          "Training cross-basis for variable '", variable,
          "' is missing effective `argvar`/`arglag` attributes.",
          call. = FALSE
        )
      }

      if (is.null(template_lag) || !is.numeric(template_lag) ||
          !length(template_lag) || anyNA(template_lag) ||
          any(!is.finite(template_lag))) {
        stop(
          "Training cross-basis for variable '", variable,
          "' is missing a valid lag attribute.",
          call. = FALSE
        )
      }

      template_max_lag <- as.integer(max(template_lag))

      if (!identical(template_max_lag, max_lag)) {
        stop(
          "`max_lag` contradicts the training cross-basis for variable '",
          variable, "': requested ", max_lag,
          ", template maximum lag ", template_max_lag, ".",
          call. = FALSE
        )
      }

      current_spec <- stored_spec[[variable]]

      if (!is.list(current_spec) || is.null(current_spec$max_lag) ||
          is.null(current_spec$argvar) || is.null(current_spec$arglag)) {
        stop(
          "Incomplete training exposure specification for variable '",
          variable, "'.",
          call. = FALSE
        )
      }

      if (!is.numeric(current_spec$max_lag) ||
          length(current_spec$max_lag) != 1L ||
          is.na(current_spec$max_lag) ||
          !is.finite(current_spec$max_lag) ||
          current_spec$max_lag < 0 ||
          current_spec$max_lag != as.integer(current_spec$max_lag)) {
        stop(
          "Training exposure specification for variable '", variable,
          "' must contain one non-negative integer `max_lag`.",
          call. = FALSE
        )
      }

      spec_max_lag <- as.integer(current_spec$max_lag)

      if (!identical(spec_max_lag, template_max_lag) ||
          !identical(spec_max_lag, max_lag)) {
        stop(
          "Exposure specification, cross-basis template, and requested ",
          "`max_lag` disagree for variable '", variable, "'.",
          call. = FALSE
        )
      }

      if (ncol(template) < 1L) {
        stop(
          "Training cross-basis for variable '", variable,
          "' contains no columns.",
          call. = FALSE
        )
      }

      native_names <- colnames(template)
      if (!is.null(native_names) &&
          (length(native_names) != ncol(template) ||
           anyNA(native_names) || any(!nzchar(native_names)) ||
           anyDuplicated(native_names))) {
        stop(
          "Training cross-basis for variable '", variable,
          "' has invalid native column names.",
          call. = FALSE
        )
      }

      canonical_names <- paste0(
        "cb_", variable, "_", seq_len(ncol(template))
      )

      if (anyDuplicated(canonical_names)) {
        stop(
          "Canonical cross-basis names are not unique for variable '",
          variable, "'.",
          call. = FALSE
        )
      }

      # The returned crossbasis attributes are authoritative because onebasis()
      # may normalize/modify arguments relative to the original call.
      effective_spec[[variable]]$max_lag <- template_max_lag
      effective_spec[[variable]]$argvar <- template_argvar
      effective_spec[[variable]]$arglag <- template_arglag

      template_details[[variable]] <- list(
        max_lag = template_max_lag,
        argvar = template_argvar,
        arglag = template_arglag,
        ncol = ncol(template),
        native_names = native_names,
        canonical_names = canonical_names
      )

      cb_columns <- c(cb_columns, canonical_names)
    }

    if (anyDuplicated(cb_columns)) {
      stop(
        "Canonical cross-basis column names are duplicated across candidate ",
        "variables.",
        call. = FALSE
      )
    }

    # Keep the effective specification attached to the ordered template list as
    # well. The crossbasis objects themselves are unchanged.
    attr(templates, "spec") <- effective_spec

    list(
      templates = templates,
      spec = effective_spec,
      details = template_details,
      cb_cols = cb_columns,
      variables = variables,
      max_lag = max_lag,
      history_length = required_history_length,
      history_contract = history_contract
    )
  }

  .build_design_from_templates <- function(
    input_data,
    template_info,
    include_response = TRUE
  ) {

    if (!is.data.frame(input_data) || !nrow(input_data)) {
      stop(
        "Internal design construction requires non-empty long-format data.",
        call. = FALSE
      )
    }

    required <- c("epi_id", "time", template_info$variables)
    if (include_response) required <- c(required, "y")

    missing <- setdiff(required, names(input_data))
    if (length(missing)) {
      stop(
        "Internal design data are missing required column(s): ",
        paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }

    .epix_validate_regular_time(
      data = input_data,
      group_index = as.character(input_data$epi_id),
      time_col = "time"
    )

    counts <- table(as.character(input_data$epi_id))
    invalid_counts <- names(counts)[
      counts != template_info$history_length
    ]

    if (length(invalid_counts)) {
      examples <- paste0(
        invalid_counts,
        "=",
        as.integer(counts[invalid_counts])
      )

      stop(
        "Every group used to construct a candidate design must contain exactly ",
        template_info$history_length,
        " observations for template max_lag = ",
        template_info$max_lag,
        ". Non-matching group(s) include: ",
        paste(utils::head(examples, 5L), collapse = ", "),
        if (length(examples) > 5L) "; ..." else ".",
        call. = FALSE
      )
    }

    group_ids <- unique(input_data$epi_id)
    n_groups <- length(group_ids)

    out <- data.frame(
      epi_id = group_ids,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    for (variable in template_info$variables) {
      detail <- template_info$details[[variable]]
      matrix_var <- matrix(
        NA_real_,
        nrow = n_groups,
        ncol = detail$ncol
      )
      colnames(matrix_var) <- detail$canonical_names

      for (i in seq_along(group_ids)) {
        idx <- which(input_data$epi_id == group_ids[i])
        idx <- idx[order(input_data$time[idx])]
        x <- input_data[[variable]][idx]

        if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
          stop(
            "Exposure variable '", variable,
            "' contains invalid values while constructing a candidate design.",
            call. = FALSE
          )
        }

        cb_new <- dlnm::crossbasis(
          x,
          lag = detail$max_lag,
          argvar = detail$argvar,
          arglag = detail$arglag
        )

        if (ncol(cb_new) != detail$ncol) {
          stop(
            "Reconstructed cross-basis dimension mismatch for variable '",
            variable, "': training template has ", detail$ncol,
            " column(s), reconstructed basis has ", ncol(cb_new), ".",
            call. = FALSE
          )
        }

        reconstructed_native_names <- colnames(cb_new)

        if (!is.null(detail$native_names)) {
          if (is.null(reconstructed_native_names) ||
              !identical(reconstructed_native_names, detail$native_names)) {
            stop(
              "Reconstructed cross-basis native column names/order do not ",
              "match the training template for variable '", variable, "'.",
              call. = FALSE
            )
          }
        }

        cb_lag <- attr(cb_new, "lag", exact = TRUE)
        if (is.null(cb_lag) ||
            as.integer(max(cb_lag)) != detail$max_lag) {
          stop(
            "Reconstructed cross-basis lag does not match the training ",
            "template for variable '", variable, "'.",
            call. = FALSE
          )
        }

        final_row <- as.numeric(
          cb_new[nrow(cb_new), , drop = TRUE]
        )

        if (length(final_row) != detail$ncol ||
            any(!is.finite(final_row))) {
          stop(
            "Could not obtain a finite final cross-basis row for variable '",
            variable, "'.",
            call. = FALSE
          )
        }

        matrix_var[i, ] <- final_row
      }

      for (column in detail$canonical_names) {
        out[[column]] <- matrix_var[, column]
      }
    }

    actual_cb <- grep("^cb_", names(out), value = TRUE)

    if (!identical(actual_cb, template_info$cb_cols)) {
      stop(
        "Internal cross-basis design columns are inconsistent with the ",
        "training-template column contract.",
        call. = FALSE
      )
    }

    if (include_response) {
      y_values <- numeric(n_groups)

      for (i in seq_along(group_ids)) {
        yy <- unique(input_data$y[input_data$epi_id == group_ids[i]])

        if (length(yy) != 1L || !is.finite(yy)) {
          stop(
            "Exactly one finite response value is required per epidemic.",
            call. = FALSE
          )
        }

        y_values[i] <- yy
      }

      out$y <- y_values
    }

    # Store the actual ORIGINAL training templates. Do not attach a crossbasis
    # reconstructed from the first group.
    attr(out, "cb_templates") <- template_info$templates
    attr(out, "epiexposure_spec") <- template_info$spec
    attr(out, "epiexposure_cb_cols") <- template_info$cb_cols
    attr(out, "epiexposure_max_lag") <- template_info$max_lag
    attr(out, "epiexposure_history_length") <- template_info$history_length
    attr(out, "epiexposure_history_contract") <- template_info$history_contract

    out
  }

  .build_bdlnm_basis_objects <- function(
    input_data,
    template_info,
    design
  ) {
    if (!identical(model_engine, "bdlnm")) {
      return(template_info$templates)
    }

    if (!requireNamespace("dlnm", quietly = TRUE)) {
      stop(
        "Package 'dlnm' is required to construct bdlnm basis objects.",
        call. = FALSE
      )
    }

    group_ids <- design$epi_id
    basis_objects <- vector(
      "list",
      length(template_info$variables)
    )
    names(basis_objects) <- template_info$variables

    for (variable in template_info$variables) {
      detail <- template_info$details[[variable]]
      L <- detail$max_lag

      # Matrix-form crossbasis input represents one complete exposure history
      # per epidemic. Columns are ordered by retrospective lag 0,...,L, whereas
      # the package's long input is chronological oldest,...,most recent.
      histories <- matrix(
        NA_real_,
        nrow = length(group_ids),
        ncol = L + 1L
      )

      for (i in seq_along(group_ids)) {
        idx <- which(input_data$epi_id == group_ids[i])
        idx <- idx[order(input_data$time[idx])]
        x <- input_data[[variable]][idx]

        if (length(x) != L + 1L || any(!is.finite(x))) {
          stop(
            "bdlnm exposure history for variable '", variable,
            "' must contain exactly ", L + 1L,
            " finite observations for max_lag = ", L, ".",
            call. = FALSE
          )
        }

        histories[i, ] <- rev(x)
      }

      cb_epidemic <- dlnm::crossbasis(
        histories,
        lag = c(0L, L),
        argvar = detail$argvar,
        arglag = detail$arglag
      )

      if (!inherits(cb_epidemic, "crossbasis") ||
          nrow(cb_epidemic) != nrow(design) ||
          ncol(cb_epidemic) != detail$ncol) {
        stop(
          "Epidemic-level bdlnm cross-basis dimension mismatch for variable '",
          variable, "'.",
          call. = FALSE
        )
      }

      if (!is.null(detail$native_names)) {
        if (is.null(colnames(cb_epidemic)) ||
            !identical(colnames(cb_epidemic), detail$native_names)) {
          stop(
            "Epidemic-level bdlnm cross-basis column names/order do not match ",
            "the training template for variable '", variable, "'.",
            call. = FALSE
          )
        }
      }

      canonical <- detail$canonical_names
      design_values <- as.matrix(
        design[, canonical, drop = FALSE]
      )
      basis_values <- as.matrix(cb_epidemic)

      difference <- abs(design_values - basis_values)
      tolerance <- 1e-8 * pmax(
        1,
        abs(design_values),
        abs(basis_values)
      )

      if (any(!is.finite(basis_values)) ||
          any(difference > tolerance)) {
        stop(
          "The epidemic-level bdlnm cross-basis is not numerically equivalent ",
          "to the candidate design for variable '", variable,
          "'. Check exposure-history orientation and fitted basis metadata.",
          call. = FALSE
        )
      }

      # Native dlnm cross-basis names, such as v1.l1, are repeated when
      # multiple exposure bases are included in the same bdlnm model.
      # INLA requires unique internal keys. Prefix the native column names
      # with the exposure name after numerical equivalence has been checked.

      native_basis_names <- colnames(
        cb_epidemic
      )

      if (is.null(native_basis_names)) {
        native_basis_names <- paste0(
          "basis_",
          seq_len(
            ncol(cb_epidemic)
          )
        )
      }

      unique_basis_names <- paste0(
        variable,
        "_",
        native_basis_names
      )

      if (anyNA(unique_basis_names) ||
          any(!nzchar(unique_basis_names)) ||
          anyDuplicated(unique_basis_names)) {
        stop(
          "Could not create unique bdlnm cross-basis column names for ",
          "variable '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      colnames(
        cb_epidemic
      ) <- unique_basis_names

      basis_objects[[variable]] <- cb_epidemic
    }

    # Validate global uniqueness across all bdlnm cross-basis objects.
    # Each individual basis has already received the exposure-specific
    # prefix above. This final check ensures that no internal column name
    # is duplicated across different exposure bases.

    all_bdlnm_basis_names <- unlist(
      lapply(
        basis_objects,
        colnames
      ),
      use.names = FALSE
    )

    if (!length(all_bdlnm_basis_names)) {
      stop(
        "No internal column names were found in the candidate bdlnm ",
        "cross-basis objects.",
        call. = FALSE
      )
    }

    if (anyNA(all_bdlnm_basis_names) ||
        any(!nzchar(all_bdlnm_basis_names)) ||
        anyDuplicated(all_bdlnm_basis_names)) {

      duplicated_names <- unique(
        all_bdlnm_basis_names[
          duplicated(all_bdlnm_basis_names)
        ]
      )

      if (length(duplicated_names)) {

        stop(
          "The candidate bdlnm cross-basis objects contain duplicated ",
          "internal column names: ",
          paste(
            duplicated_names,
            collapse = ", "
          ),
          ".",
          call. = FALSE
        )

      } else {

        stop(
          "The candidate bdlnm cross-basis objects contain missing or ",
          "empty internal column names.",
          call. = FALSE
        )
      }
    }

    attr(
      basis_objects,
      "spec"
    ) <- template_info$spec

    basis_objects
  }

  .attach_spatial_metadata <- function(design, source_data) {
    if (is.null(spatial_effect)) return(design)

    metadata_columns <- unique(c("epi_id", spatial_effect, spatial_group))
    missing <- setdiff(metadata_columns, names(source_data))
    if (length(missing)) {
      stop("Spatial metadata missing from candidate data: ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }

    # Validate each metadata column independently and report the precise id.
    for (column in setdiff(metadata_columns, "epi_id")) {
      unique_pairs <- unique(source_data[, c("epi_id", column), drop = FALSE])
      bad <- duplicated(as.character(unique_pairs$epi_id))
      if (any(bad)) {
        stop("Spatial metadata column '", column,
             "' is not constant within epi_id '",
             as.character(unique_pairs$epi_id[which(bad)[1L]]), "'.",
             call. = FALSE)
      }
    }

    metadata_rows <- unique(source_data[, metadata_columns, drop = FALSE])
    index <- match(as.character(design$epi_id),
                   as.character(metadata_rows$epi_id))
    if (anyNA(index) || anyDuplicated(as.character(metadata_rows$epi_id))) {
      stop("Could not uniquely align spatial metadata with epidemic-level design rows.",
           call. = FALSE)
    }
    for (column in setdiff(metadata_columns, "epi_id")) {
      design[[column]] <- metadata_rows[[column]][index]
    }
    design
  }

  .attach_random_effect_metadata <- function(design, source_data) {
    if (is.null(random_effect)) {
      return(design)
    }

    # When the random intercept is the same as the validation group, the canonical
    # epidemic id already contains exactly the required grouping information.
    if (identical(random_effect, group)) {
      if (!identical(group, "epi_id")) {
        design[[random_effect]] <- design$epi_id
      }
      return(design)
    }

    if (identical(random_effect, "epi_id")) {
      # This case is only reachable when group == "epi_id" because ambiguity was
      # rejected during validation.
      return(design)
    }

    metadata_rows <- unique(
      source_data[, c("epi_id", random_effect), drop = FALSE]
    )

    if (anyDuplicated(metadata_rows$epi_id)) {
      stop(
        "`random_effect` must be constant within each epidemic.",
        call. = FALSE
      )
    }

    idx <- match(design$epi_id, metadata_rows$epi_id)

    if (anyNA(idx)) {
      stop(
        "Could not align random-effect metadata with epidemic-level design rows.",
        call. = FALSE
      )
    }

    design[[random_effect]] <- metadata_rows[[random_effect]][idx]
    design
  }

  build_candidate <- function(
    input_data,
    variables,
    exposure_df,
    lag_df
  ) {

    templates <- define_exposures(
      data = input_data,
      vars = variables,
      max_lag = max_lag,
      df_var = exposure_df,
      df_lag = lag_df,
      fun_var = fun_var,
      fun_lag = fun_lag
    )

    template_info <- .validate_training_templates(
      templates = templates,
      variables = variables
    )

    design <- .build_design_from_templates(
      input_data = input_data,
      template_info = template_info,
      include_response = TRUE
    )

    design <- .attach_random_effect_metadata(
      design = design,
      source_data = input_data
    )
    design <- .attach_spatial_metadata(
      design = design,
      source_data = input_data
    )

    fit_basis_objects <- .build_bdlnm_basis_objects(
      input_data = input_data,
      template_info = template_info,
      design = design
    )

    prepared <- prepare_response(
      data = design,
      y_var = "y",
      family = family_name
    )

    if (!is.data.frame(prepared) || nrow(prepared) != nrow(design) ||
        !"y_model" %in% names(prepared)) {
      stop(
        "`prepare_response()` returned an invalid candidate fitting data set.",
        call. = FALSE
      )
    }

    if (!is.null(spatial_effect) &&
        !all(c(spatial_effect, spatial_group) %in% names(prepared))) {
      stop("`prepare_response()` removed required spatial fitting columns.",
           call. = FALSE)
    }

    list(
      templates = template_info$templates,
      basis_objects = fit_basis_objects,
      template_info = template_info,
      spec = template_info$spec,
      design = design,
      prepared = prepared,
      variables = variables
    )
  }

  fit_candidate <- function(candidate) {
    fit_call <- c(
      list(
        data = candidate$prepared,
        model_engine = model_engine,
        family = family,
        random_effect = random_effect,
        random_effect_prior = random_effect_prior,
        spatial_effect = spatial_effect,
        spatial_structure = if (is.null(spatial_effect)) "matern" else spatial_structure,
        spatial_group = spatial_group,
        epiexposure_spec = candidate$spec,
        basis_objects = candidate$basis_objects
      ),
      fit_dots
    )

    fitted <- do.call(fit_epidlnm, fit_call)

    # Enforce the strict metadata/prediction contract immediately, rather than
    # discovering an outdated fit only when a held-out prediction is attempted.
    meta <- .get_epiexposure_metadata(fitted)

    fitted_random_effect_prior <- attr(
      fitted,
      "epiexposure_random_effect_prior",
      exact = TRUE
    )

    if (!identical(
      fitted_random_effect_prior,
      random_effect_prior
    )) {
      stop(
        "The fitted model did not preserve the requested ",
        "`random_effect_prior` metadata.",
        call. = FALSE
      )
    }

    expected_random_effect_model <- if (
      !is.null(random_effect) &&
      model_engine %in% c(
        "inla",
        "bdlnm"
      )
    ) {
      "iid"
    } else {
      NULL
    }

    fitted_random_effect_model <- attr(
      fitted,
      "epiexposure_random_effect_model",
      exact = TRUE
    )

    if (!identical(
      fitted_random_effect_model,
      expected_random_effect_model
    )) {
      stop(
        "The fitted model did not preserve the expected ",
        "`epiexposure_random_effect_model` metadata.",
        call. = FALSE
      )
    }

    if (!setequal(meta$vars, candidate$variables)) {
      stop(
        "Fitted EpiExposure metadata do not match the candidate exposure set.",
        call. = FALSE
      )
    }

    if (!identical(meta$max_lag, candidate$template_info$max_lag) ||
        !identical(
          meta$history_length,
          candidate$template_info$history_length
        ) ||
        !identical(
          meta$history_contract,
          candidate$template_info$history_contract
        )) {
      stop(
        "Fitted EpiExposure temporal metadata do not match the candidate ",
        "training-template exact-history contract.",
        call. = FALSE
      )
    }

    for (variable in candidate$variables) {
      if (!identical(
        meta$spec[[variable]]$max_lag,
        candidate$template_info$details[[variable]]$max_lag
      )) {
        stop(
          "Fitted metadata max_lag does not match the candidate training ",
          "template for variable '", variable, "'.",
          call. = FALSE
        )
      }
    }

    fitted
  }

  # Candidate variable sets

  if (is.null(var_sets)) {
    sizes <- seq.int(min_vars, max_vars)

    var_sets <- unlist(
      lapply(
        sizes,
        function(size) utils::combn(vars, size, simplify = FALSE)
      ),
      recursive = FALSE
    )
  } else {
    if (!is.list(var_sets) || !length(var_sets)) {
      stop(
        "`var_sets` must be NULL or a non-empty list.",
        call. = FALSE
      )
    }

    var_sets <- lapply(
      var_sets,
      function(set) {
        if (!is.character(set) || !length(set) || anyNA(set) ||
            any(!nzchar(set)) || anyDuplicated(set)) {
          stop(
            "Every element of `var_sets` must contain unique non-empty ",
            "variable names.",
            call. = FALSE
          )
        }

        missing <- setdiff(set, vars)
        if (length(missing)) {
          stop(
            "Variables in `var_sets` not listed in `vars`: ",
            paste(missing, collapse = ", "), ".",
            call. = FALSE
          )
        }

        # Keep the order supplied in `vars` so basis/column ordering is stable
        # across candidates, folds, and full-data refits.
        vars[vars %in% set]
      }
    )

    keys <- vapply(
      var_sets,
      paste,
      collapse = "||",
      character(1)
    )
    var_sets <- var_sets[!duplicated(keys)]
  }

  if (!length(var_sets)) {
    stop("No candidate variable sets were generated.", call. = FALSE)
  }

  # Grouped validation partition construction

  group_ids <- as.character(response_check$epi_id)
  group_response <- as.numeric(response_check$y)
  n_groups <- length(group_ids)

  if (n_groups < 2L) {
    stop(
      "At least two independent groups are required for grouped validation.",
      call. = FALSE
    )
  }

  if (anyDuplicated(group_ids)) {
    stop(
      "Internal validation error: group identifiers are not unique after ",
      "canonicalization.",
      call. = FALSE
    )
  }

  validation_stratified <- FALSE
  training_groups <- character(0)
  test_groups <- character(0)
  evaluation_groups <- character(0)
  n_training_groups <- NA_integer_
  n_test_groups <- NA_integer_
  effective_test_prop <- NA_real_

  if (identical(validation_method, "LOOCV")) {

    n_folds <- n_groups
    effective_k <- n_folds
    fold_assignment <- seq_len(n_groups)
    validation_scheme <- "leave_one_group_out"
    evaluation_groups <- group_ids

    validation_assignments <- data.frame(
      group = group_ids,
      partition = rep("evaluation", n_groups),
      fold = as.integer(fold_assignment),
      observed = group_response,
      stringsAsFactors = FALSE
    )

  } else if (identical(validation_method, "k-fold")) {

    if (k >= n_groups) {
      stop(
        "For `validation_method = 'k-fold'`, `k` must be smaller than the ",
        "number of independent groups. Received k = ", k, " and ", n_groups,
        " groups. Use `validation_method = 'LOOCV'` to leave one group out ",
        "at a time.",
        call. = FALSE
      )
    }

    n_folds <- k
    effective_k <- k
    evaluation_groups <- group_ids

    if (identical(outcome_type, "binary")) {

      class_counts <- table(
        factor(
          group_response,
          levels = c(0, 1)
        )
      )

      if (any(class_counts == 0L)) {
        stop(
          "Grouped stratified k-fold validation for a binomial outcome requires ",
          "both response classes 0 and 1 to be present.",
          call. = FALSE
        )
      }

      if (any(class_counts < k)) {
        minority_count <- min(as.integer(class_counts))
        stop(
          "Grouped stratified k-fold validation cannot place both binomial ",
          "classes in every fold because the less frequent class contains only ",
          minority_count, " group(s), while k = ", k, ". Choose `k <= ",
          minority_count, "` or use `validation_method = 'LOOCV'`.",
          call. = FALSE
        )
      }

      validation_stratified <- TRUE
      validation_scheme <- "stratified_grouped_k_fold"

      fold_assignment <- .with_local_seed(
        seed,
        {
          assignment <- integer(n_groups)
          total_fold_counts <- integer(k)

          # Process the larger stratum first. Within each class, groups are
          # randomized and greedily allocated to folds with the smallest current
          # class count, then smallest total size. This guarantees class counts
          # differing by at most one while also preserving overall fold balance.
          class_order <- c(0, 1)[
            order(
              as.integer(class_counts),
              decreasing = TRUE
            )
          ]

          for (class_value in class_order) {
            class_index <- which(group_response == class_value)
            class_index <- class_index[
              sample.int(length(class_index))
            ]

            class_fold_counts <- integer(k)

            for (idx in class_index) {
              candidates <- which(
                class_fold_counts == min(class_fold_counts)
              )

              minimum_total <- min(
                total_fold_counts[candidates]
              )

              candidates <- candidates[
                total_fold_counts[candidates] == minimum_total
              ]

              chosen <- candidates[
                sample.int(
                  length(candidates),
                  size = 1L
                )
              ]

              assignment[idx] <- chosen
              class_fold_counts[chosen] <-
                class_fold_counts[chosen] + 1L
              total_fold_counts[chosen] <-
                total_fold_counts[chosen] + 1L
            }
          }

          assignment
        }
      )

    } else {

      validation_scheme <- "grouped_k_fold"

      fold_assignment <- .with_local_seed(
        seed,
        {
          assignment <- integer(n_groups)
          shuffled_groups <- sample.int(n_groups)
          randomized_fold_order <- sample.int(k)
          labels <- rep(
            randomized_fold_order,
            length.out = n_groups
          )
          assignment[shuffled_groups] <- labels
          assignment
        }
      )
    }

    if (length(fold_assignment) != n_groups ||
        anyNA(fold_assignment) ||
        any(!fold_assignment %in% seq_len(n_folds))) {
      stop(
        "Internal validation error: invalid k-fold assignment.",
        call. = FALSE
      )
    }

    validation_assignments <- data.frame(
      group = group_ids,
      partition = rep("evaluation", n_groups),
      fold = as.integer(fold_assignment),
      observed = group_response,
      stringsAsFactors = FALSE
    )

  } else {

    # One grouped proportional holdout split is created once and reused by all
    # candidate models. Complete groups are never divided across partitions.
    n_test_target <- as.integer(round(n_groups * test_prop))
    n_test_target <- max(
      1L,
      min(n_groups - 1L, n_test_target)
    )

    if (n_test_target < 1L || n_test_target >= n_groups) {
      stop(
        "Grouped holdout validation could not retain at least one complete ",
        "group in both training and test partitions.",
        call. = FALSE
      )
    }

    n_folds <- 1L
    effective_k <- NA_integer_
    effective_test_prop <- test_prop

    if (identical(outcome_type, "binary")) {

      class_counts <- table(
        factor(
          group_response,
          levels = c(0, 1)
        )
      )

      if (any(class_counts < 2L)) {
        stop(
          "Grouped stratified holdout validation for a binomial outcome ",
          "requires at least two complete groups in each response class so ",
          "training and test partitions can both contain outcomes 0 and 1.",
          call. = FALSE
        )
      }

      validation_stratified <- TRUE
      validation_scheme <- "stratified_grouped_holdout"

      # At least two test groups and two training groups are necessary to place
      # one group from each response class in both partitions. Choose the closest
      # feasible test size to the user-requested proportional target.
      n_test_groups <- max(
        2L,
        min(n_groups - 2L, n_test_target)
      )

      n0 <- as.integer(class_counts[["0"]])
      n1 <- as.integer(class_counts[["1"]])

      feasible_t0 <- seq.int(1L, n0 - 1L)
      feasible_t1 <- seq.int(1L, n1 - 1L)
      allocation_grid <- expand.grid(
        n_test_0 = feasible_t0,
        n_test_1 = feasible_t1,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      allocation_grid$n_test_total <-
        allocation_grid$n_test_0 + allocation_grid$n_test_1
      allocation_grid$total_gap <- abs(
        allocation_grid$n_test_total - n_test_groups
      )
      allocation_grid$class_gap <-
        abs(allocation_grid$n_test_0 - n0 * test_prop) +
        abs(allocation_grid$n_test_1 - n1 * test_prop)
      allocation_grid$composition_gap <- abs(
        allocation_grid$n_test_0 / allocation_grid$n_test_total -
          n0 / n_groups
      )

      allocation_grid <- allocation_grid[
        order(
          allocation_grid$total_gap,
          allocation_grid$class_gap,
          allocation_grid$composition_gap,
          allocation_grid$n_test_0
        ),
        ,
        drop = FALSE
      ]

      selected_allocation <- allocation_grid[1L, , drop = FALSE]
      n_test_groups <- as.integer(selected_allocation$n_test_total[[1L]])
      n_test_0 <- as.integer(selected_allocation$n_test_0[[1L]])
      n_test_1 <- as.integer(selected_allocation$n_test_1[[1L]])
      n_train_groups <- n_groups - n_test_groups

      selected_test_indices <- .with_local_seed(
        seed,
        {
          index_0 <- which(group_response == 0)
          index_1 <- which(group_response == 1)
          c(
            index_0[sample.int(length(index_0), size = n_test_0)],
            index_1[sample.int(length(index_1), size = n_test_1)]
          )
        }
      )

    } else {

      validation_scheme <- "grouped_holdout"
      n_test_groups <- n_test_target
      n_train_groups <- n_groups - n_test_groups

      selected_test_indices <- .with_local_seed(
        seed,
        sample.int(
          n_groups,
          size = n_test_groups,
          replace = FALSE
        )
      )
    }

    if (n_test_groups < 1L || n_train_groups < 1L) {
      stop(
        "Grouped holdout validation could not retain at least one complete ",
        "group in both training and test partitions.",
        call. = FALSE
      )
    }

    selected_test_groups <- group_ids[selected_test_indices]
    test_groups <- group_ids[group_ids %in% selected_test_groups]
    training_groups <- group_ids[!group_ids %in% test_groups]
    evaluation_groups <- test_groups
    n_training_groups <- as.integer(length(training_groups))
    n_test_groups <- as.integer(length(test_groups))

    if (length(training_groups) + length(test_groups) != n_groups) {
      stop(
        "Internal holdout validation error: training and test group counts do ",
        "not sum to the total number of groups.",
        call. = FALSE
      )
    }

    if (length(intersect(training_groups, test_groups)) != 0L) {
      stop(
        "Internal holdout validation error: at least one group appears in both ",
        "training and test partitions.",
        call. = FALSE
      )
    }

    if (!setequal(union(training_groups, test_groups), group_ids)) {
      stop(
        "Internal holdout validation error: training and test partitions do not ",
        "cover exactly the original groups.",
        call. = FALSE
      )
    }

    if (length(test_groups) != n_test_groups) {
      stop(
        "Internal holdout validation error: the effective number of test groups ",
        "does not match the constructed test partition.",
        call. = FALSE
      )
    }

    if (identical(outcome_type, "binary")) {
      train_response <- group_response[match(training_groups, group_ids)]
      test_response <- group_response[match(test_groups, group_ids)]

      if (!setequal(unique(train_response), c(0, 1))) {
        stop(
          "Internal holdout validation error: the stratified training partition ",
          "does not contain both binomial outcomes 0 and 1.",
          call. = FALSE
        )
      }

      if (!setequal(unique(test_response), c(0, 1))) {
        stop(
          "Internal holdout validation error: the stratified test partition ",
          "does not contain both binomial outcomes 0 and 1.",
          call. = FALSE
        )
      }
    }

    validation_assignments <- data.frame(
      group = group_ids,
      partition = ifelse(
        group_ids %in% test_groups,
        "test",
        "train"
      ),
      fold = ifelse(
        group_ids %in% test_groups,
        1L,
        NA_integer_
      ),
      observed = group_response,
      stringsAsFactors = FALSE
    )
    validation_assignments$fold <- as.integer(validation_assignments$fold)
  }

  n_evaluation_groups <- as.integer(length(evaluation_groups))
  max_success_groups <- n_evaluation_groups

  if (min_success > max_success_groups) {
    stop(
      "`min_success` cannot exceed the number of groups assigned to ",
      "out-of-sample evaluation (", max_success_groups, ").",
      call. = FALSE
    )
  }

  if (identical(validation_method, "holdout")) {

    validation_balance <- data.frame(
      partition = c("train", "test"),
      fold = c(NA_integer_, 1L),
      n_groups = as.integer(c(
        length(training_groups),
        length(test_groups)
      )),
      stringsAsFactors = FALSE
    )

    if (identical(outcome_type, "binary")) {
      validation_balance$n_outcome_0 <- c(
        sum(validation_assignments$partition == "train" &
              validation_assignments$observed == 0),
        sum(validation_assignments$partition == "test" &
              validation_assignments$observed == 0)
      )
      validation_balance$n_outcome_1 <- c(
        sum(validation_assignments$partition == "train" &
              validation_assignments$observed == 1),
        sum(validation_assignments$partition == "test" &
              validation_assignments$observed == 1)
      )

      if (any(validation_balance$n_outcome_0 < 1L) ||
          any(validation_balance$n_outcome_1 < 1L)) {
        stop(
          "Internal holdout validation error: both binomial classes must be ",
          "present in both validation partitions.",
          call. = FALSE
        )
      }
    }

    if (sum(validation_balance$n_groups) != n_groups) {
      stop(
        "Internal holdout validation error: validation-balance counts do not ",
        "sum to the total number of groups.",
        call. = FALSE
      )
    }

    # Reuse the existing fold-oriented evaluation architecture. For holdout,
    # there is exactly one validation iteration containing all test groups.
    fold_groups <- list(test_groups)
    fold_rows <- list(
      which(as.character(data_long$epi_id) %in% test_groups)
    )

  } else {

    fold_sizes <- tabulate(
      validation_assignments$fold,
      nbins = n_folds
    )

    if (any(fold_sizes < 1L)) {
      stop(
        "Internal validation error: at least one evaluation fold contains no ",
        "groups.",
        call. = FALSE
      )
    }

    if (identical(validation_method, "k-fold") &&
        max(fold_sizes) - min(fold_sizes) > 1L) {
      stop(
        "Internal validation error: grouped k-fold allocation is not balanced ",
        "by number of groups.",
        call. = FALSE
      )
    }

    validation_balance <- data.frame(
      partition = rep("evaluation", n_folds),
      fold = seq_len(n_folds),
      n_groups = as.integer(fold_sizes),
      stringsAsFactors = FALSE
    )

    if (identical(outcome_type, "binary")) {
      validation_balance$n_outcome_0 <- vapply(
        seq_len(n_folds),
        function(fold_index) {
          sum(
            validation_assignments$fold == fold_index &
              validation_assignments$observed == 0
          )
        },
        integer(1)
      )

      validation_balance$n_outcome_1 <- vapply(
        seq_len(n_folds),
        function(fold_index) {
          sum(
            validation_assignments$fold == fold_index &
              validation_assignments$observed == 1
          )
        },
        integer(1)
      )

      if (identical(validation_method, "k-fold")) {
        if (max(validation_balance$n_outcome_0) -
            min(validation_balance$n_outcome_0) > 1L ||
            max(validation_balance$n_outcome_1) -
            min(validation_balance$n_outcome_1) > 1L ||
            any(validation_balance$n_outcome_0 < 1L) ||
            any(validation_balance$n_outcome_1 < 1L)) {
          stop(
            "Internal validation error: binomial k-fold stratification is not ",
            "balanced as required.",
            call. = FALSE
          )
        }
      }
    }

    row_group_match <- match(
      as.character(data_long$epi_id),
      validation_assignments$group
    )

    if (anyNA(row_group_match)) {
      stop(
        "Internal validation error: could not align long-format rows with ",
        "validation assignments.",
        call. = FALSE
      )
    }

    row_fold <- validation_assignments$fold[row_group_match]

    fold_rows <- lapply(
      seq_len(n_folds),
      function(fold_index) which(row_fold == fold_index)
    )

    fold_groups <- lapply(
      seq_len(n_folds),
      function(fold_index) {
        validation_assignments$group[
          !is.na(validation_assignments$fold) &
            validation_assignments$fold == fold_index
        ]
      }
    )
  }

  test_data_list <- lapply(
    fold_rows,
    function(idx) data_long[idx, , drop = FALSE]
  )

  if (verbose) {
    if (identical(validation_method, "LOOCV")) {
      message(
        "Validation: LOOCV with ", n_groups,
        " complete group(s); one complete group is held out per iteration."
      )
    } else if (identical(validation_method, "k-fold")) {
      message(
        "Validation: grouped ", k, "-fold with ", n_groups,
        " complete group(s); fold sizes = ",
        paste(validation_balance$n_groups, collapse = ", "),
        if (validation_stratified) {
          paste0(
            "; binomial class counts 0/1 by fold = ",
            paste(
              paste0(
                validation_balance$n_outcome_0,
                "/",
                validation_balance$n_outcome_1
              ),
              collapse = ", "
            )
          )
        } else {
          ""
        },
        "."
      )
    } else {
      message(
        "Validation: grouped proportional holdout with ", n_groups,
        " complete group(s); training = ", n_training_groups,
        ", test = ", n_test_groups,
        ", requested test_prop = ", format(test_prop),
        if (validation_stratified) {
          paste0(
            "; binomial class counts 0/1: train = ",
            validation_balance$n_outcome_0[
              validation_balance$partition == "train"
            ],
            "/",
            validation_balance$n_outcome_1[
              validation_balance$partition == "train"
            ],
            ", test = ",
            validation_balance$n_outcome_0[
              validation_balance$partition == "test"
            ],
            "/",
            validation_balance$n_outcome_1[
              validation_balance$partition == "test"
            ]
          )
        } else {
          ""
        },
        "."
      )
    }
  }

  # Warning classification

  warning_type_levels <- c(
    "knot",
    "convergence",
    "hessian",
    "other"
  )

  classify_warning <- function(message) {
    message_lower <- tolower(
      as.character(message)[1L]
    )

    # Check Hessian warnings before the broader convergence class because some
    # glmmTMB messages include both "convergence problem" and "Hessian".
    if (grepl("hessian", message_lower, fixed = TRUE) ||
        grepl("positive-definite", message_lower, fixed = TRUE) ||
        grepl("positive definite", message_lower, fixed = TRUE)) {
      return("hessian")
    }

    if (grepl(
      "shoving 'interior' knots matching boundary knots to inside",
      message_lower,
      fixed = TRUE
    ) ||
    (grepl("knot", message_lower, fixed = TRUE) &&
     grepl("boundary", message_lower, fixed = TRUE))) {
      return("knot")
    }

    if (grepl("converg", message_lower) ||
        grepl(
          "function evaluation limit reached",
          message_lower,
          fixed = TRUE
        ) ||
        grepl("iteration limit", message_lower, fixed = TRUE) ||
        grepl(
          "maximum number of iterations",
          message_lower,
          fixed = TRUE
        )) {
      return("convergence")
    }

    "other"
  }

  # Candidate grid

  # Preserve the original candidate order:
  # df_var (outer) -> df_lag -> variable set (inner).
  candidate_grid <- expand.grid(
    var_set_id = seq_along(var_sets),
    df_lag = df_lag_grid,
    df_var = df_var_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  candidate_grid <- candidate_grid[
    ,
    c("df_var", "df_lag", "var_set_id"),
    drop = FALSE
  ]

  candidate_grid$model_id <- seq_len(
    nrow(candidate_grid)
  )

  total_candidates <- nrow(candidate_grid)

  if (verbose) {
    message(
      "EpiExposure grouped validation search: engine = ", model_engine,
      ", family = ", family_name,
      ", link = ", link_name,
      ", candidates = ", total_candidates,
      ", ranking metric = ", rank_metric,
      " (", rank_direction, ")."
    )

    message("Candidate models:")
    for (candidate_index in seq_len(total_candidates)) {
      candidate_row <- candidate_grid[candidate_index, , drop = FALSE]
      variables <- var_sets[[candidate_row$var_set_id[[1L]]]]
      message(
        "  [", candidate_index, "/", total_candidates, "] df_var = ",
        candidate_row$df_var[[1L]],
        ", df_lag = ", candidate_row$df_lag[[1L]],
        ", vars = ", paste(variables, collapse = " + ")
      )
    }
  }

  # Candidate evaluator

  evaluate_candidate <- function(candidate_index) {

    candidate_start_time <- proc.time()[["elapsed"]]

    candidate_row <- candidate_grid[
      candidate_index,
      ,
      drop = FALSE
    ]

    model_id <- candidate_row$model_id[[1L]]
    exposure_df <- candidate_row$df_var[[1L]]
    lag_df <- candidate_row$df_lag[[1L]]
    variables <- var_sets[[candidate_row$var_set_id[[1L]]]]
    variables_label <- paste(
      variables,
      collapse = " + "
    )

    # Prediction bookkeeping remains indexed by every independent group.
    # Metrics are later restricted explicitly to `evaluation_groups`, so holdout
    # training groups are never counted as failed validation predictions.
    observed <- rep(
      NA_real_,
      n_groups
    )

    predicted <- rep(
      NA_real_,
      n_groups
    )

    success <- rep(
      FALSE,
      n_groups
    )

    fold_success <- rep(
      FALSE,
      n_folds
    )

    # Warning bookkeeping remains validation-iteration level. LOOCV and k-fold
    # use folds; holdout uses one training/test validation iteration.
    warning_fold <- rep(
      FALSE,
      n_folds
    )

    warning_events <- integer(
      n_folds
    )

    candidate_warning_messages <- vector(
      "list",
      n_folds
    )

    warning_type_fold <- matrix(
      FALSE,
      nrow = n_folds,
      ncol = length(warning_type_levels),
      dimnames = list(
        NULL,
        warning_type_levels
      )
    )

    warning_type_events <- matrix(
      0L,
      nrow = n_folds,
      ncol = length(warning_type_levels),
      dimnames = list(
        NULL,
        warning_type_levels
      )
    )

    # One list slot per validation iteration. Under k-fold and holdout a slot
    # may contain multiple group-level prediction/failure rows.
    candidate_predictions <- vector(
      "list",
      n_folds
    )

    candidate_failures <- vector(
      "list",
      n_folds + as.integer(keep_fits)
    )

    # One additional warning slot is reserved for performance metrics; another
    # is used for the full-data refit when keep_fits = TRUE.
    candidate_warnings <- vector(
      "list",
      n_folds + 1L + as.integer(keep_fits)
    )

    validation_stage <- validation_method

    # IMPORTANT: validation iterations remain sequential. One future owns one
    # candidate and executes its complete validation workflow.
    for (fold_index in seq_len(n_folds)) {

      test_groups <- fold_groups[[fold_index]]
      test_idx <- fold_rows[[fold_index]]

      # Negative integer subsetting preserves the existing architecture while
      # ensuring no held-out group contributes to training basis estimation.
      train_data <- data_long[
        -test_idx,
        ,
        drop = FALSE
      ]

      test_data <- test_data_list[[fold_index]]

      train_groups_iteration <- unique(as.character(train_data$epi_id))
      test_groups_iteration <- unique(as.character(test_data$epi_id))

      if (length(intersect(train_groups_iteration, test_groups_iteration)) != 0L) {
        stop(
          "Internal validation error: at least one group appears in both ",
          "training and evaluation data within the same iteration.",
          call. = FALSE
        )
      }

      if (!setequal(
        union(train_groups_iteration, test_groups_iteration),
        group_ids
      )) {
        stop(
          "Internal validation error: training and evaluation data do not cover ",
          "exactly all original groups within the validation iteration.",
          call. = FALSE
        )
      }

      if (!setequal(test_groups_iteration, test_groups)) {
        stop(
          "Internal validation error: evaluation rows do not match the groups ",
          "assigned to the current validation iteration.",
          call. = FALSE
        )
      }

      if (identical(validation_method, "holdout") &&
          !setequal(train_groups_iteration, training_groups)) {
        stop(
          "Internal holdout validation error: training rows do not match the ",
          "fixed holdout training partition.",
          call. = FALSE
        )
      }

      fold_warning_messages <- character(0)

      fold_result <- tryCatch(
        withCallingHandlers(
          {
            # Templates are estimated ONLY on training groups.
            training_candidate <- build_candidate(
              train_data,
              variables,
              exposure_df,
              lag_df
            )

            fitted_model <- fit_candidate(
              training_candidate
            )

            # Apply the exact effective training cross-basis definition to all
            # held-out groups in this fold. No df/knots/boundaries are estimated
            # from held-out exposure values.
            test_design <- .build_design_from_templates(
              input_data = test_data,
              template_info = training_candidate$template_info,
              include_response = TRUE
            )

            test_design <- .attach_random_effect_metadata(
              design = test_design,
              source_data = test_data
            )
            test_design <- .attach_spatial_metadata(
              design = test_design,
              source_data = test_data
            )

            prepared_test <- prepare_response(
              data = test_design,
              y_var = "y",
              family = family_name
            )

            if (!is.null(spatial_effect) &&
                !all(c(spatial_effect, spatial_group) %in% names(prepared_test))) {
              stop("`prepare_response()` removed spatial columns from test design.")
            }

            if (!is.data.frame(prepared_test) ||
                nrow(prepared_test) != length(test_groups) ||
                !"epi_id" %in% names(prepared_test) ||
                !"y_model" %in% names(prepared_test)) {
              stop(
                "Internal validation error: the held-out design did not ",
                "return exactly one prediction row per test group."
              )
            }

            prepared_groups <- as.character(
              prepared_test$epi_id
            )

            if (anyNA(prepared_groups) ||
                anyDuplicated(prepared_groups) ||
                !setequal(prepared_groups, test_groups)) {
              stop(
                "Internal validation error: held-out group identifiers ",
                "could not be aligned with the fold assignment."
              )
            }

            # Validation ranking uses only the central population/fixed component.
            # Coefficient/posterior uncertainty, conventional random effects,
            # and spatially autocorrelated effects are all excluded.
            prediction <- .predict_point_population(
              fit = fitted_model,
              newdata = prepared_test,
              type = "response",
              warn_fallback = TRUE
            )

            if (length(prediction) != nrow(prepared_test)) {
              stop(
                "Population-level held-out prediction returned ",
                length(prediction), " value(s) for ",
                nrow(prepared_test), " held-out group(s)."
              )
            }

            alignment <- match(
              test_groups,
              prepared_groups
            )

            if (anyNA(alignment)) {
              stop(
                "Internal validation error: failed to align predicted ",
                "rows with test groups."
              )
            }

            list(
              ok = TRUE,
              groups = test_groups,
              obs = as.numeric(
                prepared_test$y_model[alignment]
              ),
              pred = as.numeric(
                prediction[alignment]
              )
            )
          },
          warning = function(w) {
            fold_warning_messages <<- c(
              fold_warning_messages,
              conditionMessage(w)
            )
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          list(
            ok = FALSE,
            error = conditionMessage(e)
          )
        }
      )

      if (length(fold_warning_messages)) {

        warning_fold[fold_index] <- TRUE
        warning_events[fold_index] <-
          length(fold_warning_messages)

        candidate_warning_messages[[fold_index]] <-
          fold_warning_messages

        fold_warning_types <- vapply(
          fold_warning_messages,
          classify_warning,
          character(1)
        )

        fold_type_counts <- table(
          factor(
            fold_warning_types,
            levels = warning_type_levels
          )
        )

        warning_type_fold[fold_index, ] <-
          as.integer(fold_type_counts) > 0L

        warning_type_events[fold_index, ] <-
          as.integer(fold_type_counts)

        candidate_warnings[[fold_index]] <- data.frame(
          model_id = model_id,
          fold = fold_index,
          group = if (length(test_groups) == 1L) {
            as.character(test_groups)
          } else {
            NA_character_
          },
          stage = validation_stage,
          df_var = exposure_df,
          df_lag = lag_df,
          vars = variables_label,
          n_warning_events = length(fold_warning_messages),
          warning_types = paste(
            unique(fold_warning_types),
            collapse = " | "
          ),
          n_knot_warning_events =
            as.integer(fold_type_counts[["knot"]]),
          n_convergence_warning_events =
            as.integer(fold_type_counts[["convergence"]]),
          n_hessian_warning_events =
            as.integer(fold_type_counts[["hessian"]]),
          n_other_warning_events =
            as.integer(fold_type_counts[["other"]]),
          warning = paste(
            unique(fold_warning_messages),
            collapse = " | "
          ),
          stringsAsFactors = FALSE
        )
      }

      if (!isTRUE(fold_result$ok)) {

        candidate_failures[[fold_index]] <- data.frame(
          model_id = rep(
            model_id,
            length(test_groups)
          ),
          fold = rep(
            fold_index,
            length(test_groups)
          ),
          group = as.character(
            test_groups
          ),
          df_var = rep(
            exposure_df,
            length(test_groups)
          ),
          df_lag = rep(
            lag_df,
            length(test_groups)
          ),
          vars = rep(
            variables_label,
            length(test_groups)
          ),
          error = rep(
            fold_result$error,
            length(test_groups)
          ),
          stringsAsFactors = FALSE
        )

        next
      }

      group_position <- match(
        fold_result$groups,
        group_ids
      )

      if (anyNA(group_position) ||
          anyDuplicated(group_position)) {
        candidate_failures[[fold_index]] <- data.frame(
          model_id = rep(
            model_id,
            length(test_groups)
          ),
          fold = rep(
            fold_index,
            length(test_groups)
          ),
          group = as.character(
            test_groups
          ),
          df_var = rep(
            exposure_df,
            length(test_groups)
          ),
          df_lag = rep(
            lag_df,
            length(test_groups)
          ),
          vars = rep(
            variables_label,
            length(test_groups)
          ),
          error = rep(
            "Internal validation error: duplicate or unknown held-out group.",
            length(test_groups)
          ),
          stringsAsFactors = FALSE
        )

        next
      }

      valid_prediction <- is.finite(
        fold_result$obs
      ) & is.finite(
        fold_result$pred
      )

      if (any(valid_prediction)) {

        valid_position <- group_position[
          valid_prediction
        ]

        observed[valid_position] <-
          fold_result$obs[valid_prediction]

        predicted[valid_position] <-
          fold_result$pred[valid_prediction]

        success[valid_position] <- TRUE

        prediction_row <- data.frame(
          model_id = rep(
            model_id,
            sum(valid_prediction)
          ),
          fold = rep(
            fold_index,
            sum(valid_prediction)
          ),
          group = as.character(
            fold_result$groups[valid_prediction]
          ),
          df_var = rep(
            exposure_df,
            sum(valid_prediction)
          ),
          df_lag = rep(
            lag_df,
            sum(valid_prediction)
          ),
          vars = rep(
            variables_label,
            sum(valid_prediction)
          ),
          observed = fold_result$obs[valid_prediction],
          predicted = fold_result$pred[valid_prediction],
          stringsAsFactors = FALSE
        )

        # For binomial outcomes, predicted remains the out-of-sample validation probability.
        # Thresholded classes are stored separately.
        if (identical(outcome_type, "binary")) {
          prediction_row$predicted_class <-
            as.integer(
              prediction_row$predicted >= threshold
            )
        }

        candidate_predictions[[fold_index]] <-
          prediction_row
      }

      if (any(!valid_prediction)) {

        invalid_groups <- fold_result$groups[
          !valid_prediction
        ]

        candidate_failures[[fold_index]] <- data.frame(
          model_id = rep(
            model_id,
            length(invalid_groups)
          ),
          fold = rep(
            fold_index,
            length(invalid_groups)
          ),
          group = as.character(
            invalid_groups
          ),
          df_var = rep(
            exposure_df,
            length(invalid_groups)
          ),
          df_lag = rep(
            lag_df,
            length(invalid_groups)
          ),
          vars = rep(
            variables_label,
            length(invalid_groups)
          ),
          error = rep(
            "Held-out prediction or observed value was non-finite.",
            length(invalid_groups)
          ),
          stringsAsFactors = FALSE
        )
      }

      # A fold is successful only when all groups assigned to that fold produced
      # valid out-of-sample validation predictions.
      fold_success[fold_index] <-
        length(valid_prediction) == length(test_groups) &&
        all(valid_prediction)
    }

    evaluation_mask <- group_ids %in% evaluation_groups
    successful_evaluation <- success & evaluation_mask

    n_predictions <- length(evaluation_groups)
    n_success <- sum(successful_evaluation)
    n_failed <- n_predictions - n_success

    n_success_folds <- sum(fold_success)
    n_failed_folds <- n_folds - n_success_folds

    n_warning_folds <- sum(warning_fold)
    n_warning_events <- sum(warning_events)
    warning_rate <- n_warning_folds / n_folds

    warning_folds_by_type <- colSums(
      warning_type_fold
    )

    warning_events_by_type <- colSums(
      warning_type_events
    )

    all_validation_warning_messages <- unlist(
      candidate_warning_messages,
      recursive = FALSE,
      use.names = FALSE
    )

    all_validation_warning_messages <- as.character(
      all_validation_warning_messages
    )

    unique_validation_warning_messages <- unique(
      all_validation_warning_messages
    )

    candidate_result <- NULL
    metric_warning_messages <- character(0)
    full_fit <- NULL
    full_fit_warning_messages <- character(0)

    full_fit_warning_events_by_type <- setNames(
      integer(
        length(warning_type_levels)
      ),
      warning_type_levels
    )

    # min_success is deliberately evaluated on successful GROUP predictions.
    if (n_success >= min_success) {

      metrics <- withCallingHandlers(
        .compute_performance_metrics(
          observed = observed[successful_evaluation],
          predicted = predicted[successful_evaluation],
          family = family_name,
          threshold = threshold
        ),
        warning = function(w) {
          metric_warning_messages <<- c(
            metric_warning_messages,
            conditionMessage(w)
          )
          invokeRestart("muffleWarning")
        }
      )

      metric_warning_messages <- unique(
        metric_warning_messages
      )

      if (length(metric_warning_messages)) {

        candidate_warnings[[n_folds + 1L]] <- data.frame(
          model_id = model_id,
          fold = NA_integer_,
          group = NA_character_,
          stage = "metrics",
          df_var = exposure_df,
          df_lag = lag_df,
          vars = variables_label,
          n_warning_events =
            length(metric_warning_messages),
          warning_types = "other",
          n_knot_warning_events = 0L,
          n_convergence_warning_events = 0L,
          n_hessian_warning_events = 0L,
          n_other_warning_events =
            length(metric_warning_messages),
          warning = paste(
            metric_warning_messages,
            collapse = " | "
          ),
          stringsAsFactors = FALSE
        )
      }

      candidate_identity <- data.frame(
        model_id = model_id,
        df_var = exposure_df,
        df_lag = lag_df,
        vars = variables_label,
        n_vars = length(variables),
        stringsAsFactors = FALSE
      )

      candidate_diagnostics <- data.frame(
        n_groups_total = n_groups,
        n_folds = n_folds,
        n_success_folds = n_success_folds,
        n_failed_folds = n_failed_folds,
        n_predictions = n_predictions,
        n_success = n_success,
        n_failed = n_failed,
        n_warning_folds = n_warning_folds,
        warning_rate = warning_rate,
        n_warning_events = n_warning_events,
        n_knot_warning_folds =
          as.integer(warning_folds_by_type[["knot"]]),
        n_convergence_warning_folds =
          as.integer(warning_folds_by_type[["convergence"]]),
        n_hessian_warning_folds =
          as.integer(warning_folds_by_type[["hessian"]]),
        n_other_warning_folds =
          as.integer(warning_folds_by_type[["other"]]),
        n_knot_warning_events =
          as.integer(warning_events_by_type[["knot"]]),
        n_convergence_warning_events =
          as.integer(warning_events_by_type[["convergence"]]),
        n_hessian_warning_events =
          as.integer(warning_events_by_type[["hessian"]]),
        n_other_warning_events =
          as.integer(warning_events_by_type[["other"]]),
        stringsAsFactors = FALSE
      )

      candidate_result <- cbind(
        candidate_identity,
        metrics,
        candidate_diagnostics
      )

      # Preserve historical behavior: every eligible candidate is refitted on
      # all data when keep_fits = TRUE. Ranking/top_n filtering occurs later.
      if (keep_fits) {

        full_fit_result <- tryCatch(
          withCallingHandlers(
            {
              list(
                ok = TRUE,
                fit = fit_candidate(
                  build_candidate(
                    data_long,
                    variables,
                    exposure_df,
                    lag_df
                  )
                )
              )
            },
            warning = function(w) {
              full_fit_warning_messages <<- c(
                full_fit_warning_messages,
                conditionMessage(w)
              )
              invokeRestart("muffleWarning")
            }
          ),
          error = function(e) {
            list(
              ok = FALSE,
              error = conditionMessage(e)
            )
          }
        )

        if (length(full_fit_warning_messages)) {

          full_fit_warning_types <- vapply(
            full_fit_warning_messages,
            classify_warning,
            character(1)
          )

          full_fit_type_counts <- table(
            factor(
              full_fit_warning_types,
              levels = warning_type_levels
            )
          )

          full_fit_warning_events_by_type[] <-
            as.integer(full_fit_type_counts)

          candidate_warnings[[n_folds + 2L]] <- data.frame(
            model_id = model_id,
            fold = NA_integer_,
            group = NA_character_,
            stage = "full_refit",
            df_var = exposure_df,
            df_lag = lag_df,
            vars = variables_label,
            n_warning_events =
              length(full_fit_warning_messages),
            warning_types = paste(
              unique(full_fit_warning_types),
              collapse = " | "
            ),
            n_knot_warning_events =
              as.integer(full_fit_type_counts[["knot"]]),
            n_convergence_warning_events =
              as.integer(full_fit_type_counts[["convergence"]]),
            n_hessian_warning_events =
              as.integer(full_fit_type_counts[["hessian"]]),
            n_other_warning_events =
              as.integer(full_fit_type_counts[["other"]]),
            warning = paste(
              unique(full_fit_warning_messages),
              collapse = " | "
            ),
            stringsAsFactors = FALSE
          )
        }

        if (isTRUE(full_fit_result$ok)) {
          full_fit <- full_fit_result$fit
        } else {
          candidate_failures[[n_folds + 1L]] <- data.frame(
            model_id = model_id,
            fold = NA_integer_,
            group = NA_character_,
            df_var = exposure_df,
            df_lag = lag_df,
            vars = variables_label,
            error = paste0(
              "Full-data refit failed: ",
              full_fit_result$error
            ),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    warning_summary <- data.frame(
      model_id = model_id,
      df_var = exposure_df,
      df_lag = lag_df,
      vars = variables_label,
      n_groups_total = n_groups,
      n_folds = n_folds,
      n_success_folds = n_success_folds,
      n_failed_folds = n_failed_folds,
      n_predictions = n_predictions,
      n_success = n_success,
      n_failed = n_failed,
      n_warning_folds = n_warning_folds,
      warning_rate = warning_rate,
      n_warning_events = n_warning_events,
      n_knot_warning_folds =
        as.integer(warning_folds_by_type[["knot"]]),
      n_convergence_warning_folds =
        as.integer(warning_folds_by_type[["convergence"]]),
      n_hessian_warning_folds =
        as.integer(warning_folds_by_type[["hessian"]]),
      n_other_warning_folds =
        as.integer(warning_folds_by_type[["other"]]),
      n_knot_warning_events =
        as.integer(warning_events_by_type[["knot"]]),
      n_convergence_warning_events =
        as.integer(warning_events_by_type[["convergence"]]),
      n_hessian_warning_events =
        as.integer(warning_events_by_type[["hessian"]]),
      n_other_warning_events =
        as.integer(warning_events_by_type[["other"]]),
      n_warning_types =
        length(unique_validation_warning_messages),
      n_warning_classes =
        length(
          unique(
            vapply(
              unique_validation_warning_messages,
              classify_warning,
              character(1)
            )
          )
        ),
      warning_classes = if (
        length(unique_validation_warning_messages)
      ) {
        paste(
          unique(
            vapply(
              unique_validation_warning_messages,
              classify_warning,
              character(1)
            )
          ),
          collapse = " | "
        )
      } else {
        NA_character_
      },
      warning_messages = if (
        length(unique_validation_warning_messages)
      ) {
        paste(
          unique_validation_warning_messages,
          collapse = " | "
        )
      } else {
        NA_character_
      },
      metric_warning_events =
        length(metric_warning_messages),
      metric_warning_messages = if (
        length(metric_warning_messages)
      ) {
        paste(
          metric_warning_messages,
          collapse = " | "
        )
      } else {
        NA_character_
      },
      full_refit_warning_events =
        length(full_fit_warning_messages),
      full_refit_knot_warning_events =
        as.integer(
          full_fit_warning_events_by_type[["knot"]]
        ),
      full_refit_convergence_warning_events =
        as.integer(
          full_fit_warning_events_by_type[["convergence"]]
        ),
      full_refit_hessian_warning_events =
        as.integer(
          full_fit_warning_events_by_type[["hessian"]]
        ),
      full_refit_other_warning_events =
        as.integer(
          full_fit_warning_events_by_type[["other"]]
        ),
      full_refit_warning_classes = if (
        length(full_fit_warning_messages)
      ) {
        paste(
          unique(
            vapply(
              full_fit_warning_messages,
              classify_warning,
              character(1)
            )
          ),
          collapse = " | "
        )
      } else {
        NA_character_
      },
      full_refit_warning_messages = if (
        length(full_fit_warning_messages)
      ) {
        paste(
          unique(full_fit_warning_messages),
          collapse = " | "
        )
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )

    candidate_elapsed <-
      proc.time()[["elapsed"]] -
      candidate_start_time

    list(
      model_id = model_id,
      result = candidate_result,
      predictions = Filter(
        Negate(is.null),
        candidate_predictions
      ),
      failures = Filter(
        Negate(is.null),
        candidate_failures
      ),
      warnings = Filter(
        Negate(is.null),
        candidate_warnings
      ),
      warning_summary = warning_summary,
      fit = full_fit,
      elapsed_seconds =
        as.numeric(candidate_elapsed)
    )
  }

  # Package functions should not change the user's global future::plan().
  # Therefore parallelism is used whenever the caller has configured >1 worker;
  # otherwise the exact same evaluator runs sequentially.
  use_future <-
    total_candidates > 1L &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE) &&
    future::nbrOfWorkers() > 1L

  if (use_future) {
    if (verbose) {
      message(
        "Evaluating ", total_candidates,
        " candidate models across ", future::nbrOfWorkers(),
        " future workers; validation iterations remain sequential within each ",
        "candidate."
      )
    }

    candidate_outputs <- future.apply::future_lapply(
      X = seq_len(total_candidates),
      FUN = evaluate_candidate,
      future.seed = TRUE,
      # One future per candidate: validation iterations remain sequential inside each future.
      future.scheduling = Inf
    )
  } else {
    if (verbose && total_candidates > 1L) {
      message(
        "Evaluating candidates sequentially. To enable candidate-level ",
        "parallelism, set a future plan with >1 worker before calling ",
        "find_bestfit()."
      )
    }
    candidate_outputs <- lapply(
      seq_len(total_candidates),
      evaluate_candidate
    )
  }

  results <- Filter(
    Negate(is.null),
    lapply(candidate_outputs, function(x) x$result)
  )
  predictions <- unlist(
    lapply(candidate_outputs, function(x) x$predictions),
    recursive = FALSE,
    use.names = FALSE
  )
  failures <- unlist(
    lapply(candidate_outputs, function(x) x$failures),
    recursive = FALSE,
    use.names = FALSE
  )
  warnings <- unlist(
    lapply(candidate_outputs, function(x) x$warnings),
    recursive = FALSE,
    use.names = FALSE
  )

  warning_summary_data <- do.call(
    rbind,
    lapply(candidate_outputs, function(x) x$warning_summary)
  )
  rownames(warning_summary_data) <- NULL

  candidate_times <- data.frame(
    model_id = vapply(
      candidate_outputs,
      function(x) as.integer(x$model_id),
      integer(1)
    ),
    elapsed_seconds = vapply(
      candidate_outputs,
      function(x) as.numeric(x$elapsed_seconds),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  # Warnings have already been muffled inside each validation iteration. Emit at most one
  # concise summary per affected candidate, in deterministic model_id order.
  if (verbose && nrow(warning_summary_data)) {
    affected <- warning_summary_data$n_warning_folds > 0L |
      warning_summary_data$metric_warning_events > 0L |
      warning_summary_data$full_refit_warning_events > 0L

    for (i in which(affected)) {
      ws <- warning_summary_data[i, , drop = FALSE]
      validation_warning_text <- if (ws$n_warning_folds > 0L) {
        if (identical(validation_method, "holdout")) {
          paste0(
            ws$n_warning_folds, "/", ws$n_folds,
            " holdout validation iteration(s) (",
            ws$n_warning_events, " warning event(s))"
          )
        } else {
          paste0(
            ws$n_warning_folds, "/", ws$n_folds,
            " validation fold(s) (",
            ws$n_warning_events, " warning event(s))"
          )
        }
      } else if (identical(validation_method, "holdout")) {
        "0 holdout validation iterations"
      } else {
        "0 validation folds"
      }

      type_fold_counts <- c(
        knot = ws$n_knot_warning_folds,
        convergence = ws$n_convergence_warning_folds,
        hessian = ws$n_hessian_warning_folds,
        other = ws$n_other_warning_folds
      )
      type_fold_counts <- type_fold_counts[type_fold_counts > 0L]
      type_text <- if (length(type_fold_counts)) {
        paste0(
          "; warning classes (affected validation iterations): ",
          paste(
            paste0(names(type_fold_counts), "=", as.integer(type_fold_counts)),
            collapse = ", "
          )
        )
      } else {
        ""
      }

      metric_text <- if (ws$metric_warning_events > 0L) {
        paste0(
          "; performance metrics: ", ws$metric_warning_events,
          " warning event(s)"
        )
      } else {
        ""
      }

      full_fit_text <- if (ws$full_refit_warning_events > 0L) {
        paste0(
          "; full-data refit: ", ws$full_refit_warning_events,
          " warning event(s)"
        )
      } else {
        ""
      }

      message(
        "Warning summary [model ", ws$model_id, "/", total_candidates,
        "] df_var = ", ws$df_var,
        ", df_lag = ", ws$df_lag,
        ", vars = ", ws$vars,
        ": ", validation_warning_text, type_text, metric_text, full_fit_text,
        ". See attr(result, \"warning_summary\") and attr(result, \"warnings\") for details."
      )
    }
  }

  fits_list <- list()
  if (keep_fits) {
    for (candidate_output in candidate_outputs) {
      if (!is.null(candidate_output$fit)) {
        fits_list[[as.character(candidate_output$model_id)]] <- candidate_output$fit
      }
    }
  }

  if (!length(results)) {
    stop(
      "No candidate model produced at least `min_success` successful ",
      "group-level out-of-sample validation predictions.",
      call. = FALSE
    )
  }
  results_data <- do.call(rbind, results)

  # Rank primarily by the user-selected metric. Remaining family-appropriate
  # metrics are deterministic tie-breakers, each using its own known direction.
  ranking_metrics <- c(
    rank_metric,
    setdiff(available_metrics, rank_metric)
  )

  if (all(is.na(results_data[[rank_metric]]))) {
    warning(
      "`rank_metric = ", rank_metric,
      "` is NA for all eligible candidates. Ranking will therefore be ",
      "resolved by the remaining available metrics in canonical order.",
      call. = FALSE
    )
  }

  ranking_vectors <- lapply(
    ranking_metrics,
    function(metric_name) {
      values <- results_data[[metric_name]]
      if (identical(.metric_direction(metric_name), "maximize")) {
        -values
      } else {
        values
      }
    }
  )

  ranking_order <- do.call(
    order,
    c(
      ranking_vectors,
      list(results_data$model_id),
      list(na.last = TRUE)
    )
  )

  results_data <- results_data[ranking_order, , drop = FALSE]
  results_data$rank <- seq_len(nrow(results_data))
  results_data <- results_data[, c("rank", setdiff(names(results_data), "rank")), drop = FALSE]
  if (is.finite(top_n)) results_data <- utils::head(results_data, top_n)
  rownames(results_data) <- NULL

  prediction_data <- if (length(predictions)) do.call(rbind, predictions) else data.frame()
  failure_data <- if (length(failures)) do.call(rbind, failures) else data.frame()
  warning_data <- if (length(warnings)) do.call(rbind, warnings) else data.frame()

  # Preserve the original flat prediction data frame and additionally provide
  # a convenient list split by model_id for direct model-specific extraction.
  predictions_by_model <- if (nrow(prediction_data)) {
    split_predictions <- split(
      prediction_data,
      as.character(prediction_data$model_id),
      drop = TRUE
    )
    lapply(split_predictions, function(x) {
      rownames(x) <- NULL
      x
    })
  } else {
    list()
  }

  total_elapsed_seconds <- as.numeric(
    proc.time()[["elapsed"]] - function_start_time
  )

  format_elapsed <- function(seconds) {
    seconds <- max(0, as.numeric(seconds))
    hours <- floor(seconds / 3600)
    minutes <- floor((seconds %% 3600) / 60)
    secs <- seconds %% 60
    sprintf("%02d:%02d:%05.2f", as.integer(hours), as.integer(minutes), secs)
  }

  timing_data <- data.frame(
    total_elapsed_seconds = total_elapsed_seconds,
    total_elapsed = format_elapsed(total_elapsed_seconds),
    total_candidates = total_candidates,
    validation_method = validation_method,
    validation_scheme = validation_scheme,
    n_groups_total = n_groups,
    n_evaluation_groups = n_evaluation_groups,
    n_training_groups = if (identical(validation_method, "holdout")) {
      n_training_groups
    } else {
      NA_integer_
    },
    n_test_groups = if (identical(validation_method, "holdout")) {
      n_test_groups
    } else {
      NA_integer_
    },
    test_prop = if (identical(validation_method, "holdout")) {
      effective_test_prop
    } else {
      NA_real_
    },
    n_folds = n_folds,
    parallel = use_future,
    workers = if (use_future) future::nbrOfWorkers() else 1L,
    stringsAsFactors = FALSE
  )

  attr(results_data, "predictions") <- prediction_data
  attr(results_data, "predictions_by_model") <- predictions_by_model
  attr(results_data, "failures") <- failure_data
  attr(results_data, "warnings") <- warning_data
  attr(results_data, "warning_summary") <- warning_summary_data
  attr(results_data, "timing") <- timing_data
  attr(results_data, "candidate_times") <- candidate_times
  attr(results_data, "family") <- family_name
  attr(results_data, "outcome_type") <- outcome_type
  attr(results_data, "rank_metric") <- rank_metric
  attr(results_data, "threshold") <- threshold
  attr(results_data, "prediction_level") <- "population"
  attr(results_data, "prediction_estimand") <- "expected_response"
  attr(results_data, "prediction_contract") <- "central_expected_response"
  attr(results_data, "random_effect") <- random_effect
  attr(results_data, "random_effect_model") <- if (!is.null(random_effect) &&
                                                   model_engine %in% c("inla","bdlnm")) {"iid"} else {NA_character_}
  attr(results_data, "random_effect_prior") <- random_effect_prior
  attr(results_data, "spatial_effect") <- spatial_effect
  attr(results_data, "spatial_structure") <- spatial_structure
  attr(results_data, "spatial_group") <- spatial_group
  attr(results_data, "spatial_term") <- spatial_term
  attr(results_data, "has_spatial_effect") <- !is.null(spatial_effect)
  attr(results_data, "validation_method") <- validation_method
  attr(results_data, "validation_scheme") <- validation_scheme
  attr(results_data, "k") <- effective_k
  attr(results_data, "validation_seed") <- if (
    validation_method %in% c("k-fold", "holdout") && !is.null(seed)
  ) {
    seed
  } else {
    NA_integer_
  }
  attr(results_data, "validation_stratified") <- validation_stratified
  attr(results_data, "test_prop") <- if (
    identical(validation_method, "holdout")
  ) {
    effective_test_prop
  } else {
    NA_real_
  }
  attr(results_data, "n_groups_total") <- n_groups
  attr(results_data, "n_evaluation_groups") <- n_evaluation_groups
  attr(results_data, "training_groups") <- training_groups
  attr(results_data, "test_groups") <- test_groups
  attr(results_data, "evaluation_groups") <- evaluation_groups
  attr(results_data, "validation_assignments") <- validation_assignments
  attr(results_data, "validation_balance") <- validation_balance
  attr(results_data, "retained_fit_scope") <- if (keep_fits) {
    "full_data_refit_after_validation"
  } else {
    NA_character_
  }
  attr(results_data, "basis_training_only") <- TRUE
  attr(results_data, "max_lag") <- max_lag
  attr(results_data, "history_length") <- required_history_length
  attr(results_data, "history_contract") <- history_contract
  attr(results_data, "time_step") <- time_step

  if (keep_fits) {
    retained_ids <- as.character(results_data$model_id)
    attr(results_data, "fits") <- fits_list[names(fits_list) %in% retained_ids]
  }

  if (verbose) {
    message(
      "find_bestfit completed in ",
      timing_data$total_elapsed,
      " (", format(round(total_elapsed_seconds, 2), nsmall = 2),
      " s). Validation = ", validation_method,
      if (identical(validation_method, "k-fold")) {
        paste0(" (k=", effective_k, ")")
      } else if (identical(validation_method, "holdout")) {
        paste0(" (test_prop=", format(effective_test_prop), ")")
      } else {
        ""
      },
      "; ranked by ", rank_metric, " (", rank_direction, ")."
    )
  }

  results_data
}
