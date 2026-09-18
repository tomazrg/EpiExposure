#' Find the best DLNM model structure using grouped cross-validation
#'
#' Evaluates combinations of exposure variables and cross-basis dimensions using
#' grouped cross-validation, then ranks candidate models with family-aware
#' predictive-performance metrics.
#'
#' Two cross-validation schemes are available:
#'
#' - `"LOOCV"`: leave one complete group out at a time;
#' - `"k-fold"`: divide complete groups into `k` approximately equal folds.
#'
#' The unit of cross-validation is always the complete `group`. Exposure-history
#' rows belonging to one group are never split between training and test data.
#' For binomial outcomes, grouped k-fold allocation is stratified so the numbers
#' of outcome-0 and outcome-1 groups are distributed as evenly as possible
#' across folds.
#'
#' The statistical target used for every held-out prediction is the same target
#' used by the EpiExposure v1 prediction layer: the population/fixed-component
#' expected response, with fitted random effects set to zero. Cross-validation
#' ranking uses deterministic predictions from the harmonized central parameter
#' estimate and does not propagate coefficient, posterior, residual, or
#' future-observation uncertainty.
#'
#' @param data Long-format data frame containing the response, grouping column,
#'   chronological time column, and candidate exposure variables.
#' @param response Character scalar naming the response column in `data`. The
#'   response may be repeated over exposure-history rows but must be constant
#'   within each cross-validation group because one outcome is predicted per
#'   group.
#' @param group Character scalar naming the independent cross-validation unit.
#'   Every row belonging to one group is kept together in either training or
#'   test data within a fold.
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
#'   before cross-validation and passed unchanged to `fit_epidlnm()`.
#'   To specify a non-default link, supply a supported family object,
#'   such as family = stats::Gamma(link = "inverse");
#'   character family names use the EpiExposure default link.
#' @param random_effect Optional character scalar naming one grouping column used
#'   as a random intercept by `fit_epidlnm()`. It must be constant within each
#'   cross-validation group. Random effects may contribute to model fitting, but
#'   held-out predictions are always population-level and therefore set fitted
#'   random effects to zero. `glm` and `gls` do not accept `random_effect`;
#'   `gamm` requires it under the current EpiExposure v1 fitting contract.
#' @param min_success Positive integer of at least 2 giving the minimum number of
#'   **groups with a finite out-of-fold prediction** required for a candidate to
#'   be eligible for metric calculation and, when requested, full-data refitting.
#'
#'   `min_success` counts successful group-level predictions, not successful
#'   folds. Under `"LOOCV"` these quantities coincide because each fold contains
#'   one held-out group. Under `"k-fold"` one fold contains multiple groups, so a
#'   failed fold can remove several out-of-fold predictions at once.
#'
#'   `min_success` is a technical eligibility threshold, not a recommendation
#'   that large numbers of failed predictions are acceptable. Candidate output
#'   also reports `n_success`, `n_failed`, `n_success_folds`, and
#'   `n_failed_folds` so cross-validation completeness can be inspected.
#' @param rank_metric Optional character scalar selecting the primary ranking
#'   metric. If `NULL`, the default is `"CCC"` for non-binary outcomes and
#'   `"ROC_AUC"` for binary outcomes. Non-binary metrics are `"CCC"`, `"Cb"`,
#'   `"rho"`, `"RMSE"`, and `"MAE"`. Binary metrics are `"ROC_AUC"`, `"Brier"`,
#'   `"LogLoss"`, `"Accuracy"`, `"Balanced_Accuracy"`, `"Sensitivity"`,
#'   `"Specificity"`, `"F1"`, `"MCC"`, and `"Precision"`. Metric direction
#'   (maximize/minimize) is resolved automatically.
#' @param threshold Numeric probability strictly between 0 and 1 used only for
#'   binomial classification metrics. The default is `0.5`. Predicted
#'   probabilities are retained; thresholding is used only to create the
#'   auxiliary 0/1 class required by classification metrics. ROC AUC, Brier
#'   score, and Log Loss do not depend on the threshold.
#' @param top_n Positive integer or `Inf`; number of ranked candidates returned.
#' @param keep_fits Logical. If `TRUE`, every candidate that reaches
#'   `min_success` is refitted on the complete data using a full-data basis.
#'   After ranking and `top_n` filtering, only retained candidate fits are stored
#'   in `attr(result, "fits")`. This preserves the historical behavior of
#'   `find_bestfit()`.
#' @param verbose Logical. Print progress and concise warning summaries.
#' @param cv_method Character. Cross-validation scheme: `"LOOCV"` (default) or
#'   `"k-fold"`. `"LOOCV"` leaves one complete group out at a time. `"k-fold"`
#'   assigns complete groups to `k` folds and never splits exposure-history rows
#'   from one group across folds.
#' @param k Positive integer number of folds used only when
#'   `cv_method = "k-fold"`. It must satisfy `2 <= k < number of groups`.
#'
#'   Fold sizes are made as equal as possible; therefore the number of groups in
#'   any two folds differs by at most one. For binomial outcomes, allocation is
#'   additionally stratified by the group-level 0/1 response. Each outcome class
#'   must contain at least `k` groups so every fold can contain both classes.
#' @param seed Optional non-negative integer controlling only the random
#'   assignment of groups to folds when `cv_method = "k-fold"`. Supplying a seed
#'   makes `attr(result, "fold_assignments")` reproducible. The caller's existing
#'   global random-number state is restored after fold construction. `seed` does
#'   not alter LOOCV fold membership and is not passed to `fit_epidlnm()`.
#' @param ... Named additional arguments passed to `fit_epidlnm()`. Core
#'   arguments managed by `find_bestfit()` (`data`, `model_engine`, `family`,
#'   `random_effect`, `epiexposure_spec`, and `basis_objects`) cannot be supplied
#'   again through `...`. For INLA-backed engines, likelihood availability,
#'   `control.family$control.link$model` consistency, and
#'   `control.compute$config = TRUE` compatibility are checked before CV.
#'   `fit_epidlnm()` then supplies required INLA controls during each fit.
#'
#' @return A data frame ranked by `rank_metric`. Non-binary families report
#'   `CCC`, `Cb`, `rho`, `RMSE`, and `MAE`. Binomial models report `ROC_AUC`,
#'   `Brier`, `LogLoss`, `Accuracy`, `Balanced_Accuracy`, `Sensitivity`,
#'   `Specificity`, `F1`, `MCC`, and `Precision`.
#'
#'   Candidate diagnostics distinguish prediction completeness from fold
#'   completeness. `n_success` and `n_failed` count held-out groups, whereas
#'   `n_success_folds` and `n_failed_folds` count folds. A fold is considered
#'   successful only when every group assigned to that fold receives a finite
#'   out-of-fold prediction.
#'
#'   The result retains the diagnostic attributes used by earlier versions:
#'   `"predictions"`, `"predictions_by_model"`, `"failures"`, `"warnings"`,
#'   `"warning_summary"`, `"timing"`, and `"candidate_times"`. With
#'   `keep_fits = TRUE`, `"fits"` contains retained full-data refits.
#'
#'   Cross-validation metadata include:
#'
#'   - `"cv_method"`: `"LOOCV"` or `"k-fold"`;
#'   - `"cv_scheme"`: normalized descriptive scheme;
#'   - `"k"`: effective number of folds;
#'   - `"cv_seed"`: supplied k-fold seed or `NA`;
#'   - `"cv_stratified"`: whether binomial stratification was used;
#'   - `"fold_assignments"`: data frame with `group`, `fold`, and `observed`,
#'     giving the single fold assigned to every independent group;
#'   - `"fold_balance"`: number of groups per fold and, for binomial outcomes,
#'     the numbers of outcome-0 and outcome-1 groups.
#'
#'   Standard downstream metadata are stored in `"family"`, `"outcome_type"`,
#'   `"rank_metric"`, and `"threshold"`. Additional attributes document the
#'   prediction contract: `"prediction_level" = "population"`,
#'   `"prediction_estimand" = "expected_response"`,
#'   `"prediction_contract" = "central_expected_response"`,
#'   `"basis_training_only" = TRUE`, `"max_lag"`, `"history_length"`,
#'   `"history_contract"`, and `"time_step"`.
#'
#' @details
#' ## Grouped cross-validation
#'
#' Fold membership is constructed once before candidate evaluation and reused
#' unchanged for every candidate model. This ensures that candidate metrics are
#' compared on exactly the same training/test partitions.
#'
#' With `cv_method = "LOOCV"`, a data set containing \eqn{G} groups produces
#' \eqn{G} folds, each containing one test group.
#'
#' With `cv_method = "k-fold"`, complete groups are randomized into `k`
#' approximately equal folds. The total number of groups does not need to be
#' divisible by `k`; for example, 521 groups with `k = 5` produce fold sizes
#' 105, 104, 104, 104, and 104 in some fold order.
#'
#' For binomial outcomes, k-fold allocation is grouped and stratified. Outcome-0
#' groups and outcome-1 groups are each distributed as evenly as possible while
#' maintaining overall fold balance. Because every fold is required to contain
#' both classes, `k` cannot exceed the number of groups in the less frequent
#' class.
#'
#' Engine, family, and link compatibility is validated before candidate
#' evaluation. Unsupported combinations fail immediately rather than producing
#' repeated fitting failures across folds.
#'
#' ## Training-only basis construction
#'
#' For every fold, `define_exposures()` is called **only on the training groups**.
#' The effective `argvar`, `arglag`, lag range, spline knots, and boundary knots
#' stored in those returned training cross-basis objects are then reused to
#' transform every held-out exposure history. Held-out values therefore never
#' determine the training basis parameterization.
#'
#' For predictors with highly skewed distributions or many repeated values,
#' `splines::ns()` may issue a knot-placement warning when interior knots
#' coincide with boundary values. This adjustment is handled automatically
#' and does not prevent model fitting.
#'
#' Reusing the training basis does not mean reusing the same numerical
#' cross-basis matrix. Held-out exposure values produce new matrix values, but
#' they are transformed with the **same training basis parameterization**.
#'
#' The function validates basis transport explicitly. For every candidate, the
#' requested `max_lag`, the lag stored in the training cross-basis, and the lag
#' stored in the exposure specification must agree. Reconstructed training and
#' held-out cross-bases must have the same number of columns and, when native
#' cross-basis column names are available, the same native column names/order as
#' the training template. Canonical EpiExposure columns
#' `cb_<variable>_<index>` are then assigned deterministically.
#'
#' For `model_engine = "bdlnm"`, an additional epidemic-level matrix-form
#' cross-basis is built from those same training parameters so the cross-basis
#' included in the `bdlnm` formula has exactly one row per epidemic-level
#' outcome. Its numerical values are checked against the canonical candidate
#' design before fitting.
#'
#' The design-matrix attribute `cb_templates` stores the **original training
#' templates**, not cross-bases reconstructed from the first epidemic. The
#' effective specification passed to `fit_epidlnm()` is synchronized with the
#' returned training-template attributes.
#'
#' ## Temporal requirements
#'
#' A vector supplied to `dlnm::crossbasis()` represents one complete, ordered,
#' equally spaced exposure history. Accordingly, this function rejects
#' duplicated or irregular time values, requires the same time step across
#' groups, and requires **exactly `max_lag + 1` observations in every group**.
#'
#' Histories with fewer observations are rejected because the fitted lag window
#' is incomplete. Histories with more observations are also rejected: they are
#' not truncated to a trailing window and the function never chooses silently
#' which observations define the epidemiological history.
#'
#' All candidate exposure variables are columns of these same validated
#' long-format rows and must contain only finite values. Therefore every
#' candidate variable uses exactly the same number of time points and the same
#' temporal positions within each group.
#'
#' ## Held-out prediction target
#'
#' Every successful held-out group receives one deterministic out-of-fold
#' prediction of the expected response from the harmonized central
#' fixed/population parameter estimate. Group-specific random effects are set to
#' zero, including when a random intercept was fitted and held-out groups are
#' unseen. Bayesian engines therefore use posterior-mean fixed parameters for
#' deterministic CV prediction; candidates are not ranked using medians of
#' posterior expected predictions.
#'
#' Coefficient/posterior uncertainty and residual/future-observation noise are
#' intentionally excluded from model ranking. `find_bestfit()` evaluates point
#' predictive performance; uncertainty belongs to downstream prediction and
#' effect functions.
#'
#' ## Performance, warnings, and ranking
#'
#' Performance metrics are computed after all successful out-of-fold group
#' predictions for a candidate have been collected. Both `find_bestfit()` and
#' `ensemble_bestfit()` use `.compute_performance_metrics()` so metric
#' definitions are shared across the package.
#'
#' For non-binary outcomes (`beta`, `poisson`, `gamma`, `gaussian`, and
#' `negative_binomial`), performance is summarized by CCC, Cb, Pearson
#' correlation, RMSE, and MAE. For binomial outcomes, out-of-fold probabilities
#' are retained for ROC AUC, Brier score, and Log Loss; `threshold` is applied
#' only when classification metrics are calculated.
#'
#' Warnings raised during fold fitting/prediction are captured and classified as
#' `"knot"`, `"convergence"`, `"hessian"`, or `"other"`. Warning occurrence does
#' not by itself make a fold fail. Performance-metric warnings are recorded at
#' candidate level. When `verbose = TRUE`, only concise candidate-level warning
#' summaries are emitted after evaluation.
#'
#' Candidate ranking uses `rank_metric` first. Remaining family-appropriate
#' metrics are deterministic tie-breakers in canonical order, each with its own
#' maximize/minimize direction.
#'
#' Candidate-level parallelism is preserved: when the caller has configured a
#' `future` plan with more than one worker, candidates are evaluated in parallel
#' while all cross-validation folds belonging to one candidate remain
#' sequential.
#'
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
    min_success = 2,
    rank_metric = NULL,
    threshold = 0.5,
    top_n = Inf,
    keep_fits = FALSE,
    verbose = TRUE,
    cv_method = c("LOOCV", "k-fold"),
    k = 5,
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

  cv_method <- match.arg(
    cv_method,
    choices = c("LOOCV", "k-fold")
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

  if (identical(cv_method, "k-fold")) {
    if (!is_whole_scalar(k) || k < 2L) {
      stop(
        "`k` must be an integer greater than or equal to 2 when ",
        "`cv_method = 'k-fold'`.",
        call. = FALSE
      )
    }
    k <- as.integer(k)
  }

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

  # --------------------------------------------------------------------------
  # Basic columns and user inputs
  # --------------------------------------------------------------------------

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
    dot_names <- names(fit_dots)
    if (is.null(dot_names) || anyNA(dot_names) || any(!nzchar(dot_names))) {
      stop(
        "All arguments supplied through `...` must be explicitly named.",
        call. = FALSE
      )
    }

    reserved_fit_args <- c(
      "data", "model_engine", "family", "random_effect",
      "epiexposure_spec", "basis_objects"
    )
    duplicated_fit_args <- intersect(dot_names, reserved_fit_args)

    if (length(duplicated_fit_args)) {
      stop(
        "Argument(s) managed by `find_bestfit()` cannot be supplied again ",
        "through `...`: ",
        paste(duplicated_fit_args, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # --------------------------------------------------------------------------
  # Canonical family and early response/engine validation
  # --------------------------------------------------------------------------

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

  # Reject unsupported engine/family/link combinations BEFORE starting any CV
  # candidate. fit_epidlnm() repeats this validation before native fitting.
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

    # Validate conflicting INLA controls once, before creating CV folds. Do not
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

  # --------------------------------------------------------------------------
  # Canonical working data
  # --------------------------------------------------------------------------

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
      "cross-validation group. One group-level outcome is required.",
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

  # --------------------------------------------------------------------------
  # Temporal completeness and common lag unit
  # --------------------------------------------------------------------------

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
      "Every cross-validation group must contain exactly ",
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
  # be unique/constant inside each cross-validation group.
  if (!is.null(random_effect) && random_effect != group) {
    re_check <- unique(
      data_long[, c("epi_id", random_effect), drop = FALSE]
    )
    if (anyDuplicated(re_check$epi_id)) {
      stop(
        "`random_effect` must be constant within each cross-validation group because ",
        "the fitted design contains one row per group.",
        call. = FALSE
      )
    }
  }

  # --------------------------------------------------------------------------
  # Ranking metric contract
  # --------------------------------------------------------------------------

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
        ". It is applied only to out-of-fold predicted probabilities when ",
        "classification metrics are computed; ROC_AUC, Brier, and LogLoss ",
        "remain probability-scale metrics."
      )
    }
  }

  # --------------------------------------------------------------------------
  # Strict training-template transport
  # --------------------------------------------------------------------------

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

      basis_objects[[variable]] <- cb_epidemic
    }

    attr(basis_objects, "spec") <- template_info$spec
    basis_objects
  }

  .attach_random_effect_metadata <- function(design, source_data) {
    if (is.null(random_effect)) {
      return(design)
    }

    # When the random intercept is the same as the CV group, the canonical
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
        epiexposure_spec = candidate$spec,
        basis_objects = candidate$basis_objects
      ),
      fit_dots
    )

    fitted <- do.call(fit_epidlnm, fit_call)

    # Enforce the strict metadata/prediction contract immediately, rather than
    # discovering an outdated fit only when a held-out prediction is attempted.
    meta <- .get_epiexposure_metadata(fitted)

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

  # --------------------------------------------------------------------------
  # Candidate variable sets
  # --------------------------------------------------------------------------

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

  # --------------------------------------------------------------------------
  # Cross-validation fold construction
  # --------------------------------------------------------------------------

  group_ids <- as.character(response_check$epi_id)
  group_response <- as.numeric(response_check$y)
  n_groups <- length(group_ids)

  if (n_groups < 2L) {
    stop(
      "At least two independent groups are required for cross-validation.",
      call. = FALSE
    )
  }

  if (anyDuplicated(group_ids)) {
    stop(
      "Internal cross-validation error: group identifiers are not unique after ",
      "canonicalization.",
      call. = FALSE
    )
  }

  if (min_success > n_groups) {
    stop(
      "`min_success` cannot exceed the number of independent groups (",
      n_groups, ").",
      call. = FALSE
    )
  }

  cv_stratified <- FALSE

  if (identical(cv_method, "LOOCV")) {

    n_folds <- n_groups
    effective_k <- n_folds
    fold_assignment <- seq_len(n_groups)
    cv_scheme <- "leave_one_group_out"

  } else {

    if (k >= n_groups) {
      stop(
        "For `cv_method = 'k-fold'`, `k` must be smaller than the number of ",
        "independent groups. Received k = ", k, " and ", n_groups,
        " groups. Use `cv_method = 'LOOCV'` to leave one group out at a time.",
        call. = FALSE
      )
    }

    n_folds <- k
    effective_k <- k

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
          minority_count, "` or use `cv_method = 'LOOCV'`.",
          call. = FALSE
        )
      }

      cv_stratified <- TRUE
      cv_scheme <- "stratified_grouped_k_fold"

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

      cv_scheme <- "grouped_k_fold"

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
  }

  if (length(fold_assignment) != n_groups ||
      anyNA(fold_assignment) ||
      any(!fold_assignment %in% seq_len(n_folds))) {
    stop(
      "Internal cross-validation error: invalid fold assignment.",
      call. = FALSE
    )
  }

  fold_assignments <- data.frame(
    group = group_ids,
    fold = as.integer(fold_assignment),
    observed = group_response,
    stringsAsFactors = FALSE
  )

  fold_sizes <- tabulate(
    fold_assignments$fold,
    nbins = n_folds
  )

  if (any(fold_sizes < 1L)) {
    stop(
      "Internal cross-validation error: at least one fold contains no groups.",
      call. = FALSE
    )
  }

  if (identical(cv_method, "k-fold") &&
      max(fold_sizes) - min(fold_sizes) > 1L) {
    stop(
      "Internal cross-validation error: grouped k-fold allocation is not ",
      "balanced by number of groups.",
      call. = FALSE
    )
  }

  fold_balance <- data.frame(
    fold = seq_len(n_folds),
    n_groups = as.integer(fold_sizes),
    stringsAsFactors = FALSE
  )

  if (identical(outcome_type, "binary")) {
    fold_balance$n_outcome_0 <- vapply(
      seq_len(n_folds),
      function(fold_index) {
        sum(
          fold_assignments$fold == fold_index &
            fold_assignments$observed == 0
        )
      },
      integer(1)
    )

    fold_balance$n_outcome_1 <- vapply(
      seq_len(n_folds),
      function(fold_index) {
        sum(
          fold_assignments$fold == fold_index &
            fold_assignments$observed == 1
        )
      },
      integer(1)
    )

    if (identical(cv_method, "k-fold")) {
      if (max(fold_balance$n_outcome_0) -
          min(fold_balance$n_outcome_0) > 1L ||
          max(fold_balance$n_outcome_1) -
          min(fold_balance$n_outcome_1) > 1L ||
          any(fold_balance$n_outcome_0 < 1L) ||
          any(fold_balance$n_outcome_1 < 1L)) {
        stop(
          "Internal cross-validation error: binomial k-fold stratification is ",
          "not balanced as required.",
          call. = FALSE
        )
      }
    }
  }

  if (verbose) {
    if (identical(cv_method, "LOOCV")) {
      message(
        "Cross-validation: LOOCV with ", n_groups,
        " complete group(s); one held-out group per fold."
      )
    } else {
      message(
        "Cross-validation: grouped ", k, "-fold with ", n_groups,
        " group(s); fold sizes = ",
        paste(fold_balance$n_groups, collapse = ", "),
        if (cv_stratified) {
          paste0(
            "; binomial stratification 0/1 = ",
            paste(
              paste0(
                fold_balance$n_outcome_0,
                "/",
                fold_balance$n_outcome_1
              ),
              collapse = ", "
            )
          )
        } else {
          ""
        },
        "."
      )
    }
  }

  # Map every long-format row to the single fold assigned to its complete group.
  row_group_match <- match(
    as.character(data_long$epi_id),
    fold_assignments$group
  )

  if (anyNA(row_group_match)) {
    stop(
      "Internal cross-validation error: could not align long-format rows with ",
      "fold assignments.",
      call. = FALSE
    )
  }

  row_fold <- fold_assignments$fold[row_group_match]

  fold_rows <- lapply(
    seq_len(n_folds),
    function(fold_index) which(row_fold == fold_index)
  )

  fold_groups <- lapply(
    seq_len(n_folds),
    function(fold_index) {
      fold_assignments$group[
        fold_assignments$fold == fold_index
      ]
    }
  )

  test_data_list <- lapply(
    fold_rows,
    function(idx) data_long[idx, , drop = FALSE]
  )

  # --------------------------------------------------------------------------
  # Warning classification
  # --------------------------------------------------------------------------

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

  # --------------------------------------------------------------------------
  # Candidate grid
  # --------------------------------------------------------------------------

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

  # --------------------------------------------------------------------------
  # Candidate evaluator
  # --------------------------------------------------------------------------

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

    if (verbose) {
      message(
        "[", candidate_index, "/", total_candidates,
        "] df_var = ", exposure_df,
        ", df_lag = ", lag_df,
        ", vars = ", variables_label
      )
    }

    # Prediction bookkeeping is indexed by independent group, not fold.
    # This distinction is essential for k-fold CV because one fold contains
    # multiple held-out groups.
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

    # Warning bookkeeping remains fold-level because fitting and prediction are
    # executed once per training/test fold.
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

    # One list slot per fold. Under k-fold each slot may contain multiple
    # group-level prediction/failure rows.
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

    cv_stage <- if (identical(cv_method, "LOOCV")) {
      "LOOCV"
    } else {
      "k-fold"
    }

    # IMPORTANT: this loop remains sequential. One future owns one candidate and
    # executes the complete CV workflow for that candidate.
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

            prepared_test <- prepare_response(
              data = test_design,
              y_var = "y",
              family = family_name
            )

            if (!is.data.frame(prepared_test) ||
                nrow(prepared_test) != length(test_groups) ||
                !"epi_id" %in% names(prepared_test) ||
                !"y_model" %in% names(prepared_test)) {
              stop(
                "Internal cross-validation error: the held-out design did not ",
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
                "Internal cross-validation error: held-out group identifiers ",
                "could not be aligned with the fold assignment."
              )
            }

            # CV ranking is based on the harmonized central expected response,
            # population/fixed component only. Coefficient/posterior uncertainty
            # and group-specific random effects are excluded.
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
                "Internal cross-validation error: failed to align predicted ",
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
          stage = cv_stage,
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
            "Internal cross-validation error: duplicate or unknown held-out group.",
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

        # For binomial outcomes, predicted remains the out-of-fold probability.
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
      # valid out-of-fold predictions.
      fold_success[fold_index] <-
        length(valid_prediction) == length(test_groups) &&
        all(valid_prediction)
    }

    n_predictions <- n_groups
    n_success <- sum(success)
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

    all_cv_warning_messages <- unlist(
      candidate_warning_messages,
      recursive = FALSE,
      use.names = FALSE
    )

    all_cv_warning_messages <- as.character(
      all_cv_warning_messages
    )

    unique_cv_warning_messages <- unique(
      all_cv_warning_messages
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
          observed = observed[success],
          predicted = predicted[success],
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
        length(unique_cv_warning_messages),
      n_warning_classes =
        length(
          unique(
            vapply(
              unique_cv_warning_messages,
              classify_warning,
              character(1)
            )
          )
        ),
      warning_classes = if (
        length(unique_cv_warning_messages)
      ) {
        paste(
          unique(
            vapply(
              unique_cv_warning_messages,
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
        length(unique_cv_warning_messages)
      ) {
        paste(
          unique_cv_warning_messages,
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
        " future workers; cross-validation folds remain sequential within each ",
        "candidate."
      )
    }

    candidate_outputs <- future.apply::future_lapply(
      X = seq_len(total_candidates),
      FUN = evaluate_candidate,
      future.seed = TRUE,
      # One future per candidate: folds remain sequential inside each future.
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

  # Warnings have already been muffled inside each fold. Emit at most one
  # concise summary per affected candidate, in deterministic model_id order.
  if (verbose && nrow(warning_summary_data)) {
    affected <- warning_summary_data$n_warning_folds > 0L |
      warning_summary_data$metric_warning_events > 0L |
      warning_summary_data$full_refit_warning_events > 0L

    for (i in which(affected)) {
      ws <- warning_summary_data[i, , drop = FALSE]
      cv_warning_text <- if (ws$n_warning_folds > 0L) {
        paste0(
          ws$n_warning_folds, "/", ws$n_folds,
          " CV folds (", ws$n_warning_events, " warning event(s))"
        )
      } else {
        "0 CV folds"
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
          "; warning classes (affected folds): ",
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
        ": ", cv_warning_text, type_text, metric_text, full_fit_text,
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
      "group-level out-of-fold predictions.",
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
    cv_method = cv_method,
    n_groups = n_groups,
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
  attr(results_data, "cv_method") <- cv_method
  attr(results_data, "cv_scheme") <- cv_scheme
  attr(results_data, "k") <- effective_k
  attr(results_data, "cv_seed") <- if (
    identical(cv_method, "k-fold") && !is.null(seed)
  ) {
    seed
  } else {
    NA_integer_
  }
  attr(results_data, "cv_stratified") <- cv_stratified
  attr(results_data, "fold_assignments") <- fold_assignments
  attr(results_data, "fold_balance") <- fold_balance
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
      " s). CV = ", cv_method,
      if (identical(cv_method, "k-fold")) paste0(" (k=", effective_k, ")") else "",
      "; ranked by ", rank_metric, " (", rank_direction, ")."
    )
  }

  results_data
}
