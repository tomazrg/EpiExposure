# Internal performance-metric helpers for EpiExposure
#
# These helpers are intentionally not exported. Both `find_bestfit()` and
# `ensemble_bestfit()` should call `.compute_performance_metrics()` so that
# model-level and ensemble-level performance are evaluated identically.
#
# EpiExposure v1 performance contract:
#
# - supported families:
#     beta, binomial, poisson, gamma, gaussian, negative_binomial (NB2 only);
# - ordinal outcomes are not supported in EpiExposure v1;
# - NB1 is not supported in EpiExposure v1;
# - non-binary metrics:
#     CCC, Cb, rho, RMSE, MAE;
# - binomial probability metrics:
#     ROC_AUC, Brier, LogLoss;
# - binomial classification metrics:
#     Accuracy, Balanced_Accuracy, Sensitivity, Specificity, F1, MCC, Precision;
# - thresholding is used only for classification metrics;
# - predicted binomial probabilities are retained as probabilities;
# - performance is evaluated on the exact finite support supplied by the
#   caller; non-finite pairs are never discarded silently.
#
# This file does not consume raw exposure histories and therefore has no
# `max_lag + 1` validation responsibility.


# =============================================================================
# FAMILY NORMALIZATION
# =============================================================================

#' Resolve an EpiExposure family to its canonical name
#'
#' Internal helper that mirrors the family contract used by `fit_epidlnm()`,
#' `find_bestfit()`, and `ensemble_bestfit()`.
#'
#' @param family A supported family name or an engine-compatible family object
#'   containing a `family` field.
#'
#' @return One canonical EpiExposure v1 family name:
#'   `"beta"`, `"binomial"`, `"poisson"`, `"gamma"`, `"gaussian"`, or
#'   `"negative_binomial"`.
#'
#' @details
#' EpiExposure v1 supports negative binomial outcomes under the NB2
#' parameterization only. Explicit NB1 family labels are rejected.
#'
#' Ordinal outcomes are not supported in EpiExposure v1 and therefore produce
#' the package-wide explicit error immediately.
#'
#' @keywords internal
#' @noRd
.resolve_family_name <- function(family) {

  family_raw <- NULL

  if (is.character(family) &&
      length(family) == 1L &&
      !is.na(family) &&
      nzchar(trimws(family))) {

    family_raw <- family

  } else if (is.list(family) &&
             !is.null(family$family) &&
             length(family$family) >= 1L &&
             !is.na(family$family[[1L]])) {

    family_raw <- as.character(
      family$family[[1L]]
    )

  } else if (inherits(family, "family") &&
             !is.null(family$family) &&
             length(family$family) >= 1L &&
             !is.na(family$family[[1L]])) {

    family_raw <- as.character(
      family$family[[1L]]
    )

  } else {
    stop(
      "Unsupported `family` specification. Provide a supported EpiExposure ",
      "v1 family name or an engine-compatible family object containing a ",
      "`family` field.",
      call. = FALSE
    )
  }

  if (!length(family_raw) ||
      is.na(family_raw) ||
      !nzchar(trimws(family_raw))) {
    stop(
      "The supplied `family` specification does not contain a usable family ",
      "name.",
      call. = FALSE
    )
  }

  family_raw <- trimws(
    as.character(
      family_raw[[1L]]
    )
  )

  normalized <- tolower(
    family_raw
  )

  compact <- gsub(
    "[^a-z0-9]+",
    "",
    normalized
  )

  # ---------------------------------------------------------------------------
  # Explicitly unsupported families in EpiExposure v1
  # ---------------------------------------------------------------------------

  if (compact %in% c(
    "ordinal",
    "cumulative",
    "cumulativeordinal",
    "ordinalregression"
  ) ||
  grepl(
    "^cumulative",
    compact
  )) {
    stop(
      "Ordinal outcomes are not supported in EpiExposure v1.",
      call. = FALSE
    )
  }

  if (compact %in% c(
    "nbinom1",
    "negativebinomial1",
    "negbinomial1",
    "nb1"
  )) {
    stop(
      "NB1 is not supported in EpiExposure v1. Use the NB2 ",
      "`negative_binomial` parameterization.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Supported canonical families
  # ---------------------------------------------------------------------------

  if (compact %in% c(
    "beta",
    "betafamily",
    "betaproportion",
    "betar",
    "betaregression"
  )) {
    return("beta")
  }

  if (compact %in% c(
    "binomial",
    "bernoulli"
  )) {
    return("binomial")
  }

  if (identical(
    compact,
    "poisson"
  )) {
    return("poisson")
  }

  if (identical(
    compact,
    "gamma"
  )) {
    return("gamma")
  }

  if (compact %in% c(
    "gaussian",
    "normal"
  )) {
    return("gaussian")
  }

  # NB2 is the only negative-binomial contract in EpiExposure v1.
  #
  # Engine objects can expose labels such as:
  #   "nbinom2"
  #   "negative_binomial"
  #   "Negative Binomial(2.3)"
  #   "negbinomial"
  #
  # A number inside "Negative Binomial(...)" is the dispersion parameter,
  # not an NB1/NB2 selector, so these engine labels remain compatible with
  # the canonical NB2 contract.
  if (compact %in% c(
    "nbinom",
    "nbinom2",
    "negativebinomial",
    "negativebinomial2",
    "negbinomial",
    "negbinomial2",
    "nb2"
  ) ||
  grepl(
    "^negativebinomial[0-9]",
    compact
  )) {
    return("negative_binomial")
  }

  stop(
    "Unsupported family: '",
    family_raw,
    "'. Supported EpiExposure v1 canonical families are: beta, binomial, ",
    "poisson, gamma, gaussian, and negative_binomial (NB2 only).",
    call. = FALSE
  )
}


#' Resolve the performance-outcome type from a family
#'
#' @param family A supported EpiExposure v1 family name or family object.
#'
#' @return `"binary"` for binomial outcomes and `"non_binary"` otherwise.
#'
#' @keywords internal
#' @noRd
.resolve_outcome_type <- function(family) {

  family_name <- .resolve_family_name(
    family
  )

  if (identical(
    family_name,
    "binomial"
  )) {
    return("binary")
  }

  "non_binary"
}


# =============================================================================
# METRIC AVAILABILITY AND DIRECTION
# =============================================================================

#' List performance metrics available for a family
#'
#' @param family A supported EpiExposure v1 family name or family object.
#'
#' @return Character vector of canonical performance-metric names.
#'
#' @keywords internal
#' @noRd
.available_metrics <- function(family) {

  outcome_type <- .resolve_outcome_type(
    family
  )

  if (identical(
    outcome_type,
    "non_binary"
  )) {
    return(c(
      "CCC",
      "Cb",
      "rho",
      "RMSE",
      "MAE"
    ))
  }

  c(
    "ROC_AUC",
    "Brier",
    "LogLoss",
    "Accuracy",
    "Balanced_Accuracy",
    "Sensitivity",
    "Specificity",
    "F1",
    "MCC",
    "Precision"
  )
}


#' Determine whether a performance metric is maximized or minimized
#'
#' @param metric Character scalar naming one supported performance metric.
#'
#' @return `"maximize"` or `"minimize"`.
#'
#' @keywords internal
#' @noRd
.metric_direction <- function(metric) {

  if (!is.character(metric) ||
      length(metric) != 1L ||
      is.na(metric) ||
      !nzchar(metric)) {
    stop(
      "`metric` must be one non-empty metric name.",
      call. = FALSE
    )
  }

  maximize <- c(
    "CCC",
    "Cb",
    "rho",
    "ROC_AUC",
    "Accuracy",
    "Balanced_Accuracy",
    "Sensitivity",
    "Specificity",
    "F1",
    "MCC",
    "Precision"
  )

  minimize <- c(
    "RMSE",
    "MAE",
    "Brier",
    "LogLoss"
  )

  if (metric %in%
      maximize) {
    return("maximize")
  }

  if (metric %in%
      minimize) {
    return("minimize")
  }

  stop(
    "Unknown performance metric: '",
    metric,
    "'. Supported metrics are: ",
    paste(
      c(
        maximize,
        minimize
      ),
      collapse = ", "
    ),
    ".",
    call. = FALSE
  )
}


# =============================================================================
# COMMON VECTOR VALIDATION
# =============================================================================

#' Validate observed/predicted vectors for performance evaluation
#'
#' @param observed Numeric vector of observed outcomes.
#' @param predicted Numeric vector of predicted outcomes.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @details
#' Performance is evaluated on the exact supplied support. Non-finite values
#' are not removed pairwise because doing so silently could let a candidate or
#' ensemble be evaluated on fewer observations than the caller supplied.
#'
#' `find_bestfit()` already records failed/non-finite held-out predictions
#' before invoking the central metric helper, and `ensemble_bestfit()` aligns
#' and validates finite OOF support before metric calculation.
#'
#' @keywords internal
#' @noRd
.validate_metric_vectors <- function(
    observed,
    predicted
) {

  if (!is.numeric(observed) ||
      !is.numeric(predicted)) {
    stop(
      "`observed` and `predicted` must be numeric vectors.",
      call. = FALSE
    )
  }

  if (length(observed) !=
      length(predicted)) {
    stop(
      "`observed` and `predicted` must have the same length.",
      call. = FALSE
    )
  }

  if (!length(observed)) {
    stop(
      "`observed` and `predicted` cannot be empty.",
      call. = FALSE
    )
  }

  if (anyNA(observed) ||
      anyNA(predicted) ||
      any(!is.finite(observed)) ||
      any(!is.finite(predicted))) {
    stop(
      "`observed` and `predicted` must contain only finite, non-missing ",
      "values. Performance metrics are evaluated on the exact supplied ",
      "support; non-finite pairs are not removed silently.",
      call. = FALSE
    )
  }

  invisible(
    TRUE
  )
}


# =============================================================================
# NON-BINARY METRICS
# =============================================================================

#' Compute Lin's concordance and numeric prediction-error metrics
#'
#' @param observed Finite numeric vector of observed values.
#' @param predicted Finite numeric vector of predicted values.
#'
#' @return A one-row data frame containing `CCC`, `Cb`, `rho`, `RMSE`, and
#'   `MAE`.
#'
#' @details
#' Lin's concordance correlation coefficient is computed using the standard
#' Lin concordance definition with a plug-in empirical-moment formulation:
#'
#' \deqn{
#'   CCC =
#'   \frac{2s_{xy}}
#'   {s_x^2 + s_y^2 + (\bar{x} - \bar{y})^2}.
#' }
#'
#' The empirical variances and covariance are calculated with divisor `n`.
#' This is the same Lin CCC definition; the moment convention is stated
#' explicitly here so that the finite-sample implementation is transparent and
#' reproducible.
#'
#' With this convention, the standard decomposition
#'
#' \deqn{
#'   CCC = \rho C_b
#' }
#'
#' is evaluated from the same empirical moments, with
#'
#' \deqn{
#'   C_b =
#'   \frac{2s_xs_y}
#'   {s_x^2 + s_y^2 + (\bar{x} - \bar{y})^2}.
#' }
#'
#' `Cb` is calculated directly rather than as `CCC / rho`, which avoids an
#' unstable ratio when `rho` is close to zero.
#'
#' If one series is constant, Pearson `rho` is undefined and returned as `NA`.
#' `CCC` and `Cb` are zero when the concordance denominator is positive. If
#' both vectors are the same constant, the concordance denominator is zero and
#' `CCC`, `Cb`, and `rho` are returned as `NA`; RMSE and MAE remain zero.
#'
#' With fewer than two observations, the three concordance/correlation
#' statistics are `NA`, but RMSE and MAE remain available.
#'
#' @keywords internal
#' @noRd
.ccc_lins <- function(
    observed,
    predicted
) {

  .validate_metric_vectors(
    observed,
    predicted
  )

  rmse <- sqrt(
    mean(
      (observed - predicted)^2
    )
  )

  mae <- mean(
    abs(
      observed - predicted
    )
  )

  if (length(observed) <
      2L) {
    return(
      data.frame(
        CCC = NA_real_,
        Cb = NA_real_,
        rho = NA_real_,
        RMSE = as.numeric(rmse),
        MAE = as.numeric(mae),
        stringsAsFactors = FALSE
      )
    )
  }

  mean_obs <- mean(
    observed
  )

  mean_pred <- mean(
    predicted
  )

  centered_obs <- observed -
    mean_obs

  centered_pred <- predicted -
    mean_pred

  # Empirical population moments (divisor n).
  var_obs <- mean(
    centered_obs^2
  )

  var_pred <- mean(
    centered_pred^2
  )

  covariance <- mean(
    centered_obs *
      centered_pred
  )

  mean_difference_sq <- (
    mean_obs -
      mean_pred
  )^2

  denominator <- var_obs +
    var_pred +
    mean_difference_sq

  ccc <- if (
    is.finite(denominator) &&
    denominator > 0
  ) {
    2 *
      covariance /
      denominator
  } else {
    NA_real_
  }

  cb <- if (
    is.finite(denominator) &&
    denominator > 0 &&
    is.finite(var_obs) &&
    is.finite(var_pred) &&
    var_obs >= 0 &&
    var_pred >= 0
  ) {
    2 *
      sqrt(
        var_obs *
          var_pred
      ) /
      denominator
  } else {
    NA_real_
  }

  rho <- if (
    is.finite(var_obs) &&
    is.finite(var_pred) &&
    var_obs > 0 &&
    var_pred > 0
  ) {
    covariance /
      sqrt(
        var_obs *
          var_pred
      )
  } else {
    NA_real_
  }

  data.frame(
    CCC = as.numeric(ccc),
    Cb = as.numeric(cb),
    rho = as.numeric(rho),
    RMSE = as.numeric(rmse),
    MAE = as.numeric(mae),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# BINARY PROBABILITY METRICS
# =============================================================================

#' Compute ROC AUC from binary observations and predicted probabilities
#'
#' Uses the Mann-Whitney/rank-sum definition of ROC AUC and handles tied
#' probabilities with average ranks.
#'
#' @param observed Finite numeric vector coded as 0/1.
#' @param predicted Finite numeric vector of predicted probabilities in
#'   `[0, 1]`.
#'
#' @return Numeric scalar. Returns `NA_real_` when only one observed class is
#'   present.
#'
#' @keywords internal
#' @noRd
.roc_auc <- function(
    observed,
    predicted
) {

  .validate_metric_vectors(
    observed,
    predicted
  )

  if (!all(
    observed %in%
    c(
      0,
      1
    )
  )) {
    stop(
      "ROC AUC requires binary observed outcomes coded as 0/1.",
      call. = FALSE
    )
  }

  if (any(
    predicted < 0 |
    predicted > 1
  )) {
    stop(
      "ROC AUC requires predicted probabilities in [0, 1].",
      call. = FALSE
    )
  }

  classes <- sort(
    unique(
      observed
    )
  )

  if (length(classes) <
      2L) {
    return(
      NA_real_
    )
  }

  n_pos <- sum(
    observed ==
      1
  )

  n_neg <- sum(
    observed ==
      0
  )

  ranks <- rank(
    predicted,
    ties.method = "average"
  )

  auc <- (
    sum(
      ranks[
        observed ==
          1
      ]
    ) -
      n_pos *
      (
        n_pos +
          1
      ) /
      2
  ) /
    (
      n_pos *
        n_neg
    )

  as.numeric(
    auc
  )
}


#' Return an all-NA binary probability metric row
#'
#' @return One-row data frame containing `ROC_AUC`, `Brier`, and `LogLoss`.
#'
#' @keywords internal
#' @noRd
.empty_binary_probability_metrics <- function() {

  data.frame(
    ROC_AUC = NA_real_,
    Brier = NA_real_,
    LogLoss = NA_real_,
    stringsAsFactors = FALSE
  )
}


#' Compute probability-based metrics for a binary outcome
#'
#' @param observed Finite numeric vector coded as 0/1.
#' @param predicted Finite numeric vector of predicted probabilities in
#'   `[0, 1]`.
#'
#' @return A one-row data frame containing `ROC_AUC`, `Brier`, and `LogLoss`.
#'
#' @details
#' Under the EpiExposure v1 edge-case contract, if the observed vector contains
#' only one class, all binary metrics are returned as `NA`. The central helper
#' emits the warning.
#'
#' Log Loss is evaluated without epsilon clipping. Exact impossible
#' probabilities therefore produce `Inf`, which is the mathematical Log Loss
#' rather than an arbitrary finite value determined by machine precision.
#' Correct exact boundary probabilities contribute zero loss.
#'
#' @keywords internal
#' @noRd
.binary_probability_metrics <- function(
    observed,
    predicted
) {

  .validate_metric_vectors(
    observed,
    predicted
  )

  if (!all(
    observed %in%
    c(
      0,
      1
    )
  )) {
    stop(
      "Binary performance metrics require observed outcomes coded as 0/1.",
      call. = FALSE
    )
  }

  if (any(
    predicted < 0 |
    predicted > 1
  )) {
    stop(
      "Binary probability metrics require predicted probabilities between ",
      "0 and 1.",
      call. = FALSE
    )
  }

  if (length(
    unique(
      observed
    )
  ) <
  2L) {
    return(
      .empty_binary_probability_metrics()
    )
  }

  auc <- .roc_auc(
    observed,
    predicted
  )

  brier <- mean(
    (
      predicted -
        observed
    )^2
  )

  # Direct Bernoulli log-score contributions avoid 0 * log(0) while keeping
  # the mathematically correct infinite penalty for an exact impossible
  # probability.
  log_loss_contribution <- numeric(
    length(
      observed
    )
  )

  positive <- observed ==
    1

  log_loss_contribution[
    positive
  ] <- -log(
    predicted[
      positive
    ]
  )

  log_loss_contribution[
    !positive
  ] <- -log1p(
    -predicted[
      !positive
    ]
  )

  log_loss <- mean(
    log_loss_contribution
  )

  data.frame(
    ROC_AUC = as.numeric(auc),
    Brier = as.numeric(brier),
    LogLoss = as.numeric(log_loss),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# BINARY CLASSIFICATION METRICS
# =============================================================================

#' Return an all-NA binary classification metric row
#'
#' @return One-row data frame containing the canonical threshold-based binary
#'   metrics.
#'
#' @keywords internal
#' @noRd
.empty_binary_classification_metrics <- function() {

  data.frame(
    Accuracy = NA_real_,
    Balanced_Accuracy = NA_real_,
    Sensitivity = NA_real_,
    Specificity = NA_real_,
    F1 = NA_real_,
    MCC = NA_real_,
    Precision = NA_real_,
    stringsAsFactors = FALSE
  )
}


#' Compute threshold-based classification metrics for a binary outcome
#'
#' @param observed Finite numeric vector coded as 0/1.
#' @param predicted Finite numeric vector of predicted probabilities in
#'   `[0, 1]`.
#' @param threshold Finite probability strictly between 0 and 1. It is used
#'   only to convert predicted probabilities to classes.
#'
#' @return A one-row data frame containing `Accuracy`, `Balanced_Accuracy`,
#'   `Sensitivity`, `Specificity`, `F1`, `MCC`, and `Precision`.
#'
#' @details
#' If the observed vector contains only one class, all binary metrics are
#' undefined under the EpiExposure v1 edge-case contract and this helper
#' returns an all-`NA` row.
#'
#' Any metric whose own denominator is zero is returned as `NA`, not `NaN` or
#' `Inf`. For example, `Precision` is undefined when there are no predicted
#' positives, and MCC is undefined when its denominator is zero.
#'
#' @keywords internal
#' @noRd
.binary_classification_metrics <- function(
    observed,
    predicted,
    threshold = 0.5
) {

  .validate_metric_vectors(
    observed,
    predicted
  )

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

  if (!all(
    observed %in%
    c(
      0,
      1
    )
  )) {
    stop(
      "Binary classification metrics require observed outcomes coded as 0/1.",
      call. = FALSE
    )
  }

  if (any(
    predicted < 0 |
    predicted > 1
  )) {
    stop(
      "Binary classification metrics require predicted probabilities between ",
      "0 and 1.",
      call. = FALSE
    )
  }

  if (length(
    unique(
      observed
    )
  ) <
  2L) {
    return(
      .empty_binary_classification_metrics()
    )
  }

  predicted_class <- as.integer(
    predicted >=
      threshold
  )

  tp <- sum(
    observed ==
      1 &
      predicted_class ==
      1
  )

  tn <- sum(
    observed ==
      0 &
      predicted_class ==
      0
  )

  fp <- sum(
    observed ==
      0 &
      predicted_class ==
      1
  )

  fn <- sum(
    observed ==
      1 &
      predicted_class ==
      0
  )

  safe_ratio <- function(
    numerator,
    denominator
  ) {

    if (!is.finite(denominator) ||
        denominator <= 0) {
      return(
        NA_real_
      )
    }

    as.numeric(
      numerator /
        denominator
    )
  }

  accuracy <- safe_ratio(
    tp +
      tn,
    tp +
      tn +
      fp +
      fn
  )

  sensitivity <- safe_ratio(
    tp,
    tp +
      fn
  )

  specificity <- safe_ratio(
    tn,
    tn +
      fp
  )

  precision <- safe_ratio(
    tp,
    tp +
      fp
  )

  f1 <- safe_ratio(
    2 *
      tp,
    2 *
      tp +
      fp +
      fn
  )

  balanced_accuracy <- if (
    is.finite(sensitivity) &&
    is.finite(specificity)
  ) {
    (
      sensitivity +
        specificity
    ) /
      2
  } else {
    NA_real_
  }

  mcc_denominator <- sqrt(
    (
      tp +
        fp
    ) *
      (
        tp +
          fn
      ) *
      (
        tn +
          fp
      ) *
      (
        tn +
          fn
      )
  )

  mcc <- if (
    is.finite(mcc_denominator) &&
    mcc_denominator > 0
  ) {
    (
      (
        tp *
          tn
      ) -
        (
          fp *
            fn
        )
    ) /
      mcc_denominator
  } else {
    NA_real_
  }

  data.frame(
    Accuracy = as.numeric(accuracy),
    Balanced_Accuracy = as.numeric(
      balanced_accuracy
    ),
    Sensitivity = as.numeric(
      sensitivity
    ),
    Specificity = as.numeric(
      specificity
    ),
    F1 = as.numeric(f1),
    MCC = as.numeric(mcc),
    Precision = as.numeric(
      precision
    ),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# CENTRAL PERFORMANCE ENGINE
# =============================================================================

#' Compute EpiExposure performance metrics
#'
#' Central internal performance-metric engine shared by `find_bestfit()` and
#' `ensemble_bestfit()`.
#'
#' @param observed Numeric vector of observed outcomes.
#' @param predicted Numeric vector of model predictions. For
#'   `family = "binomial"`, these must be predicted probabilities in `[0, 1]`,
#'   not already-thresholded class labels.
#' @param family A supported EpiExposure v1 family name or family object.
#' @param threshold Numeric probability threshold strictly between 0 and 1.
#'   Default is `0.5`. It is used only for binary classification metrics.
#'   ROC AUC, Brier score, and Log Loss always use the original probabilities.
#'
#' @return A one-row data frame.
#'
#'   Non-binary outcomes return `CCC`, `Cb`, `rho`, `RMSE`, and `MAE`.
#'
#'   Binomial outcomes return `ROC_AUC`, `Brier`, `LogLoss`, `Accuracy`,
#'   `Balanced_Accuracy`, `Sensitivity`, `Specificity`, `F1`, `MCC`, and
#'   `Precision`.
#'
#' @details
#' ## Supported families
#'
#' Supported non-binary families are `beta`, `poisson`, `gamma`, `gaussian`,
#' and `negative_binomial`. EpiExposure v1 uses NB2 only. `binomial` is treated
#' as a binary outcome and requires observed values coded exactly as 0/1.
#'
#' Ordinal outcomes are not supported in EpiExposure v1. NB1 is also not
#' supported. Both stop explicitly in the family-normalization layer.
#'
#' ## Exact evaluation support
#'
#' `observed` and `predicted` must have the same positive length and contain
#' only finite, non-missing values. The metric engine does not silently drop
#' incomplete pairs. This keeps performance comparisons tied to the exact
#' support supplied by the caller.
#'
#' `find_bestfit()` separately records failed held-out predictions and passes
#' only its explicitly successful support to this helper.
#' `ensemble_bestfit()` requires aligned finite OOF support before calling the
#' same helper.
#'
#' ## Non-binary metrics
#'
#' Lin's CCC, its bias-correction factor `Cb`, and Pearson correlation `rho`
#' evaluate concordance/correlation, whereas RMSE and MAE quantify prediction
#' error. CCC uses the standard Lin definition with empirical variances and
#' covariance calculated using divisor `n`; this convention is documented
#' explicitly for finite-sample reproducibility. `Cb` is computed directly
#' rather than as `CCC / rho`.
#'
#' ## Binary probability and classification metrics
#'
#' Binomial predictions remain probabilities. ROC AUC, Brier score, and Log
#' Loss are calculated directly from those probabilities. `threshold` is
#' applied only to construct the auxiliary 0/1 predicted class used by
#' threshold-dependent metrics.
#'
#' If supplied binomial probabilities contain only exact 0/1 values, a warning
#' is emitted because probability-resolution information has been lost,
#' although 0 and 1 remain mathematically valid probabilities.
#'
#' Log Loss is not epsilon-clipped. Exact impossible probabilities therefore
#' produce `Inf`, preserving the mathematical scoring rule.
#'
#' ## Binary edge cases
#'
#' If `observed` contains only class 0 or only class 1, EpiExposure v1 returns
#' all binary metrics as `NA` and emits one warning. This prevents model ranking
#' from a validation support that does not contain both outcome classes.
#'
#' If both observed classes are present but the selected threshold yields only
#' one predicted class, metrics that remain mathematically defined are retained;
#' metrics with zero denominators return `NA`. A warning records this condition.
#'
#' @keywords internal
#' @noRd
.compute_performance_metrics <- function(
    observed,
    predicted,
    family,
    threshold = 0.5
) {

  .validate_metric_vectors(
    observed,
    predicted
  )

  family_name <- .resolve_family_name(
    family
  )

  outcome_type <- .resolve_outcome_type(
    family_name
  )

  if (identical(
    outcome_type,
    "non_binary"
  )) {
    return(
      .ccc_lins(
        observed,
        predicted
      )
    )
  }

  # ---------------------------------------------------------------------------
  # Binomial validation
  # ---------------------------------------------------------------------------

  if (!all(
    observed %in%
    c(
      0,
      1
    )
  )) {
    stop(
      "`family = 'binomial'` requires observed outcomes coded exactly as 0/1 ",
      "for performance evaluation.",
      call. = FALSE
    )
  }

  if (any(
    predicted < 0 |
    predicted > 1
  )) {
    stop(
      "`family = 'binomial'` requires predicted probabilities between 0 and ",
      "1 for performance evaluation.",
      call. = FALSE
    )
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

  observed_classes <- sort(
    unique(
      observed
    )
  )

  if (length(
    observed_classes
  ) <
  2L) {

    warning(
      "The observed binary outcome contains only class ",
      observed_classes[[1L]],
      ". EpiExposure v1 requires both observed classes for binary performance ",
      "evaluation; all binary metrics are returned as NA.",
      call. = FALSE
    )

    return(
      cbind(
        .empty_binary_probability_metrics(),
        .empty_binary_classification_metrics()
      )
    )
  }

  if (all(
    predicted %in%
    c(
      0,
      1
    )
  )) {
    warning(
      "The supplied binomial predictions contain only exact 0/1 values. ",
      "They are treated as probabilities, but probability-based metrics no ",
      "longer retain probability-resolution information.",
      call. = FALSE
    )
  }

  predicted_class <- as.integer(
    predicted >=
      threshold
  )

  if (length(
    unique(
      predicted_class
    )
  ) <
  2L) {
    warning(
      "At threshold = ",
      format(
        threshold
      ),
      ", all predicted classes are ",
      unique(
        predicted_class
      )[[1L]],
      ". Threshold-based metrics with zero denominators are returned as NA ",
      "where appropriate.",
      call. = FALSE
    )
  }

  probability_metrics <-
    .binary_probability_metrics(
      observed = observed,
      predicted = predicted
    )

  classification_metrics <-
    .binary_classification_metrics(
      observed = observed,
      predicted = predicted,
      threshold = threshold
    )

  cbind(
    probability_metrics,
    classification_metrics
  )
}
