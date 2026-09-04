# Internal performance-metric helpers for EpiExposure
#
# These helpers are intentionally not exported. Both `find_bestfit()` and
# `ensemble_bestfit()` should call `.compute_performance_metrics()` so that
# model-level and ensemble-level performance are evaluated identically.

#' Resolve an EpiExposure family to its canonical name
#'
#' Internal helper that mirrors the family normalization used by
#' `fit_epidlnm()`.
#'
#' @param family A supported family name or an engine-compatible family object
#'   containing a `family` field.
#'
#' @return One canonical family name.
#' @keywords internal
#' @noRd
.resolve_family_name <- function(family) {
  family_raw <- NULL
  
  if (is.character(family) &&
      length(family) == 1L &&
      !is.na(family) &&
      nzchar(family)) {
    family_raw <- family
  } else if (is.list(family) &&
             !is.null(family$family) &&
             length(family$family) >= 1L) {
    family_raw <- as.character(family$family[[1]])
  } else if (inherits(family, "family") &&
             !is.null(family$family)) {
    family_raw <- as.character(family$family[[1]])
  } else {
    stop(
      "Unsupported `family` specification. Provide a supported family name ",
      "or an engine-compatible family object containing a `family` field.",
      call. = FALSE
    )
  }
  
  normalized <- tolower(trimws(family_raw))
  normalized <- gsub("[[:space:]-]+", "_", normalized)
  normalized <- gsub("[^a-z0-9_]", "", normalized)
  
  if (normalized %in% c("beta", "beta_family", "beta_proportion")) {
    return("beta")
  }
  
  if (normalized %in% c("binomial", "bernoulli")) {
    return("binomial")
  }
  
  if (normalized == "poisson") {
    return("poisson")
  }
  
  if (normalized == "gamma") {
    return("gamma")
  }
  
  if (normalized %in% c("gaussian", "normal")) {
    return("gaussian")
  }
  
  if (
    normalized %in% c(
      "negbin", "nbinom", "nbinom1", "nbinom2",
      "negative_binomial", "negative_binomial_1",
      "negative_binomial_2"
    ) ||
    grepl("negative.*binomial", normalized)
  ) {
    return("negative_binomial")
  }
  
  if (normalized %in% c("ordinal", "cumulative")) {
    return("ordinal")
  }
  
  stop(
    "Unsupported family: '", family_raw, "'. Supported canonical families are: ",
    "beta, binomial, poisson, gamma, gaussian, negative_binomial, and ordinal.",
    call. = FALSE
  )
}


#' Resolve the performance-outcome type from a family
#'
#' @param family A supported family name or family object.
#'
#' @return One of `"binary"`, `"non_binary"`, or `"ordinal"`.
#' @keywords internal
#' @noRd
.resolve_outcome_type <- function(family) {
  family_name <- .resolve_family_name(family)
  
  if (identical(family_name, "binomial")) {
    return("binary")
  }
  
  if (identical(family_name, "ordinal")) {
    return("ordinal")
  }
  
  "non_binary"
}


#' List performance metrics available for a family
#'
#' @param family A supported family name or family object.
#'
#' @return Character vector of metric names. Ordinal performance metrics are
#'   intentionally not defined yet and therefore return `character(0)`.
#' @keywords internal
#' @noRd
.available_metrics <- function(family) {
  outcome_type <- .resolve_outcome_type(family)
  
  if (identical(outcome_type, "non_binary")) {
    return(c(
      "CCC",
      "Cb",
      "rho",
      "RMSE",
      "MAE"
    ))
  }
  
  if (identical(outcome_type, "binary")) {
    return(c(
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
    ))
  }
  
  character(0)
}


#' Determine whether a performance metric is maximized or minimized
#'
#' @param metric Character scalar naming one supported performance metric.
#'
#' @return `"maximize"` or `"minimize"`.
#' @keywords internal
#' @noRd
.metric_direction <- function(metric) {
  if (!is.character(metric) ||
      length(metric) != 1L ||
      is.na(metric) ||
      !nzchar(metric)) {
    stop("`metric` must be one non-empty metric name.", call. = FALSE)
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
  
  if (metric %in% maximize) {
    return("maximize")
  }
  
  if (metric %in% minimize) {
    return("minimize")
  }
  
  stop(
    "Unknown performance metric: '", metric, "'. Supported metrics are: ",
    paste(c(maximize, minimize), collapse = ", "),
    ".",
    call. = FALSE
  )
}


#' Compute Lin's concordance and numeric prediction-error metrics
#'
#' @param observed Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#'
#' @return A one-row data frame containing `CCC`, `Cb`, `rho`, `RMSE`, and
#'   `MAE`.
#' @keywords internal
#' @noRd
.ccc_lins <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  
  if (length(observed) < 2L) {
    return(data.frame(
      CCC = NA_real_,
      Cb = NA_real_,
      rho = NA_real_,
      RMSE = NA_real_,
      MAE = NA_real_
    ))
  }
  
  mean_obs <- mean(observed)
  mean_pred <- mean(predicted)
  var_obs <- stats::var(observed)
  var_pred <- stats::var(predicted)
  covariance <- stats::cov(observed, predicted)
  rho <- suppressWarnings(stats::cor(observed, predicted))
  
  denominator <- var_obs + var_pred + (mean_obs - mean_pred)^2
  
  ccc <- if (is.finite(denominator) &&
             denominator > .Machine$double.eps) {
    (2 * covariance) / denominator
  } else {
    NA_real_
  }
  
  cb <- if (is.finite(rho) &&
            abs(rho) > .Machine$double.eps &&
            is.finite(ccc)) {
    ccc / rho
  } else {
    NA_real_
  }
  
  data.frame(
    CCC = as.numeric(ccc),
    Cb = as.numeric(cb),
    rho = as.numeric(rho),
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted))
  )
}


#' Compute ROC AUC from binary observations and predicted probabilities
#'
#' Uses the rank-sum definition of the area under the ROC curve and handles
#' tied probabilities through average ranks.
#'
#' @param observed Numeric vector coded as 0/1.
#' @param predicted Numeric vector of predicted probabilities.
#'
#' @return A numeric scalar. Returns `NA_real_` when only one observed class is
#'   present.
#' @keywords internal
#' @noRd
.roc_auc <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  
  if (!length(observed)) {
    return(NA_real_)
  }
  
  classes <- sort(unique(observed))
  
  if (!all(classes %in% c(0, 1))) {
    stop(
      "ROC AUC requires binary observed outcomes coded as 0/1.",
      call. = FALSE
    )
  }
  
  if (length(classes) < 2L) {
    return(NA_real_)
  }
  
  n_pos <- sum(observed == 1)
  n_neg <- sum(observed == 0)
  
  ranks <- rank(predicted, ties.method = "average")
  
  auc <- (
    sum(ranks[observed == 1]) -
      n_pos * (n_pos + 1) / 2
  ) / (n_pos * n_neg)
  
  as.numeric(auc)
}


#' Compute probability-based metrics for a binary outcome
#'
#' @param observed Numeric vector coded as 0/1.
#' @param predicted Numeric vector of predicted probabilities in `[0, 1]`.
#'
#' @return A one-row data frame containing `ROC_AUC`, `Brier`, and `LogLoss`.
#' @keywords internal
#' @noRd
.binary_probability_metrics <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  
  if (!length(observed)) {
    return(data.frame(
      ROC_AUC = NA_real_,
      Brier = NA_real_,
      LogLoss = NA_real_
    ))
  }
  
  if (!all(unique(observed) %in% c(0, 1))) {
    stop(
      "Binary performance metrics require observed outcomes coded as 0/1.",
      call. = FALSE
    )
  }
  
  if (any(predicted < 0 | predicted > 1)) {
    stop(
      "Binary probability metrics require predicted probabilities between 0 and 1.",
      call. = FALSE
    )
  }
  
  auc <- .roc_auc(observed, predicted)
  brier <- mean((predicted - observed)^2)
  
  eps <- .Machine$double.eps
  p_log <- pmin(pmax(predicted, eps), 1 - eps)
  
  log_loss <- -mean(
    observed * log(p_log) +
      (1 - observed) * log(1 - p_log)
  )
  
  data.frame(
    ROC_AUC = as.numeric(auc),
    Brier = as.numeric(brier),
    LogLoss = as.numeric(log_loss)
  )
}


#' Compute threshold-based classification metrics for a binary outcome
#'
#' @param observed Numeric vector coded as 0/1.
#' @param predicted Numeric vector of predicted probabilities in `[0, 1]`.
#' @param threshold Numeric probability threshold strictly between 0 and 1 used
#'   to convert predicted probabilities to classes 0/1. The threshold is meant
#'   for predicted probabilities, not predictions that have already been
#'   reduced to class labels 0/1.
#'
#' @return A one-row data frame containing `Accuracy`, `Balanced_Accuracy`,
#'   `Sensitivity`, `Specificity`, `F1`, `MCC`, and `Precision`.
#' @keywords internal
#' @noRd
.binary_classification_metrics <- function(
    observed,
    predicted,
    threshold = 0.5
) {
  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      !is.finite(threshold) ||
      threshold <= 0 ||
      threshold >= 1) {
    stop(
      "`threshold` must be one finite numeric probability strictly between 0 and 1.",
      call. = FALSE
    )
  }
  
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  
  if (!length(observed)) {
    return(data.frame(
      Accuracy = NA_real_,
      Balanced_Accuracy = NA_real_,
      Sensitivity = NA_real_,
      Specificity = NA_real_,
      F1 = NA_real_,
      MCC = NA_real_,
      Precision = NA_real_
    ))
  }
  
  if (!all(unique(observed) %in% c(0, 1))) {
    stop(
      "Binary classification metrics require observed outcomes coded as 0/1.",
      call. = FALSE
    )
  }
  
  if (any(predicted < 0 | predicted > 1)) {
    stop(
      "Binary classification metrics require predicted probabilities between 0 and 1.",
      call. = FALSE
    )
  }
  
  predicted_class <- as.integer(predicted >= threshold)
  
  tp <- sum(observed == 1 & predicted_class == 1)
  tn <- sum(observed == 0 & predicted_class == 0)
  fp <- sum(observed == 0 & predicted_class == 1)
  fn <- sum(observed == 1 & predicted_class == 0)
  
  safe_ratio <- function(numerator, denominator) {
    if (!is.finite(denominator) ||
        denominator <= 0) {
      return(NA_real_)
    }
    as.numeric(numerator / denominator)
  }
  
  accuracy <- safe_ratio(tp + tn, tp + tn + fp + fn)
  sensitivity <- safe_ratio(tp, tp + fn)
  specificity <- safe_ratio(tn, tn + fp)
  precision <- safe_ratio(tp, tp + fp)
  f1 <- safe_ratio(2 * tp, 2 * tp + fp + fn)
  
  balanced_accuracy <- if (is.finite(sensitivity) &&
                           is.finite(specificity)) {
    (sensitivity + specificity) / 2
  } else {
    NA_real_
  }
  
  mcc_denominator <- sqrt(
    (tp + fp) *
      (tp + fn) *
      (tn + fp) *
      (tn + fn)
  )
  
  mcc <- if (is.finite(mcc_denominator) &&
             mcc_denominator > 0) {
    ((tp * tn) - (fp * fn)) / mcc_denominator
  } else {
    NA_real_
  }
  
  data.frame(
    Accuracy = as.numeric(accuracy),
    Balanced_Accuracy = as.numeric(balanced_accuracy),
    Sensitivity = as.numeric(sensitivity),
    Specificity = as.numeric(specificity),
    F1 = as.numeric(f1),
    MCC = as.numeric(mcc),
    Precision = as.numeric(precision)
  )
}


#' Compute EpiExposure performance metrics
#'
#' Central internal performance-metric engine shared by `find_bestfit()` and
#' `ensemble_bestfit()`.
#'
#' @param observed Numeric vector of observed outcomes.
#' @param predicted Numeric vector of model predictions. For
#'   `family = "binomial"`, these must be predicted probabilities in `[0, 1]`,
#'   not already-thresholded class labels.
#' @param family A supported EpiExposure family name or family object.
#' @param threshold Numeric probability threshold strictly between 0 and 1.
#'   It is used only for binary classification metrics to convert predicted
#'   probabilities into classes 0/1. It does not affect ROC AUC, Brier score,
#'   or Log Loss. The threshold is intended for probabilities and should not be
#'   interpreted as a rule to be applied to predictions that are already class
#'   labels 0/1.
#'
#' @return A one-row data frame. Non-binary outcomes return `CCC`, `Cb`, `rho`,
#'   `RMSE`, and `MAE`. Binary outcomes return `ROC_AUC`, `Brier`, `LogLoss`,
#'   `Accuracy`, `Balanced_Accuracy`, `Sensitivity`, `Specificity`, `F1`, `MCC`,
#'   and `Precision`.
#'
#' @details
#' Supported non-binary families are `beta`, `poisson`, `gamma`, `gaussian`,
#' and `negative_binomial`. The `binomial` family is treated as a binary
#' outcome and requires observed values coded as 0/1.
#'
#' The `ordinal` family is recognized by the family-normalization layer but
#' ordinal-specific performance metrics have not yet been defined. Calling this
#' function with `family = "ordinal"` therefore produces an explicit error
#' rather than silently applying numeric or binary metrics.
#'
#' If a binary observed outcome contains only one class, metrics requiring both
#' classes are undefined. The function returns `NA` where appropriate and emits
#' one warning explaining why. Likewise, metrics whose denominator is zero
#' return `NA` rather than `NaN` or `Inf`.
#'
#' @keywords internal
#' @noRd
.compute_performance_metrics <- function(
    observed,
    predicted,
    family,
    threshold = 0.5
) {
  if (!is.numeric(observed) ||
      !is.numeric(predicted)) {
    stop(
      "`observed` and `predicted` must be numeric vectors.",
      call. = FALSE
    )
  }
  
  if (length(observed) != length(predicted)) {
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
  
  ok <- is.finite(observed) & is.finite(predicted)
  observed <- observed[ok]
  predicted <- predicted[ok]
  
  if (!length(observed)) {
    stop(
      "No finite observed/predicted pairs are available for performance evaluation.",
      call. = FALSE
    )
  }
  
  family_name <- .resolve_family_name(family)
  outcome_type <- .resolve_outcome_type(family_name)
  
  if (identical(outcome_type, "ordinal")) {
    stop(
      "Performance metrics for `family = 'ordinal'` have not yet been defined. ",
      "Ordinal outcomes are recognized by EpiExposure but are not silently ",
      "evaluated with continuous or binary metrics.",
      call. = FALSE
    )
  }
  
  if (identical(outcome_type, "non_binary")) {
    return(.ccc_lins(observed, predicted))
  }
  
  if (!all(unique(observed) %in% c(0, 1))) {
    stop(
      "`family = 'binomial'` requires observed outcomes coded as 0/1 ",
      "for performance evaluation.",
      call. = FALSE
    )
  }
  
  if (any(predicted < 0 | predicted > 1)) {
    stop(
      "`family = 'binomial'` requires predicted probabilities between 0 and 1 ",
      "for performance evaluation.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      !is.finite(threshold) ||
      threshold <= 0 ||
      threshold >= 1) {
    stop(
      "`threshold` must be one finite numeric probability strictly between 0 and 1.",
      call. = FALSE
    )
  }
  
  observed_classes <- sort(unique(observed))
  
  if (length(observed_classes) < 2L) {
    warning(
      "The observed binary outcome contains only class ",
      observed_classes[[1]],
      ". ROC AUC and classification metrics requiring the missing class ",
      "are undefined and are returned as NA where appropriate.",
      call. = FALSE
    )
  }
  
  if (all(predicted %in% c(0, 1))) {
    warning(
      "`threshold` is intended to convert predicted probabilities into 0/1 ",
      "classes, but the supplied predictions contain only 0/1 values. ",
      "Classification metrics can still be computed, but probability-based ",
      "metrics no longer retain probability-resolution information.",
      call. = FALSE
    )
  }
  
  probability_metrics <- .binary_probability_metrics(
    observed = observed,
    predicted = predicted
  )
  
  classification_metrics <- .binary_classification_metrics(
    observed = observed,
    predicted = predicted,
    threshold = threshold
  )
  
  cbind(
    probability_metrics,
    classification_metrics
  )
}
