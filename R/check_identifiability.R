#' Diagnose identifiability and numerical stability of DLNM cross-basis designs
#'
#' Diagnoses the epidemic-level DLNM design that EpiExposure would construct for
#' model fitting from one or more exposure variables. The function deliberately
#' reuses the canonical EpiExposure design path:
#'
#' ```
#' templates <- define_exposures(...)
#' design    <- build_design(data, templates, include_response = FALSE)
#' ```
#'
#' This avoids maintaining a second, potentially divergent implementation of
#' cross-basis construction inside the diagnostic function. Consequently, the
#' basis definitions, training-data-dependent knots, lag convention, marginal
#' intercept rules, canonical cross-basis column order, and final-row design
#' used here are the same as those used by the current EpiExposure fitting
#' workflow.
#'
#' @param data Non-empty long-format exposure data frame containing `epi_id`,
#'   `time`, and every exposure named in `vars`.
#'
#'   EpiExposure v1 uses a strict exact-history contract. If `max_lag = L`,
#'   every epidemic must contain exactly `L + 1` rows. Histories with fewer or
#'   more rows are rejected; they are never truncated, padded, or silently
#'   realigned. All requested exposure variables occupy those same validated
#'   rows and therefore necessarily have the same temporal support and history
#'   length.
#'
#'   `epi_id` must not contain missing values. `time` and all exposure variables
#'   must contain only finite numeric values. When `max_lag > 0`, time must be
#'   strictly increasing after ordering, equally spaced within each epidemic,
#'   and use the same spacing across epidemics. Missing or non-finite exposure
#'   values are not dropped in EpiExposure v1; they are errors.
#' @param vars Character vector with one or more unique exposure-variable names.
#'   The reserved EpiExposure names `"epi_id"`, `"time"`, and `"y"` cannot be
#'   used as exposure names.
#' @param max_lag Non-negative integer maximum lag. The legacy form `c(0, L)` is
#'   accepted for backward compatibility. EpiExposure v1 requires the fitted
#'   lag range to begin at lag 0 and all fitted exposures to use the same
#'   `max_lag`.
#' @param df_var Positive integer controlling the exposure-response basis when
#'   `fun_var` is `"ns"` or `"bs"`. For `fun_var = "poly"`, it is used as the
#'   polynomial degree. It is ignored by the linear basis.
#' @param df_lag Positive integer controlling the lag-response basis when
#'   `fun_lag` is `"ns"` or `"bs"`. For `fun_lag = "poly"`, it is used as the
#'   polynomial degree. It is ignored by the linear basis.
#' @param fun_var Character exposure-basis function. Supported unpenalized
#'   EpiExposure v1 choices are `"ns"`, `"bs"`, `"poly"`, and `"lin"`.
#' @param fun_lag Character lag-basis function. Supported unpenalized
#'   EpiExposure v1 choices are `"ns"`, `"bs"`, `"poly"`, and `"lin"`.
#'
#'   Penalized dlnm basis functions `"ps"` and `"cr"` are deliberately not
#'   supported in EpiExposure v1. A penalized spline transformation alone is
#'   not a penalized DLNM: the associated penalty matrices and smoothing
#'   parameters must also be propagated to model fitting. The current
#'   engine-agnostic fitting contract therefore stops explicitly instead of
#'   diagnosing a basis that EpiExposure would later fit without its penalty.
#' @param include_intercept Logical scalar. If `TRUE` (default), include a model
#'   intercept when evaluating numerical rank of the combined epidemic-level
#'   design. This matches the standard EpiExposure fitting parameterization.
#'   `FALSE` is retained as an advanced diagnostic option.
#' @param corr_threshold Finite number in `(0, 1]`. Heuristic absolute
#'   correlation threshold used to flag potentially important between-exposure
#'   cross-basis collinearity and to report high same-time raw-exposure
#'   correlation. Default is `0.90`.
#' @param condition_warn Positive finite number. Heuristic scaled condition
#'   number at or above which a numerical warning is reported. Default is `30`.
#' @param condition_severe Positive finite number strictly larger than
#'   `condition_warn`. At or above this value the design is classified as
#'   numerically unstable. Default is `100`.
#' @param vif_threshold Positive finite number. Threshold for supplementary VIF
#'   reporting. Default is `10`. VIFs for individual spline columns can be high
#'   by construction and do not independently determine the global status.
#' @param n_per_parameter_warn Positive finite number. Heuristic minimum ratio
#'   of complete epidemics to columns in the rank design. Default is `10`.
#'   This is a design-complexity diagnostic, not a universal sample-size rule.
#' @param tol Positive finite number used for numerical QR-rank decisions,
#'   near-zero-variance checks, and numerical zero checks in the supplementary
#'   VIF calculation. Default is `1e-7`. Temporal validation itself follows the
#'   canonical tolerances in `define_exposures()` and `build_design()`.
#' @param keep_design Logical scalar. If `TRUE`, include the complete
#'   epidemic-level numerical matrix used for the rank diagnostic in the
#'   returned object. Default is `FALSE`.
#'
#' @return An object of class `epiexposure_identifiability`, a list containing:
#'
#'   - `identifiable`: whether the combined rank design has full numerical rank
#'     at `tol`;
#'   - `numerically_stable`: whether it is full rank, has no near-zero-variance
#'     cross-basis columns, has a scaled condition number below
#'     `condition_severe`, and contains more epidemics than rank-design columns;
#'   - `status`: `"ok"`, `"warning"`, or `"problem"`;
#'   - `overall`: one-row summary of the combined design;
#'   - `by_variable`: diagnostics for each exposure's cross-basis block;
#'   - `pairwise_exposure_correlation`: descriptive row-level Pearson and
#'     Spearman correlations between raw exposures;
#'   - `pairwise_crossbasis_correlation`: maximum and mean absolute pairwise
#'     correlations between columns belonging to different exposure blocks;
#'   - `vif`: supplementary column-wise VIF diagnostics;
#'   - `singular_values`: singular values of the centered/scaled usable
#'     cross-basis predictor matrix;
#'   - `dependent_columns`: pivot-based set of columns not needed to span a
#'     rank-deficient combined design;
#'   - `near_zero_variance_columns`: flagged cross-basis columns;
#'   - `complete_epi_id`: epidemic IDs represented in the design;
#'   - `basis_specification`: effective basis metadata transported from
#'     `define_exposures()`/`build_design()`;
#'   - `diagnostic_flags`: named logical flags used to construct the status;
#'   - `recommendations`: concise interpretation and follow-up suggestions;
#'   - `settings`: normalized diagnostic settings and temporal metadata;
#'   - `design_matrix`: included only when `keep_design = TRUE`.
#'
#' @details
#' ## What is diagnosed
#'
#' EpiExposure fits one epidemic-level row after transforming each complete
#' exposure history through its training cross-basis and retaining the final
#' cross-basis row. `check_identifiability()` diagnoses that same numerical
#' design. It does not diagnose every intermediate row returned by
#' `dlnm::crossbasis()` within a history.
#'
#' The function first calls `define_exposures()` to estimate the training basis
#' definitions jointly from all supplied epidemics, then calls `build_design()`
#' to reconstruct one final cross-basis row per epidemic using those fixed
#' training definitions. This is the authoritative EpiExposure v1 design path.
#'
#' ## Exact-history and common-window contract
#'
#' If the fitted maximum lag is `L`, each epidemic must contain exactly
#'
#' \deqn{
#'   L + 1
#' }
#'
#' chronological exposure observations, corresponding to lag `L` through lag
#' `0`. All fitted exposures use the same `L`. A longer history is not silently
#' reduced to its last `L + 1` rows, and a shorter history is not padded.
#'
#' For `max_lag = 0`, each epidemic contains exactly one observation and the lag
#' basis must be linear. A one-point history has no estimable temporal step, so
#' `time_step` is stored as `NA`.
#'
#' ## Identifiability versus numerical stability
#'
#' `identifiable` is based on the **numerical rank** of the epidemic-level
#' design at the user-selected tolerance; it is not an algebraic symbolic-rank
#' proof. When `include_intercept = TRUE`, the intercept is included in this
#' rank calculation.
#'
#' Numerical stability is assessed separately. The condition number is computed
#' from usable cross-basis predictor columns after centering and scaling each
#' column to unit standard deviation. The intercept is not included in this
#' scaled condition number because centering makes it orthogonal to the scaled
#' predictor columns. If the scaled predictor matrix is rank deficient, or if
#' the number of usable columns cannot be supported by its rows, the condition
#' number is reported as infinite.
#'
#' Condition-number cutoffs are heuristics, not universal inferential laws.
#' Similarly, the epidemics-per-column ratio is a descriptive warning about
#' design complexity and should not be interpreted as a formal sample-size
#' calculation.
#'
#' ## Correlations and VIFs
#'
#' Correlations among columns from the same spline basis are expected and are
#' not, by themselves, evidence that the model is invalid. The combined rank
#' and scaled condition number are the primary numerical diagnostics.
#'
#' `pairwise_crossbasis_correlation` examines **different exposure blocks** and
#' can therefore highlight two exposures that contribute highly redundant
#' transformed temporal information. It reports the maximum pairwise column
#' correlation; this is not a canonical-correlation analysis.
#'
#' Raw-exposure correlations are calculated on the long-format rows and are
#' purely descriptive. Repeated/serial observations are not treated as
#' independent observations for inferential testing, and the function does not
#' report correlation p-values. High raw correlation alone does not change the
#' global `status` when the fitted cross-basis design remains well behaved.
#'
#' VIFs are also supplementary because spline-expanded columns are commonly
#' correlated by construction. A high VIF produces a recommendation but, on its
#' own, does not upgrade an otherwise acceptable design to `"warning"`.
#'
#' ## Penalized DLNMs
#'
#' EpiExposure v1 is intentionally unpenalized. `"ps"` and `"cr"` are rejected
#' here for the same reason they are rejected by `define_exposures()`: the
#' future penalized framework must carry the basis, penalty matrices, smoothing
#' parameters, fitting engine, prediction design, and uncertainty contract as a
#' coherent model. Diagnosing `"ps"` columns while fitting them later without
#' their penalty would be misleading.
#'
#' @export
check_identifiability <- function(
    data,
    vars,
    max_lag,
    df_var = 4,
    df_lag = 4,
    fun_var = "ns",
    fun_lag = "ns",
    include_intercept = TRUE,
    corr_threshold = 0.90,
    condition_warn = 30,
    condition_severe = 100,
    vif_threshold = 10,
    n_per_parameter_warn = 10,
    tol = 1e-7,
    keep_design = FALSE
) {

  call_args <- names(as.list(sys.call())[-1L])

  if ("var" %in% call_args){
    stop(
      "`var` is no longer supported. Use `vars` instead.",
      call. = FALSE
    )
  }

  # LOCAL VALIDATORS

  stopf <- function(...) {
    stop(..., call. = FALSE)
  }

  logical_scalar <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  character_scalar <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }

  positive_numeric_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x <= 0) {
      stopf("`", name, "` must be one positive finite number.")
    }
    as.numeric(x)
  }

  positive_integer_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x < 1 || x != as.integer(x)) {
      stopf("`", name, "` must be one positive integer.")
    }
    as.integer(x)
  }

  normalize_max_lag <- function(x) {
    if (!is.numeric(x) || !length(x) || length(x) > 2L ||
        anyNA(x) || any(!is.finite(x)) || any(x < 0) ||
        any(x != as.integer(x))) {
      stopf(
        "`max_lag` must be one non-negative integer or the legacy lag range ",
        "`c(0, L)`."
      )
    }

    x <- as.integer(x)

    if (!is.null(names(x)) && any(nzchar(names(x)))) {
      stopf(
        "Named/exposure-specific `max_lag` values are not supported. All ",
        "fitted exposure variables must use one common `max_lag`."
      )
    }

    if (length(x) == 1L) {
      return(x)
    }

    if (min(x) != 0L) {
      stopf(
        "EpiExposure v1 requires the fitted lag range to start at lag 0. ",
        "Use `max_lag = L` or `max_lag = c(0, L)`."
      )
    }

    as.integer(max(x))
  }

  # BASIC INPUT VALIDATION

  if (!is.data.frame(data) || !nrow(data)) {
    stopf("`data` must be a non-empty data.frame.")
  }

  missing_required <- setdiff(c("epi_id", "time"), names(data))
  if (length(missing_required)) {
    stopf(
      "`data` is missing required column(s): ",
      paste(missing_required, collapse = ", "),
      "."
    )
  }

  if (!is.character(vars) || !length(vars) || anyNA(vars) ||
      any(!nzchar(vars)) || anyDuplicated(vars)) {
    stopf("`vars` must contain unique non-empty exposure-variable names.")
  }

  if (any(vars %in% c("epi_id", "time", "y"))) {
    stopf(
      "Exposure variables cannot use the reserved EpiExposure names ",
      "'epi_id', 'time', or 'y'."
    )
  }

  missing_variables <- setdiff(vars, names(data))
  if (length(missing_variables)) {
    stopf(
      "Exposure variable(s) not found in `data`: ",
      paste(missing_variables, collapse = ", "),
      "."
    )
  }

  if (anyNA(data$epi_id)) {
    stopf("`data$epi_id` cannot contain missing values.")
  }

  if (!is.numeric(data$time) || anyNA(data$time) ||
      any(!is.finite(data$time))) {
    stopf("`data$time` must contain only finite numeric values.")
  }

  for (variable in vars) {
    values <- data[[variable]]

    if (!is.numeric(values) || anyNA(values) ||
        any(!is.finite(values))) {
      stopf(
        "Exposure variable '", variable,
        "' must contain only finite numeric values."
      )
    }

    if (length(unique(values)) < 2L) {
      stopf(
        "Exposure variable '", variable,
        "' is constant. A cross-basis effect cannot be estimated from a ",
        "constant exposure."
      )
    }
  }

  maximum_lag <- normalize_max_lag(max_lag)
  history_length <- maximum_lag + 1L
  history_contract <-
    "all_fitted_exposures_same_exact_max_lag_plus_one"

  df_var <- positive_integer_scalar(df_var, "df_var")
  df_lag <- positive_integer_scalar(df_lag, "df_lag")

  if (!character_scalar(fun_var)) {
    stopf("`fun_var` must be one non-empty character value.")
  }
  if (!character_scalar(fun_lag)) {
    stopf("`fun_lag` must be one non-empty character value.")
  }

  fun_var <- tolower(fun_var)
  fun_lag <- tolower(fun_lag)

  penalized_functions <- c("ps", "cr")
  if (fun_var %in% penalized_functions ||
      fun_lag %in% penalized_functions) {
    requested <- unique(c(
      if (fun_var %in% penalized_functions) fun_var,
      if (fun_lag %in% penalized_functions) fun_lag
    ))

    stopf(
      "Penalized dlnm basis function(s) ",
      paste(paste0("'", requested, "'"), collapse = ", "),
      " are not supported by the engine-agnostic EpiExposure v1 fitting ",
      "contract. `ps`/`cr` require their penalty matrices to be propagated ",
      "during fitting; EpiExposure does not silently diagnose or fit them as ",
      "unpenalized basis columns."
    )
  }

  allowed_unpenalized <- c("ns", "bs", "poly", "lin")

  if (!fun_var %in% allowed_unpenalized) {
    stopf(
      "Unsupported `fun_var = ", sQuote(fun_var),
      "`. Supported EpiExposure v1 exposure bases are: ",
      paste(allowed_unpenalized, collapse = ", "), "."
    )
  }

  if (!fun_lag %in% allowed_unpenalized) {
    stopf(
      "Unsupported `fun_lag = ", sQuote(fun_lag),
      "`. Supported EpiExposure v1 lag bases are: ",
      paste(allowed_unpenalized, collapse = ", "), "."
    )
  }

  if (maximum_lag == 0L && !identical(fun_lag, "lin")) {
    stopf(
      "When `max_lag = 0`, use `fun_lag = 'lin'`. A non-linear lag basis ",
      "cannot be identified from a single lag value."
    )
  }

  if (!logical_scalar(include_intercept)) {
    stopf("`include_intercept` must be TRUE or FALSE.")
  }
  if (!logical_scalar(keep_design)) {
    stopf("`keep_design` must be TRUE or FALSE.")
  }

  corr_threshold <- positive_numeric_scalar(
    corr_threshold,
    "corr_threshold"
  )
  if (corr_threshold > 1) {
    stopf("`corr_threshold` cannot exceed 1.")
  }

  condition_warn <- positive_numeric_scalar(
    condition_warn,
    "condition_warn"
  )
  condition_severe <- positive_numeric_scalar(
    condition_severe,
    "condition_severe"
  )
  if (condition_severe <= condition_warn) {
    stopf("`condition_severe` must be greater than `condition_warn`.")
  }

  vif_threshold <- positive_numeric_scalar(vif_threshold, "vif_threshold")
  n_per_parameter_warn <- positive_numeric_scalar(
    n_per_parameter_warn,
    "n_per_parameter_warn"
  )
  tol <- positive_numeric_scalar(tol, "tol")

  if (tol >= 1) {
    stopf("`tol` must be smaller than 1.")
  }

  if (!requireNamespace("dlnm", quietly = TRUE)) {
    stopf("Package 'dlnm' is required by `check_identifiability()`.")
  }

  # CANONICAL EPIEXPOSURE DESIGN CONSTRUCTION

  templates <- tryCatch(
    define_exposures(
      data = data,
      vars = vars,
      max_lag = maximum_lag,
      df_var = df_var,
      df_lag = df_lag,
      fun_var = fun_var,
      fun_lag = fun_lag
    ),
    error = function(e) {
      stopf(
        "Could not define the DLNM exposure templates for identifiability ",
        "diagnosis: ", conditionMessage(e)
      )
    }
  )

  design <- tryCatch(
    build_design(
      data = data,
      cb_templates = templates,
      max_lag = maximum_lag,
      include_response = FALSE
    ),
    error = function(e) {
      stopf(
        "Could not construct the epidemic-level DLNM design for ",
        "identifiability diagnosis: ", conditionMessage(e)
      )
    }
  )

  design_max_lag <- attr(design, "epiexposure_max_lag", exact = TRUE)
  design_history_length <- attr(
    design,
    "epiexposure_history_length",
    exact = TRUE
  )
  design_history_contract <- attr(
    design,
    "epiexposure_history_contract",
    exact = TRUE
  )
  design_time_step <- attr(design, "epiexposure_time_step", exact = TRUE)
  design_contract <- attr(
    design,
    "epiexposure_design_contract",
    exact = TRUE
  )
  design_profile_order <- attr(
    design,
    "epiexposure_profile_order",
    exact = TRUE
  )
  basis_specification <- attr(design, "epiexposure_spec", exact = TRUE)
  cb_cols <- attr(design, "epiexposure_cb_cols", exact = TRUE)
  design_vars <- attr(design, "epiexposure_vars", exact = TRUE)

  if (!is.numeric(design_max_lag) || length(design_max_lag) != 1L ||
      is.na(design_max_lag) || !is.finite(design_max_lag) ||
      design_max_lag != as.integer(design_max_lag) ||
      !identical(as.integer(design_max_lag), maximum_lag)) {
    stopf(
      "Internal design metadata disagree with the requested common `max_lag`."
    )
  }

  if (!is.numeric(design_history_length) ||
      length(design_history_length) != 1L ||
      is.na(design_history_length) ||
      !is.finite(design_history_length) ||
      design_history_length != as.integer(design_history_length) ||
      !identical(as.integer(design_history_length), history_length)) {
    stopf(
      "Internal design metadata violate `history_length = max_lag + 1`."
    )
  }

  if (!identical(design_history_contract, history_contract)) {
    stopf(
      "Internal design metadata violate the EpiExposure exact-history contract."
    )
  }

  if (!identical(design_contract, "final_crossbasis_row_per_group")) {
    stopf(
      "Internal design metadata do not use the expected epidemic-level ",
      "final-cross-basis-row contract."
    )
  }

  if (!identical(design_profile_order, "chronological")) {
    stopf("Internal design metadata do not use chronological profile order.")
  }

  if (!identical(as.character(design_vars), as.character(vars))) {
    stopf(
      "Internal design exposure order does not match the requested variables."
    )
  }

  if (!is.character(cb_cols) || !length(cb_cols) ||
      anyNA(cb_cols) || anyDuplicated(cb_cols) ||
      any(!cb_cols %in% names(design))) {
    stopf("Internal cross-basis column metadata are invalid.")
  }

  X_cb <- as.matrix(design[, cb_cols, drop = FALSE])
  storage.mode(X_cb) <- "double"
  rownames(X_cb) <- as.character(design$epi_id)

  if (!nrow(X_cb) || !ncol(X_cb) || anyNA(X_cb) ||
      any(!is.finite(X_cb))) {
    stopf(
      "The canonical epidemic-level cross-basis design must contain only ",
      "finite numeric values."
    )
  }

  n_total <- nrow(X_cb)
  n_complete <- n_total
  n_excluded <- 0L
  complete_ids <- rownames(X_cb)

  X_rank <- if (isTRUE(include_intercept)) {
    cbind(`(Intercept)` = rep(1, nrow(X_cb)), X_cb)
  } else {
    X_cb
  }
  storage.mode(X_rank) <- "double"

  # NUMERICAL-DIAGNOSTIC HELPERS

  column_sd <- function(X) {
    if (!ncol(X)) return(numeric(0))
    if (nrow(X) < 2L) {
      out <- rep(NA_real_, ncol(X))
      names(out) <- colnames(X)
      return(out)
    }
    out <- apply(X, 2L, stats::sd)
    names(out) <- colnames(X)
    out
  }

  near_zero_flags <- function(X) {
    means <- if (ncol(X)) colMeans(X) else numeric(0)
    sds <- column_sd(X)
    cutoff <- tol * pmax(1, abs(means))
    flags <- !is.finite(sds) | sds <= cutoff
    names(flags) <- colnames(X)
    list(flags = flags, mean = means, sd = sds, cutoff = cutoff)
  }

  scaled_matrix_diagnostics <- function(X) {
    if (!is.matrix(X)) X <- as.matrix(X)

    nz <- near_zero_flags(X)
    usable <- !nz$flags

    singular_values <- numeric(0)
    condition_number <- Inf
    max_abs_correlation <- NA_real_
    median_abs_correlation <- NA_real_

    if (sum(usable) == 1L && nrow(X) >= 2L) {
      condition_number <- 1
      max_abs_correlation <- 0
      median_abs_correlation <- 0

      column <- X[, usable, drop = TRUE]
      z <- (column - mean(column)) / stats::sd(column)
      singular_values <- sqrt(sum(z^2))

    } else if (sum(usable) > 1L && nrow(X) >= 2L) {
      X_use <- X[, usable, drop = FALSE]
      means <- colMeans(X_use)
      sds <- column_sd(X_use)
      X_scaled <- sweep(X_use, 2L, means, FUN = "-")
      X_scaled <- sweep(X_scaled, 2L, sds, FUN = "/")

      sv <- tryCatch(
        svd(X_scaled, nu = 0L, nv = 0L)$d,
        error = function(e) numeric(0)
      )
      singular_values <- sv

      scaled_rank <- qr(
        X_scaled,
        tol = tol,
        LAPACK = FALSE
      )$rank

      if (nrow(X_scaled) > ncol(X_scaled) &&
          scaled_rank == ncol(X_scaled) &&
          length(sv) == ncol(X_scaled) &&
          max(sv) > 0 &&
          min(sv) > tol * max(sv)) {
        condition_number <- max(sv) / min(sv)
      }

      R <- suppressWarnings(stats::cor(X_use))
      upper_values <- abs(R[upper.tri(R)])
      upper_values <- upper_values[is.finite(upper_values)]

      if (length(upper_values)) {
        max_abs_correlation <- max(upper_values)
        median_abs_correlation <- stats::median(upper_values)
      }
    }

    list(
      near_zero = nz$flags,
      sd = nz$sd,
      mean = nz$mean,
      cutoff = nz$cutoff,
      usable = usable,
      singular_values = singular_values,
      condition_number = condition_number,
      max_abs_correlation = max_abs_correlation,
      median_abs_correlation = median_abs_correlation
    )
  }

  numerical_rank <- function(X) {
    qr(X, tol = tol, LAPACK = FALSE)$rank
  }

  compute_vif <- function(X, nz_flags) {
    p <- ncol(X)
    out <- rep(NA_real_, p)
    names(out) <- colnames(X)

    if (!p || nrow(X) < 2L) return(out)

    usable <- which(!nz_flags)
    if (!length(usable)) return(out)

    for (j in usable) {
      y <- X[, j]
      y_centered <- y - mean(y)
      sst <- sum(y_centered^2)

      if (!is.finite(sst) || sst <= tol * max(1, sum(y^2))) {
        out[j] <- NA_real_
        next
      }

      others <- setdiff(usable, j)
      if (!length(others)) {
        out[j] <- 1
        next
      }

      Z <- cbind(`(Intercept)` = 1, X[, others, drop = FALSE])
      fit_j <- tryCatch(stats::lm.fit(x = Z, y = y), error = function(e) NULL)

      if (is.null(fit_j) || is.null(fit_j$residuals) ||
          any(!is.finite(fit_j$residuals))) {
        out[j] <- NA_real_
        next
      }

      sse <- sum(fit_j$residuals^2)
      ratio <- sse / sst

      if (!is.finite(ratio)) {
        out[j] <- NA_real_
      } else if (ratio <= tol) {
        out[j] <- Inf
      } else {
        out[j] <- 1 / ratio
      }
    }

    out
  }

  # COMBINED RANK AND CONDITION DIAGNOSTICS

  qr_full <- qr(X_rank, tol = tol, LAPACK = FALSE)
  rank_full <- qr_full$rank
  p_full <- ncol(X_rank)
  full_rank <- identical(rank_full, p_full)

  dependent_columns <- character(0)
  if (!full_rank) {
    dependent_positions <- qr_full$pivot[
      seq.int(rank_full + 1L, p_full)
    ]
    dependent_columns <- colnames(X_rank)[dependent_positions]
  }

  combined_scaled <- scaled_matrix_diagnostics(X_cb)
  nzv_columns <- names(combined_scaled$near_zero)[
    combined_scaled$near_zero
  ]
  singular_values <- combined_scaled$singular_values
  condition_number <- combined_scaled$condition_number
  max_abs_cb_correlation <- combined_scaled$max_abs_correlation
  median_abs_cb_correlation <- combined_scaled$median_abs_correlation

  min_singular_value <- if (length(singular_values)) {
    min(singular_values)
  } else {
    NA_real_
  }

  max_singular_value <- if (length(singular_values)) {
    max(singular_values)
  } else {
    NA_real_
  }

  n_per_parameter <- n_complete / p_full
  design_residual_df_proxy <- n_complete - rank_full

  # SUPPLEMENTARY VIFS

  cb_variable <- rep(NA_character_, length(cb_cols))
  names(cb_variable) <- cb_cols

  for (variable in vars) {
    template <- templates[[variable]]

    if (is.null(template) || !inherits(template, "crossbasis") ||
        !ncol(template)) {
      stopf(
        "Internal training-template mapping failed for exposure '",
        variable, "'."
      )
    }

    variable_cols <- paste0(
      "cb_",
      variable,
      "_",
      seq_len(ncol(template))
    )

    if (!all(variable_cols %in% cb_cols)) {
      stopf(
        "Internal cross-basis column mapping failed for exposure '",
        variable, "'."
      )
    }

    cb_variable[variable_cols] <- variable
  }

  if (anyNA(cb_variable)) {
    stopf("Some cross-basis columns could not be mapped to an exposure variable.")
  }

  vif_values <- compute_vif(
    X_cb,
    combined_scaled$near_zero
  )

  vif_table <- data.frame(
    column = cb_cols,
    variable = unname(cb_variable[cb_cols]),
    vif = unname(vif_values[cb_cols]),
    stringsAsFactors = FALSE
  )

  nonmissing_vif <- vif_table$vif[!is.na(vif_table$vif)]
  max_vif <- if (length(nonmissing_vif)) {
    max(nonmissing_vif)
  } else {
    NA_real_
  }
  median_vif <- if (length(nonmissing_vif)) {
    stats::median(nonmissing_vif)
  } else {
    NA_real_
  }

  # VARIABLE-SPECIFIC CROSS-BASIS DIAGNOSTICS

  by_variable <- do.call(
    rbind,
    lapply(vars, function(variable) {
      variable_cols <- names(cb_variable)[cb_variable == variable]
      B <- X_cb[, variable_cols, drop = FALSE]

      block_rank <- numerical_rank(B)
      block_scaled <- scaled_matrix_diagnostics(B)
      exposure_values <- data[[variable]]

      data.frame(
        variable = variable,
        n_unique_exposure = length(unique(exposure_values)),
        exposure_missing_n = 0L,
        exposure_missing_percent = 0,
        history_length = history_length,
        max_lag = maximum_lag,
        basis_columns = ncol(B),
        rank = block_rank,
        full_rank = block_rank == ncol(B),
        rank_ratio = if (ncol(B) > 0L) block_rank / ncol(B) else NA_real_,
        condition_number_scaled = block_scaled$condition_number,
        max_abs_within_basis_correlation =
          block_scaled$max_abs_correlation,
        near_zero_variance_columns = sum(block_scaled$near_zero),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(by_variable) <- NULL

  # DESCRIPTIVE RAW-EXPOSURE CORRELATIONS

  pairwise_exposure_correlation <- data.frame(
    variable_1 = character(0),
    variable_2 = character(0),
    n_complete = integer(0),
    pearson_correlation = numeric(0),
    spearman_correlation = numeric(0),
    stringsAsFactors = FALSE
  )

  if (length(vars) >= 2L) {
    pairs <- utils::combn(vars, 2L, simplify = FALSE)

    pairwise_exposure_correlation <- do.call(
      rbind,
      lapply(pairs, function(pair) {
        x <- data[[pair[1L]]]
        y <- data[[pair[2L]]]
        ok <- is.finite(x) & is.finite(y)
        n_ok <- sum(ok)

        pearson <- if (
          n_ok >= 3L &&
          stats::sd(x[ok]) > 0 &&
          stats::sd(y[ok]) > 0
        ) {
          stats::cor(x[ok], y[ok], method = "pearson")
        } else {
          NA_real_
        }

        spearman <- if (n_ok >= 3L) {
          suppressWarnings(
            stats::cor(x[ok], y[ok], method = "spearman")
          )
        } else {
          NA_real_
        }

        data.frame(
          variable_1 = pair[1L],
          variable_2 = pair[2L],
          n_complete = n_ok,
          pearson_correlation = pearson,
          spearman_correlation = spearman,
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(pairwise_exposure_correlation) <- NULL
  }

  # BETWEEN-EXPOSURE CROSS-BASIS CORRELATIONS

  pairwise_crossbasis_correlation <- data.frame(
    variable_1 = character(0),
    variable_2 = character(0),
    max_abs_correlation = numeric(0),
    mean_abs_correlation = numeric(0),
    column_1 = character(0),
    column_2 = character(0),
    stringsAsFactors = FALSE
  )

  if (length(vars) >= 2L) {
    pairs <- utils::combn(vars, 2L, simplify = FALSE)

    pairwise_crossbasis_correlation <- do.call(
      rbind,
      lapply(pairs, function(pair) {
        A_cols <- names(cb_variable)[cb_variable == pair[1L]]
        B_cols <- names(cb_variable)[cb_variable == pair[2L]]

        A <- X_cb[, A_cols, drop = FALSE]
        B <- X_cb[, B_cols, drop = FALSE]

        A_nz <- near_zero_flags(A)$flags
        B_nz <- near_zero_flags(B)$flags

        A <- A[, !A_nz, drop = FALSE]
        B <- B[, !B_nz, drop = FALSE]

        if (!ncol(A) || !ncol(B) || nrow(X_cb) < 2L) {
          return(data.frame(
            variable_1 = pair[1L],
            variable_2 = pair[2L],
            max_abs_correlation = NA_real_,
            mean_abs_correlation = NA_real_,
            column_1 = NA_character_,
            column_2 = NA_character_,
            stringsAsFactors = FALSE
          ))
        }

        C <- suppressWarnings(stats::cor(A, B))
        abs_C <- abs(C)
        finite_mask <- is.finite(abs_C)

        if (!any(finite_mask)) {
          return(data.frame(
            variable_1 = pair[1L],
            variable_2 = pair[2L],
            max_abs_correlation = NA_real_,
            mean_abs_correlation = NA_real_,
            column_1 = NA_character_,
            column_2 = NA_character_,
            stringsAsFactors = FALSE
          ))
        }

        maximum <- max(abs_C[finite_mask])
        candidates <- which(finite_mask & abs_C == maximum, arr.ind = TRUE)
        max_index <- candidates[1L, , drop = FALSE]

        data.frame(
          variable_1 = pair[1L],
          variable_2 = pair[2L],
          max_abs_correlation = maximum,
          mean_abs_correlation = mean(abs_C[finite_mask]),
          column_1 = rownames(C)[max_index[1L, 1L]],
          column_2 = colnames(C)[max_index[1L, 2L]],
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(pairwise_crossbasis_correlation) <- NULL
  }

  # STATUS FLAGS

  severe_condition <- !is.finite(condition_number) ||
    condition_number >= condition_severe

  warning_condition <- is.finite(condition_number) &&
    condition_number >= condition_warn &&
    condition_number < condition_severe

  high_crossbasis_pair <- nrow(pairwise_crossbasis_correlation) > 0L &&
    any(
      is.finite(pairwise_crossbasis_correlation$max_abs_correlation) &
        pairwise_crossbasis_correlation$max_abs_correlation >=
        corr_threshold
    )

  high_raw_pair <- nrow(pairwise_exposure_correlation) > 0L &&
    any(
      is.finite(pairwise_exposure_correlation$pearson_correlation) &
        abs(pairwise_exposure_correlation$pearson_correlation) >=
        corr_threshold
    )

  high_vif <- nrow(vif_table) > 0L &&
    any(!is.na(vif_table$vif) & vif_table$vif >= vif_threshold)

  low_n_per_parameter <- n_per_parameter < n_per_parameter_warn
  no_residual_design_df <- n_complete <= p_full

  identifiable <- full_rank
  numerically_stable <- full_rank &&
    !severe_condition &&
    !length(nzv_columns) &&
    !no_residual_design_df

  warning_status <- warning_condition ||
    high_crossbasis_pair ||
    low_n_per_parameter

  status <- if (!identifiable || !numerically_stable) {
    "problem"
  } else if (warning_status) {
    "warning"
  } else {
    "ok"
  }

  diagnostic_flags <- c(
    rank_deficient = !full_rank,
    near_zero_variance = length(nzv_columns) > 0L,
    condition_warning = warning_condition,
    condition_severe = severe_condition,
    high_between_exposure_crossbasis_correlation = high_crossbasis_pair,
    high_raw_exposure_correlation_supplementary = high_raw_pair,
    high_vif_supplementary = high_vif,
    low_epidemics_per_design_column = low_n_per_parameter,
    no_positive_design_residual_df = no_residual_design_df
  )

  # RECOMMENDATIONS

  recommendations <- character(0)

  if (!full_rank) {
    recommendations <- c(
      recommendations,
      paste0(
        "The combined epidemic-level design is numerically rank-deficient ",
        "at tol = ", format(tol, scientific = TRUE), " (rank = ",
        rank_full, ", columns = ", p_full, "). Inspect the reported ",
        "dependent columns and simplify the exposure-lag basis. Reducing ",
        "`df_var` and/or `df_lag` is usually preferable before shortening ",
        "`max_lag`; change the lag window only when scientifically justified."
      )
    )
  }

  if (length(nzv_columns)) {
    recommendations <- c(
      recommendations,
      paste0(
        "Near-zero-variance cross-basis columns were detected: ",
        paste(utils::head(nzv_columns, 8L), collapse = ", "),
        if (length(nzv_columns) > 8L) ", ..." else "",
        ". Inspect exposure support and basis complexity."
      )
    )
  }

  if (severe_condition) {
    recommendations <- c(
      recommendations,
      paste0(
        "The centered/scaled cross-basis condition number is ",
        if (is.finite(condition_number)) {
          format(condition_number, digits = 4)
        } else {
          "infinite"
        },
        ", indicating severe numerical instability under the selected ",
        "threshold. Consider reducing basis complexity and checking whether ",
        "exposures provide redundant transformed temporal information."
      )
    )
  } else if (warning_condition) {
    recommendations <- c(
      recommendations,
      paste0(
        "The centered/scaled cross-basis condition number is ",
        format(condition_number, digits = 4),
        ". This exceeds the warning threshold; inspect the between-exposure ",
        "cross-basis diagnostics and basis complexity before fitting."
      )
    )
  }

  if (high_crossbasis_pair) {
    flagged <- pairwise_crossbasis_correlation[
      is.finite(pairwise_crossbasis_correlation$max_abs_correlation) &
        pairwise_crossbasis_correlation$max_abs_correlation >= corr_threshold,
      ,
      drop = FALSE
    ]

    flagged_text <- paste(
      paste0(
        flagged$variable_1, " vs ", flagged$variable_2,
        " (max |r| = ",
        format(flagged$max_abs_correlation, digits = 3), ")"
      ),
      collapse = "; "
    )

    recommendations <- c(
      recommendations,
      paste0(
        "Strong correlation was detected between different exposure ",
        "cross-basis blocks: ", flagged_text,
        ". Consider whether these exposures provide redundant temporal ",
        "information and compare simpler or alternative variable sets."
      )
    )
  }

  if (high_raw_pair) {
    flagged <- pairwise_exposure_correlation[
      is.finite(pairwise_exposure_correlation$pearson_correlation) &
        abs(pairwise_exposure_correlation$pearson_correlation) >=
        corr_threshold,
      ,
      drop = FALSE
    ]

    flagged_text <- paste(
      paste0(
        flagged$variable_1, " vs ", flagged$variable_2,
        " (r = ",
        format(flagged$pearson_correlation, digits = 3), ")"
      ),
      collapse = "; "
    )

    recommendations <- c(
      recommendations,
      paste0(
        "High same-time correlation was detected among raw exposures: ",
        flagged_text,
        ". This is descriptive only and does not by itself determine the ",
        "global status; the fitted cross-basis rank and conditioning are more ",
        "directly relevant to DLNM identifiability."
      )
    )
  }

  if (high_vif) {
    flagged_vif <- vif_table[
      !is.na(vif_table$vif) & vif_table$vif >= vif_threshold,
      ,
      drop = FALSE
    ]

    recommendations <- c(
      recommendations,
      paste0(
        nrow(flagged_vif),
        " cross-basis column(s) have VIF >= ", vif_threshold,
        ". Treat this as supplementary because spline-basis columns are ",
        "correlated by construction. High VIF alone does not change `status`; ",
        "prioritize the combined rank and scaled condition number."
      )
    )
  }

  if (no_residual_design_df) {
    recommendations <- c(
      recommendations,
      paste0(
        "The design contains ", n_complete, " epidemic(s) and ", p_full,
        " rank-design column(s), leaving no positive design residual degrees ",
        "of freedom proxy. Even a formally full-rank square design is not ",
        "classified as numerically stable. Simplify the basis or increase the ",
        "number of independent epidemics."
      )
    )
  }

  if (low_n_per_parameter) {
    recommendations <- c(
      recommendations,
      paste0(
        "There are ", format(n_per_parameter, digits = 3),
        " epidemics per rank-design column, below the diagnostic threshold of ",
        format(n_per_parameter_warn, digits = 3),
        ". This is a heuristic overparameterization signal rather than a ",
        "universal sample-size rule. Consider reducing basis complexity or ",
        "increasing the number of independent epidemics."
      )
    )
  }

  low_unique <- if (identical(fun_var, "lin")) {
    rep(FALSE, nrow(by_variable))
  } else {
    by_variable$n_unique_exposure <= df_var
  }

  if (any(low_unique)) {
    recommendations <- c(
      recommendations,
      paste0(
        "Limited exposure support relative to the requested exposure-basis ",
        "complexity was detected for: ",
        paste(by_variable$variable[low_unique], collapse = ", "),
        ". Consider reducing exposure-basis flexibility."
      )
    )
  }

  if (!length(recommendations)) {
    recommendations <- paste0(
      "No major identifiability or numerical-stability problem was detected ",
      "using the current diagnostic thresholds."
    )
  }

  # OUTPUT

  overall <- data.frame(
    n_epidemics_total = n_total,
    n_epidemics_complete = n_complete,
    n_epidemics_excluded = n_excluded,
    n_exposures = length(vars),
    max_lag = maximum_lag,
    history_length = history_length,
    crossbasis_columns = ncol(X_cb),
    design_columns = p_full,
    rank = rank_full,
    full_rank = full_rank,
    rank_ratio = if (p_full > 0L) rank_full / p_full else NA_real_,
    design_residual_df_proxy = design_residual_df_proxy,
    condition_number_scaled = condition_number,
    min_singular_value_scaled = min_singular_value,
    max_singular_value_scaled = max_singular_value,
    max_abs_crossbasis_correlation = max_abs_cb_correlation,
    median_abs_crossbasis_correlation = median_abs_cb_correlation,
    max_vif = max_vif,
    median_vif = median_vif,
    near_zero_variance_columns = length(nzv_columns),
    epidemics_per_design_column = n_per_parameter,
    time_step = if (is.null(design_time_step)) NA_real_ else design_time_step,
    stringsAsFactors = FALSE
  )

  out <- list(
    identifiable = identifiable,
    numerically_stable = numerically_stable,
    status = status,
    overall = overall,
    by_variable = by_variable,
    pairwise_exposure_correlation = pairwise_exposure_correlation,
    pairwise_crossbasis_correlation = pairwise_crossbasis_correlation,
    vif = vif_table,
    singular_values = singular_values,
    dependent_columns = dependent_columns,
    near_zero_variance_columns = nzv_columns,
    complete_epi_id = complete_ids,
    basis_specification = basis_specification,
    diagnostic_flags = diagnostic_flags,
    recommendations = unique(recommendations),
    settings = list(
      variables = vars,
      max_lag = maximum_lag,
      history_length = history_length,
      history_contract = history_contract,
      design_contract = design_contract,
      profile_order = design_profile_order,
      time_step = if (is.null(design_time_step)) NA_real_ else design_time_step,
      df_var = df_var,
      df_lag = df_lag,
      fun_var = fun_var,
      fun_lag = fun_lag,
      penalized = FALSE,
      include_intercept = include_intercept,
      corr_threshold = corr_threshold,
      condition_warn = condition_warn,
      condition_severe = condition_severe,
      vif_threshold = vif_threshold,
      n_per_parameter_warn = n_per_parameter_warn,
      tol = tol,
      design_source = "define_exposures_plus_build_design"
    )
  )

  if (isTRUE(keep_design)) {
    out$design_matrix <- X_rank
  }

  attr(out, "epiexposure_max_lag") <- maximum_lag
  attr(out, "epiexposure_history_length") <- history_length
  attr(out, "epiexposure_history_contract") <- history_contract
  attr(out, "epiexposure_time_step") <-
    if (is.null(design_time_step)) NA_real_ else design_time_step
  attr(out, "epiexposure_design_contract") <- design_contract
  attr(out, "epiexposure_profile_order") <- design_profile_order

  class(out) <- c("epiexposure_identifiability", "list")
  out
}


#' Print DLNM identifiability diagnostics
#'
#' @param x Object returned by `check_identifiability()`.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.epiexposure_identifiability <- function(x, ...) {
  cat("EpiExposure DLNM identifiability diagnostics\n")
  cat("Status: ", toupper(x$status), "\n", sep = "")
  cat("Identifiable (full numerical rank): ", x$identifiable, "\n", sep = "")
  cat("Numerically stable: ", x$numerically_stable, "\n\n", sep = "")

  cat("Overall design\n")
  print(x$overall, row.names = FALSE)

  cat("\nBy exposure variable\n")
  print(x$by_variable, row.names = FALSE)

  if (nrow(x$pairwise_crossbasis_correlation)) {
    cat("\nBetween-exposure cross-basis correlation\n")
    print(x$pairwise_crossbasis_correlation, row.names = FALSE)
  }

  if (nrow(x$pairwise_exposure_correlation)) {
    cat("\nRaw exposure correlation (descriptive)\n")
    print(x$pairwise_exposure_correlation, row.names = FALSE)
  }

  if (length(x$dependent_columns)) {
    cat("\nPivot-based dependent columns\n")
    cat(paste0("- ", x$dependent_columns), sep = "\n")
    cat("\n")
  }

  if (length(x$near_zero_variance_columns)) {
    cat("\nNear-zero-variance cross-basis columns\n")
    cat(paste0("- ", x$near_zero_variance_columns), sep = "\n")
    cat("\n")
  }

  cat("\nRecommendations\n")
  cat(paste0("- ", x$recommendations), sep = "\n")
  cat("\n")

  invisible(x)
}
