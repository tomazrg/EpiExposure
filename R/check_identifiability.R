#' Diagnose identifiability and numerical stability of DLNM cross-basis designs
#'
#' Evaluates the epidemic-level DLNM design matrix that would be used for model
#' fitting. One or more exposure variables can be assessed simultaneously.
#'
#' The function reports exact rank, rank deficiency, scaled condition number,
#' singular values, near-zero-variance columns, cross-basis correlations,
#' descriptive VIFs, pairwise raw-exposure correlations, pairwise correlation
#' between cross-basis blocks, and the ratio between complete epidemics and
#' model columns. It also returns concise recommendations.
#'
#' @param data Long-format exposure data containing `epi_id`, `time`, and the
#'   exposure variables.
#' @param var Character vector with one or more exposure-variable names.
#' @param max_lag Maximum retrospective lag. A vector such as `c(0, 85)` is
#'   accepted and reduced to its maximum value for backward compatibility.
#' @param df_var Positive integer degrees of freedom for the exposure basis.
#' @param df_lag Positive integer degrees of freedom for the lag basis.
#' @param fun_var Exposure-basis function: `"ns"`, `"bs"`, `"poly"`, or `"lin"`.
#' @param fun_lag Lag-basis function: `"ns"`, `"ps"`, or `"lin"`.
#' @param include_intercept Logical. Include a model intercept when checking the
#'   rank of the combined design matrix. Default is `TRUE`.
#' @param corr_threshold Heuristic absolute-correlation threshold used to flag
#'   potentially important collinearity. Default is 0.90.
#' @param condition_warn Heuristic scaled condition-number threshold for a
#'   warning. Default is 30.
#' @param condition_severe Heuristic scaled condition-number threshold for a
#'   severe numerical-stability warning. Default is 100.
#' @param vif_threshold Heuristic VIF threshold used only as a supplementary
#'   diagnostic. Spline-basis columns are often correlated by construction.
#' @param n_per_parameter_warn Heuristic minimum ratio of complete epidemics to
#'   design columns. Default is 10.
#' @param tol Numerical tolerance used in QR-rank and near-zero-variance checks.
#' @param keep_design Logical. If `TRUE`, include the complete epidemic-level
#'   cross-basis design matrix in the returned object.
#'
#' @return An object of class `epiexposure_identifiability` containing:
#'   `identifiable`, `numerically_stable`, `status`, `overall`, `by_variable`,
#'   `pairwise_exposure_correlation`, `pairwise_crossbasis_correlation`, `vif`,
#'   `dependent_columns`, `recommendations`, and optionally `design_matrix`.
#'
#' @details
#' Mathematical identifiability is defined here by full column rank of the
#' epidemic-level design matrix. Numerical stability is evaluated separately.
#' The condition-number, correlation, VIF, and observations-per-column cutoffs
#' are practical diagnostics rather than universal inferential thresholds.
#'
#' Individual VIFs should be interpreted cautiously for spline-expanded terms,
#' because correlations among columns belonging to the same basis can be
#' expected. The combined rank and scaled condition number are the primary
#' numerical diagnostics.
#'
#' @export
check_identifiability <- function(
    data,
    var,
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

  # ---------------------------------------------------------------------------
  # Basic validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  required_columns <- c("epi_id", "time")
  missing_required <- setdiff(required_columns, names(data))
  if (length(missing_required)) {
    stop(
      "`data` is missing required columns: ",
      paste(missing_required, collapse = ", "), "."
    )
  }

  if (!is.character(var) || !length(var) || anyNA(var) ||
      any(!nzchar(var)) || anyDuplicated(var)) {
    stop("`var` must contain unique, non-empty exposure-variable names.")
  }

  missing_variables <- setdiff(var, names(data))
  if (length(missing_variables)) {
    stop(
      "Exposure variables not found in `data`: ",
      paste(missing_variables, collapse = ", "), "."
    )
  }

  for (v in var) {
    if (!is.numeric(data[[v]])) {
      stop("Exposure variable '", v, "' must be numeric.")
    }
    if (all(!is.finite(data[[v]]))) {
      stop("Exposure variable '", v, "' contains no finite values.")
    }
  }

  if (!is.numeric(max_lag) || !length(max_lag) || anyNA(max_lag) ||
      any(!is.finite(max_lag))) {
    stop("`max_lag` must contain finite numeric values.")
  }
  max_lag <- max(max_lag)
  if (max_lag < 0 || abs(max_lag - round(max_lag)) > tol) {
    stop("`max_lag` must resolve to one non-negative integer.")
  }
  max_lag <- as.integer(round(max_lag))

  is_positive_integer_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
      x > 0 && abs(x - round(x)) <= tol
  }

  if (!is_positive_integer_scalar(df_var)) {
    stop("`df_var` must be one positive integer.")
  }
  if (!is_positive_integer_scalar(df_lag)) {
    stop("`df_lag` must be one positive integer.")
  }
  df_var <- as.integer(round(df_var))
  df_lag <- as.integer(round(df_lag))

  if (!is.character(fun_var) || length(fun_var) != 1L ||
      !fun_var %in% c("ns", "bs", "poly", "lin")) {
    stop("`fun_var` must be one of 'ns', 'bs', 'poly', or 'lin'.")
  }
  if (!is.character(fun_lag) || length(fun_lag) != 1L ||
      !fun_lag %in% c("ns", "ps", "lin")) {
    stop("`fun_lag` must be one of 'ns', 'ps', or 'lin'.")
  }

  logical_scalar <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!logical_scalar(include_intercept)) {
    stop("`include_intercept` must be TRUE or FALSE.")
  }
  if (!logical_scalar(keep_design)) {
    stop("`keep_design` must be TRUE or FALSE.")
  }

  numeric_positive_scalar <- function(x, name, allow_zero = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
      if (allow_zero) x >= 0 else x > 0
    if (!ok) stop("`", name, "` must be one finite ",
                  if (allow_zero) "non-negative" else "positive", " number.")
    as.numeric(x)
  }

  corr_threshold <- numeric_positive_scalar(corr_threshold, "corr_threshold")
  if (corr_threshold > 1) stop("`corr_threshold` cannot exceed 1.")
  condition_warn <- numeric_positive_scalar(condition_warn, "condition_warn")
  condition_severe <- numeric_positive_scalar(condition_severe, "condition_severe")
  if (condition_severe <= condition_warn) {
    stop("`condition_severe` must be greater than `condition_warn`.")
  }
  vif_threshold <- numeric_positive_scalar(vif_threshold, "vif_threshold")
  n_per_parameter_warn <- numeric_positive_scalar(
    n_per_parameter_warn, "n_per_parameter_warn"
  )
  tol <- numeric_positive_scalar(tol, "tol")

  if (!requireNamespace("dlnm", quietly = TRUE)) {
    stop("Package 'dlnm' is required by `check_identifiability()`.")
  }

  # ---------------------------------------------------------------------------
  # Temporal integrity
  # ---------------------------------------------------------------------------
  ids <- unique(data$epi_id)
  if (!length(ids)) stop("`data` contains no epidemics.")

  to_numeric_time <- function(x) {
    if (inherits(x, "Date") || inherits(x, "POSIXt")) return(as.numeric(x))
    if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
    stop("`time` must be numeric/integer, Date, or POSIXt.")
  }

  time_steps <- rep(NA_real_, length(ids))
  names(time_steps) <- as.character(ids)
  insufficient_ids <- character(0)
  irregular_ids <- character(0)
  duplicated_ids <- character(0)
  missing_time_ids <- character(0)

  for (i in seq_along(ids)) {
    current_id <- ids[i]
    idx <- which(data$epi_id == current_id)
    tt_raw <- data$time[idx]

    if (anyNA(tt_raw)) {
      missing_time_ids <- c(missing_time_ids, as.character(current_id))
      next
    }

    tt <- to_numeric_time(tt_raw)
    if (anyDuplicated(tt)) {
      duplicated_ids <- c(duplicated_ids, as.character(current_id))
      next
    }

    tt <- sort(tt)
    if (length(tt) < max_lag + 1L) {
      insufficient_ids <- c(insufficient_ids, as.character(current_id))
    }

    if (length(tt) > 1L) {
      dd <- diff(tt)
      if (any(dd <= 0) ||
          max(abs(dd - dd[1L])) > tol * max(1, abs(dd[1L]))) {
        irregular_ids <- c(irregular_ids, as.character(current_id))
      } else {
        time_steps[i] <- dd[1L]
      }
    }
  }

  if (length(missing_time_ids)) {
    stop(
      "Missing `time` values were detected. Example epi_id: ",
      paste(utils::head(missing_time_ids, 5L), collapse = ", "), "."
    )
  }
  if (length(duplicated_ids)) {
    stop(
      "Duplicated `time` values were detected within epidemics. Example epi_id: ",
      paste(utils::head(duplicated_ids, 5L), collapse = ", "), "."
    )
  }
  if (length(insufficient_ids)) {
    stop(
      "Some epidemics do not have enough temporal coverage for max_lag = ",
      max_lag, ". At least ", max_lag + 1L,
      " observations are required. Example epi_id: ",
      paste(utils::head(insufficient_ids, 5L), collapse = ", "), "."
    )
  }
  if (length(irregular_ids)) {
    stop(
      "Irregular time spacing was detected within epidemics. DLNM lags assume ",
      "regularly spaced observations. Example epi_id: ",
      paste(utils::head(irregular_ids, 5L), collapse = ", "), "."
    )
  }

  finite_steps <- time_steps[is.finite(time_steps)]
  if (length(finite_steps) > 1L &&
      (max(finite_steps) - min(finite_steps)) >
      tol * max(1, max(abs(finite_steps)))) {
    stop("Different temporal step sizes were detected across epidemics.")
  }
  time_step <- if (length(finite_steps)) finite_steps[1L] else NA_real_

  # ---------------------------------------------------------------------------
  # Basis specification
  # ---------------------------------------------------------------------------
  argvar <- switch(
    fun_var,
    ns   = list(fun = "ns", df = df_var),
    bs   = list(fun = "bs", df = df_var),
    poly = list(fun = "poly", degree = df_var),
    lin  = list(fun = "lin")
  )

  if (!is.null(argvar$fun) && argvar$fun != "lin") {
    argvar$intercept <- FALSE
  }

  arglag <- switch(
    fun_lag,
    ns  = list(fun = "ns", df = df_lag),
    ps  = list(fun = "ps", df = df_lag),
    lin = list(fun = "lin")
  )

  # ---------------------------------------------------------------------------
  # Helpers reproducing the epidemic-level design used for model fitting
  # ---------------------------------------------------------------------------
  build_pooled_series <- function(dat, variable, separator_n) {
    out <- vector("list", length(ids))
    for (i in seq_along(ids)) {
      current_id <- ids[i]
      idx <- which(dat$epi_id == current_id)
      ord <- order(dat$time[idx])
      values <- dat[[variable]][idx][ord]
      out[[i]] <- c(values, rep(NA_real_, separator_n))
    }
    unlist(out, use.names = FALSE)
  }

  extract_last_cb <- function(values, template) {
    if (length(values) <= max_lag) {
      return(rep(NA_real_, ncol(template)))
    }

    cb <- dlnm::crossbasis(
      values,
      lag = max_lag,
      argvar = attr(template, "argvar"),
      arglag = attr(template, "arglag")
    )

    as.numeric(cb[length(values), , drop = TRUE])
  }

  templates <- vector("list", length(var))
  names(templates) <- var
  blocks <- vector("list", length(var))
  names(blocks) <- var

  for (v in var) {
    pooled <- build_pooled_series(data, v, max_lag)

    template <- tryCatch(
      dlnm::crossbasis(
        pooled,
        lag = max_lag,
        argvar = argvar,
        arglag = arglag
      ),
      error = function(e) {
        stop(
          "Could not construct the cross-basis for variable '", v,
          "': ", conditionMessage(e), call. = FALSE
        )
      }
    )

    templates[[v]] <- template

    block <- matrix(
      NA_real_,
      nrow = length(ids),
      ncol = ncol(template)
    )

    for (i in seq_along(ids)) {
      current_id <- ids[i]
      idx <- which(data$epi_id == current_id)
      ord <- order(data$time[idx])
      values <- data[[v]][idx][ord]
      block[i, ] <- extract_last_cb(values, template)
    }

    colnames(block) <- paste0("cb_", v, "_", seq_len(ncol(block)))
    rownames(block) <- as.character(ids)
    blocks[[v]] <- block
  }

  X_cb_all <- do.call(cbind, blocks)
  complete_rows <- stats::complete.cases(X_cb_all)
  n_complete <- sum(complete_rows)
  n_total <- nrow(X_cb_all)
  n_excluded <- n_total - n_complete

  if (n_complete == 0L) {
    stop(
      "No epidemic has a complete cross-basis design across all requested ",
      "variables. Inspect missing exposure values within the lag histories."
    )
  }

  X_cb <- X_cb_all[complete_rows, , drop = FALSE]
  complete_ids <- rownames(X_cb)

  X_rank <- if (include_intercept) {
    cbind(`(Intercept)` = 1, X_cb)
  } else {
    X_cb
  }

  # ---------------------------------------------------------------------------
  # Matrix diagnostics
  # ---------------------------------------------------------------------------
  qr_full <- qr(X_rank, tol = tol)
  rank_full <- qr_full$rank
  p_full <- ncol(X_rank)
  full_rank <- rank_full == p_full

  dependent_columns <- character(0)
  if (!full_rank) {
    dependent_columns <- colnames(X_rank)[
      qr_full$pivot[seq.int(rank_full + 1L, p_full)]
    ]
  }

  predictor_sd <- apply(X_cb, 2L, stats::sd)
  predictor_mean <- colMeans(X_cb)
  nzv_cutoff <- tol * pmax(1, abs(predictor_mean))
  near_zero_variance <- !is.finite(predictor_sd) | predictor_sd <= nzv_cutoff
  nzv_columns <- names(predictor_sd)[near_zero_variance]

  usable_columns <- !near_zero_variance
  X_scaled <- NULL
  singular_values <- numeric(0)
  condition_number <- Inf
  max_abs_cb_correlation <- NA_real_
  median_abs_cb_correlation <- NA_real_

  if (sum(usable_columns) == 1L) {
    X_scaled <- scale(X_cb[, usable_columns, drop = FALSE])
    singular_values <- sqrt(sum(X_scaled^2))
    condition_number <- 1
    max_abs_cb_correlation <- 0
    median_abs_cb_correlation <- 0
  } else if (sum(usable_columns) > 1L) {
    X_scaled <- scale(X_cb[, usable_columns, drop = FALSE])
    sv <- svd(X_scaled, nu = 0L, nv = 0L)$d
    singular_values <- sv

    if (nrow(X_scaled) < ncol(X_scaled) ||
        qr(X_scaled, tol = tol)$rank < ncol(X_scaled) ||
        !length(sv) || max(sv) == 0 || min(sv) <= tol * max(sv)) {
      condition_number <- Inf
    } else {
      condition_number <- max(sv) / min(sv)
    }

    R_cb <- stats::cor(X_cb[, usable_columns, drop = FALSE])
    upper_values <- abs(R_cb[upper.tri(R_cb)])
    if (length(upper_values)) {
      max_abs_cb_correlation <- max(upper_values, na.rm = TRUE)
      median_abs_cb_correlation <- stats::median(upper_values, na.rm = TRUE)
    }
  }

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

  # ---------------------------------------------------------------------------
  # Descriptive VIFs for cross-basis columns
  # ---------------------------------------------------------------------------
  vif_variable <- unlist(
    lapply(var, function(v) rep(v, ncol(blocks[[v]]))),
    use.names = FALSE
  )
  vif_table <- data.frame(
    column = colnames(X_cb),
    variable = vif_variable,
    vif = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!any(near_zero_variance)) {
    if (ncol(X_cb) == 1L) {
      vif_table$vif <- 1
    } else if (n_complete > ncol(X_cb)) {
      R <- stats::cor(X_cb)
      invR <- tryCatch(solve(R), error = function(e) NULL)
      if (!is.null(invR)) {
        vif_table$vif <- diag(invR)
      }
    }
  }

  finite_vif <- vif_table$vif[is.finite(vif_table$vif)]
  max_vif <- if (length(finite_vif)) max(finite_vif) else NA_real_
  median_vif <- if (length(finite_vif)) stats::median(finite_vif) else NA_real_

  # ---------------------------------------------------------------------------
  # Variable-specific diagnostics on the same complete epidemics
  # ---------------------------------------------------------------------------
  by_variable <- do.call(
    rbind,
    lapply(var, function(v) {
      B <- blocks[[v]][complete_rows, , drop = FALSE]
      q <- qr(B, tol = tol)
      p_v <- ncol(B)
      rank_v <- q$rank

      sd_v <- apply(B, 2L, stats::sd)
      mean_v <- colMeans(B)
      nz_v <- !is.finite(sd_v) | sd_v <= tol * pmax(1, abs(mean_v))

      B_use <- B[, !nz_v, drop = FALSE]
      kappa_v <- Inf
      max_cor_v <- NA_real_

      if (ncol(B_use) == 1L) {
        kappa_v <- 1
        max_cor_v <- 0
      } else if (ncol(B_use) > 1L) {
        B_scaled <- scale(B_use)
        sv_v <- svd(B_scaled, nu = 0L, nv = 0L)$d
        if (nrow(B_scaled) >= ncol(B_scaled) &&
            qr(B_scaled, tol = tol)$rank == ncol(B_scaled) &&
            length(sv_v) && max(sv_v) > 0 &&
            min(sv_v) > tol * max(sv_v)) {
          kappa_v <- max(sv_v) / min(sv_v)
        }
        R_v <- stats::cor(B_use)
        vals <- abs(R_v[upper.tri(R_v)])
        if (length(vals)) max_cor_v <- max(vals, na.rm = TRUE)
      }

      exposure_values <- data[[v]]
      finite_exposure <- exposure_values[is.finite(exposure_values)]

      data.frame(
        variable = v,
        n_unique_exposure = length(unique(finite_exposure)),
        exposure_missing_n = sum(!is.finite(exposure_values)),
        exposure_missing_percent = 100 * mean(!is.finite(exposure_values)),
        basis_columns = p_v,
        rank = rank_v,
        full_rank = rank_v == p_v,
        rank_ratio = if (p_v > 0L) rank_v / p_v else NA_real_,
        condition_number_scaled = kappa_v,
        max_abs_within_basis_correlation = max_cor_v,
        near_zero_variance_columns = sum(nz_v),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(by_variable) <- NULL

  # ---------------------------------------------------------------------------
  # Pairwise raw-exposure correlations
  # ---------------------------------------------------------------------------
  pairwise_exposure_correlation <- data.frame(
    variable_1 = character(0),
    variable_2 = character(0),
    n_complete = integer(0),
    pearson_correlation = numeric(0),
    spearman_correlation = numeric(0),
    stringsAsFactors = FALSE
  )

  if (length(var) >= 2L) {
    pairs <- utils::combn(var, 2L, simplify = FALSE)
    pairwise_exposure_correlation <- do.call(
      rbind,
      lapply(pairs, function(pair) {
        x <- data[[pair[1L]]]
        y <- data[[pair[2L]]]
        ok <- is.finite(x) & is.finite(y)
        n_ok <- sum(ok)

        pearson <- if (n_ok >= 3L && stats::sd(x[ok]) > 0 && stats::sd(y[ok]) > 0) {
          stats::cor(x[ok], y[ok], method = "pearson")
        } else {
          NA_real_
        }

        spearman <- if (n_ok >= 3L) {
          suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
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

  # ---------------------------------------------------------------------------
  # Pairwise correlations between different cross-basis variable blocks
  # ---------------------------------------------------------------------------
  pairwise_crossbasis_correlation <- data.frame(
    variable_1 = character(0),
    variable_2 = character(0),
    max_abs_correlation = numeric(0),
    mean_abs_correlation = numeric(0),
    column_1 = character(0),
    column_2 = character(0),
    stringsAsFactors = FALSE
  )

  if (length(var) >= 2L) {
    pairs <- utils::combn(var, 2L, simplify = FALSE)
    pairwise_crossbasis_correlation <- do.call(
      rbind,
      lapply(pairs, function(pair) {
        A <- blocks[[pair[1L]]][complete_rows, , drop = FALSE]
        B <- blocks[[pair[2L]]][complete_rows, , drop = FALSE]

        sd_A <- apply(A, 2L, stats::sd)
        sd_B <- apply(B, 2L, stats::sd)
        A <- A[, is.finite(sd_A) & sd_A > tol, drop = FALSE]
        B <- B[, is.finite(sd_B) & sd_B > tol, drop = FALSE]

        if (!ncol(A) || !ncol(B)) {
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

        C <- stats::cor(A, B)
        abs_C <- abs(C)
        max_index <- which(abs_C == max(abs_C, na.rm = TRUE), arr.ind = TRUE)[1L, ]

        data.frame(
          variable_1 = pair[1L],
          variable_2 = pair[2L],
          max_abs_correlation = max(abs_C, na.rm = TRUE),
          mean_abs_correlation = mean(abs_C, na.rm = TRUE),
          column_1 = rownames(C)[max_index[1L]],
          column_2 = colnames(C)[max_index[2L]],
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(pairwise_crossbasis_correlation) <- NULL
  }

  # ---------------------------------------------------------------------------
  # Overall status and recommendations
  # ---------------------------------------------------------------------------
  severe_condition <- !is.finite(condition_number) ||
    condition_number >= condition_severe
  warning_condition <- is.finite(condition_number) &&
    condition_number >= condition_warn &&
    condition_number < condition_severe

  high_crossbasis_pair <- nrow(pairwise_crossbasis_correlation) > 0L &&
    any(
      pairwise_crossbasis_correlation$max_abs_correlation >= corr_threshold,
      na.rm = TRUE
    )

  high_raw_pair <- nrow(pairwise_exposure_correlation) > 0L &&
    any(
      abs(pairwise_exposure_correlation$pearson_correlation) >= corr_threshold,
      na.rm = TRUE
    )

  high_vif <- is.finite(max_vif) && max_vif >= vif_threshold
  low_n_per_parameter <- n_per_parameter < n_per_parameter_warn

  identifiable <- full_rank
  numerically_stable <- full_rank && !severe_condition &&
    !length(nzv_columns) && n_complete > p_full

  status <- if (!identifiable || !numerically_stable) {
    "problem"
  } else if (warning_condition || high_crossbasis_pair || high_raw_pair ||
             high_vif || low_n_per_parameter || n_excluded > 0L) {
    "warning"
  } else {
    "ok"
  }

  recommendations <- character(0)

  if (!full_rank) {
    recommendations <- c(
      recommendations,
      paste0(
        "The combined epidemic-level design is rank-deficient (rank = ",
        rank_full, ", columns = ", p_full, "). Inspect the reported dependent ",
        "columns and simplify the exposure-lag basis. Reducing df_var and/or ",
        "df_lag is usually preferable before shortening max_lag; change max_lag ",
        "only when scientifically justified."
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
        "The scaled condition number is ",
        if (is.finite(condition_number)) format(condition_number, digits = 4) else "infinite",
        ", indicating severe numerical instability. Consider reducing basis ",
        "complexity and checking whether exposure variables provide redundant ",
        "cross-basis information."
      )
    )
  } else if (warning_condition) {
    recommendations <- c(
      recommendations,
      paste0(
        "The scaled condition number is ", format(condition_number, digits = 4),
        ". This suggests moderate-to-strong numerical collinearity; inspect the ",
        "pairwise cross-basis diagnostics before model fitting."
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
        "Strong correlation was detected between cross-basis blocks: ",
        flagged_text,
        ". Consider whether these exposures are providing redundant temporal ",
        "information, and compare simpler or alternative variable sets."
      )
    )
  }

  if (high_raw_pair) {
    flagged <- pairwise_exposure_correlation[
      is.finite(pairwise_exposure_correlation$pearson_correlation) &
        abs(pairwise_exposure_correlation$pearson_correlation) >= corr_threshold,
      ,
      drop = FALSE
    ]
    flagged_text <- paste(
      paste0(
        flagged$variable_1, " vs ", flagged$variable_2,
        " (r = ", format(flagged$pearson_correlation, digits = 3), ")"
      ),
      collapse = "; "
    )
    recommendations <- c(
      recommendations,
      paste0(
        "High same-time correlation was detected among raw exposures: ",
        flagged_text,
        ". Use this as descriptive evidence only; the cross-basis diagnostics ",
        "are more directly relevant to the fitted DLNM design."
      )
    )
  }

  if (high_vif) {
    recommendations <- c(
      recommendations,
      paste0(
        "At least one cross-basis column has VIF >= ", vif_threshold,
        ". Treat this as a supplementary signal because spline-basis columns ",
        "can be correlated by construction; prioritize rank and condition-number ",
        "diagnostics."
      )
    )
  }

  if (low_n_per_parameter) {
    recommendations <- c(
      recommendations,
      paste0(
        "There are only ", format(n_per_parameter, digits = 3),
        " complete epidemics per design column. This may indicate an ",
        "overparameterized exposure-lag specification; consider reducing df_var ",
        "and/or df_lag or increasing the number of independent epidemics."
      )
    )
  }

  if (n_excluded > 0L) {
    recommendations <- c(
      recommendations,
      paste0(
        n_excluded, " of ", n_total,
        " epidemics were excluded from the combined diagnostic because their ",
        "cross-basis design contained missing values. Inspect missing exposure ",
        "values within the relevant lag history."
      )
    )
  }

  low_unique <- if (fun_var == "lin") {
    rep(FALSE, nrow(by_variable))
  } else {
    by_variable$n_unique_exposure <= df_var
  }
  if (any(low_unique)) {
    recommendations <- c(
      recommendations,
      paste0(
        "Limited exposure support relative to df_var was detected for: ",
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

  overall <- data.frame(
    n_epidemics_total = n_total,
    n_epidemics_complete = n_complete,
    n_epidemics_excluded = n_excluded,
    n_exposures = length(var),
    crossbasis_columns = ncol(X_cb),
    design_columns = p_full,
    rank = rank_full,
    full_rank = full_rank,
    rank_ratio = if (p_full > 0L) rank_full / p_full else NA_real_,
    condition_number_scaled = condition_number,
    min_singular_value_scaled = min_singular_value,
    max_singular_value_scaled = max_singular_value,
    max_abs_crossbasis_correlation = max_abs_cb_correlation,
    median_abs_crossbasis_correlation = median_abs_cb_correlation,
    max_vif = max_vif,
    median_vif = median_vif,
    near_zero_variance_columns = length(nzv_columns),
    epidemics_per_design_column = n_per_parameter,
    time_step = time_step,
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
    dependent_columns = dependent_columns,
    near_zero_variance_columns = nzv_columns,
    complete_epi_id = complete_ids,
    recommendations = unique(recommendations),
    settings = list(
      variables = var,
      max_lag = max_lag,
      df_var = df_var,
      df_lag = df_lag,
      fun_var = fun_var,
      fun_lag = fun_lag,
      include_intercept = include_intercept,
      corr_threshold = corr_threshold,
      condition_warn = condition_warn,
      condition_severe = condition_severe,
      vif_threshold = vif_threshold,
      n_per_parameter_warn = n_per_parameter_warn,
      tol = tol
    )
  )

  if (keep_design) {
    out$design_matrix <- X_rank
  }

  class(out) <- c("epiexposure_identifiability", "list")
  out
}

#' Print DLNM identifiability diagnostics
#'
#' @param x Object returned by `check_identifiability()`.
#' @param ... Unused.
#' @export
print.epiexposure_identifiability <- function(x, ...) {
  cat("EpiExposure DLNM identifiability diagnostics\n")
  cat("Status: ", toupper(x$status), "\n", sep = "")
  cat("Identifiable (full rank): ", x$identifiable, "\n", sep = "")
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
    cat("\nRaw exposure correlation\n")
    print(x$pairwise_exposure_correlation, row.names = FALSE)
  }

  if (length(x$dependent_columns)) {
    cat("\nDependent columns\n")
    cat(paste0("- ", x$dependent_columns), sep = "\n")
    cat("\n")
  }

  cat("\nRecommendations\n")
  cat(paste0("- ", x$recommendations), sep = "\n")
  cat("\n")

  invisible(x)
}
