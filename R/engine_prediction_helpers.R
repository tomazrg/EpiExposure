# ============================================================================
# EpiExposure - internal engine/prediction helpers
# ============================================================================
#
# This file centralizes the prediction contract used by EpiExposure v1.
# It is intentionally internal: none of the functions below are exported.
#
# Core contract
# -------------
# 1. Prediction target: expected response, not a newly simulated observation.
# 2. Prediction level: population level; fitted group-specific random effects
#    are excluded from downstream predictions.
# 3. uncertainty = FALSE (downstream): central parameter estimate.
#    - frequentist engines: fitted fixed-effect estimates;
#    - Bayesian engines: posterior mean of population/fixed effects.
# 4. uncertainty = TRUE (downstream): parameter draws are propagated
#    draw-by-draw through the link transformation, then summarized by median,
#    SD, and empirical quantiles.
# 5. No residual/process noise is added by these helpers.
# 6. Family, link, engine, basis specification, and coefficient mapping come
#    from fitted EpiExposure metadata. No family/link guessing is permitted.
# 7. Cross-basis coefficients are aligned by names. For bdlnm, whose native
#    coefficient names follow the stored basis names, the mapping is derived
#    explicitly from the stored crossbasis column names; it is never aligned
#    merely by position.
# 8. All fitted exposures must share one common non-negative integer `max_lag`.
#    Every prediction history must contain exactly `max_lag + 1` observations
#    for every fitted exposure. Histories are never truncated, padded, or
#    silently aligned across different lag windows.
#
# Expected use after this file is added to R/:
#   predict_outcomes()   -> .build_newdata_basis()
#                        -> .predict_point_population() or
#                           .predict_draws_population()
#                        -> .summarise_prediction_draws()
#
#   find_bestfit()       -> same point-prediction helper
#   summarise_effects()  -> same parameter-draw helper
#   simulate_scenarios() -> same basis/prediction helpers
#
# ============================================================================


# ----------------------------------------------------------------------------
# Small generic helpers
# ----------------------------------------------------------------------------

.epix_stop <- function(...) {
  stop(..., call. = FALSE)
}

.epix_is_scalar_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.epix_is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.epix_validate_n_samples <- function(n_samples) {
  if (!.epix_is_scalar_number(n_samples) || n_samples < 1 ||
      n_samples != as.integer(n_samples)) {
    .epix_stop("`n_samples` must be one positive integer.")
  }
  as.integer(n_samples)
}

.epix_validate_probs <- function(probs) {
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs <= 0 | probs >= 1) ||
      probs[1] >= probs[2]) {
    .epix_stop(
      "`probs` must contain two finite probabilities with ",
      "0 < probs[1] < probs[2] < 1."
    )
  }
  as.numeric(probs)
}

.epix_sort_cb_cols <- function(columns) {
  if (!length(columns)) return(character(0))

  index <- suppressWarnings(
    as.integer(sub("^.*_([0-9]+)$", "\\1", columns))
  )

  if (anyNA(index)) {
    .epix_stop(
      "Cross-basis columns must end in a numeric index. Problematic column(s): ",
      paste(columns[is.na(index)], collapse = ", "), "."
    )
  }

  columns[order(index)]
}

.epix_cb_cols_for_var <- function(cb_cols, variable) {
  prefix <- paste0("cb_", variable, "_")
  .epix_sort_cb_cols(cb_cols[startsWith(cb_cols, prefix)])
}

.epix_link_object <- function(link_name) {
  if (!.epix_is_scalar_string(link_name)) {
    .epix_stop("Stored `epiexposure_link` must be one non-empty character value.")
  }

  out <- tryCatch(
    stats::make.link(link_name),
    error = function(e) NULL
  )

  if (is.null(out) || !is.function(out$linkinv) || !is.function(out$linkfun)) {
    .epix_stop(
      "Could not construct the stored model link '", link_name,
      "'. EpiExposure does not infer or replace missing links."
    )
  }

  out
}

.epix_as_numeric_matrix <- function(x, object_name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (!length(x) || any(!is.finite(x))) {
    .epix_stop("`", object_name, "` must contain only finite numeric values.")
  }

  x
}


# ----------------------------------------------------------------------------
# Strict fitted-model metadata contract
# ----------------------------------------------------------------------------

#' Read and validate EpiExposure model metadata
#'
#' Internal helper. It validates the fitted-model contract established by
#' `fit_epidlnm()` and returns a normalized metadata list. Missing family, link,
#' engine, basis, or prediction-contract information is never reconstructed
#' from unrelated model fields.
#'
#' For `bdlnm`, a canonical `cb_<variable>_<index>` naming layer can be created
#' deterministically from the stored crossbasis objects if the engine-native fit
#' did not require `cb_*` columns in its input data. This is not coefficient
#' alignment by position: later mapping still uses the explicit native
#' crossbasis column names stored in each basis object.
#'
#' @noRd
.get_epiexposure_metadata <- function(fit) {

  if (is.null(fit)) {
    .epix_stop("`fit` cannot be NULL.")
  }

  get_attr <- function(name, allow_null = FALSE) {
    value <- attr(fit, name, exact = TRUE)
    if (is.null(value) && !allow_null) {
      .epix_stop(
        "The fitted model is missing required EpiExposure metadata `", name,
        "`. Refit the model with the current `fit_epidlnm()` implementation."
      )
    }
    value
  }

  engine <- get_attr("epiexposure_engine")
  family_name <- get_attr("epiexposure_family_name")
  link <- get_attr("epiexposure_link")
  cb_cols <- get_attr("epiexposure_cb_cols", allow_null = TRUE)
  vars <- get_attr("epiexposure_vars")
  data_template <- get_attr("epiexposure_data_template")
  id_col <- get_attr("epiexposure_id_col", allow_null = TRUE)
  random_fit_col <- get_attr(
    "epiexposure_random_effect_fit_col",
    allow_null = TRUE
  )
  random_structure <- get_attr("epiexposure_random_structure")
  spec <- get_attr("epiexposure_spec")
  basis_objects <- get_attr("epiexposure_basis_objects", allow_null = TRUE)
  prediction_level <- get_attr("epiexposure_prediction_level")
  prediction_estimand <- get_attr("epiexposure_prediction_estimand")
  point_contract <- get_attr("epiexposure_point_prediction_contract")
  uncertainty_contract <- get_attr("epiexposure_uncertainty_contract")
  stored_max_lag <- get_attr("epiexposure_max_lag")
  stored_history_length <- get_attr("epiexposure_history_length")
  history_contract <- get_attr("epiexposure_history_contract")

  supported_engines <- c(
    "glm", "glmmTMB", "gam", "gamm", "gls",
    "spamm", "brms", "inla", "bdlnm"
  )
  supported_families <- c(
    "beta", "binomial", "poisson", "gamma", "gaussian",
    "negative_binomial"
  )

  if (!.epix_is_scalar_string(engine) || !engine %in% supported_engines) {
    .epix_stop(
      "Invalid stored `epiexposure_engine`. Expected one of: ",
      paste(supported_engines, collapse = ", "), "."
    )
  }

  if (!.epix_is_scalar_string(family_name) ||
      !family_name %in% supported_families) {
    .epix_stop(
      "Invalid stored `epiexposure_family_name`. EpiExposure v1 supports: ",
      paste(supported_families, collapse = ", "), "."
    )
  }

  if (!.epix_is_scalar_string(link)) {
    .epix_stop("Invalid stored `epiexposure_link`.")
  }
  invisible(.epix_link_object(link))

  if (!is.character(vars) || !length(vars) || anyNA(vars) ||
      any(!nzchar(vars)) || anyDuplicated(vars)) {
    .epix_stop("Stored `epiexposure_vars` must contain unique exposure names.")
  }

  if (!is.data.frame(data_template) || nrow(data_template) < 1L) {
    .epix_stop(
      "Stored `epiexposure_data_template` must be a non-empty data.frame."
    )
  }

  if (!is.list(spec) || is.null(names(spec)) || anyNA(names(spec)) ||
      any(!nzchar(names(spec))) || anyDuplicated(names(spec))) {
    .epix_stop("Stored `epiexposure_spec` must be a named list with unique names.")
  }

  if (!setequal(names(spec), vars)) {
    .epix_stop(
      "Stored `epiexposure_spec` names must match `epiexposure_vars` exactly."
    )
  }
  spec <- spec[vars]

  fitted_max_lags <- vapply(
    vars,
    function(variable) {
      current <- spec[[variable]]

      if (!is.list(current) || is.null(current$max_lag) ||
          is.null(current$argvar) || is.null(current$arglag)) {
        .epix_stop(
          "Incomplete stored exposure specification for variable '", variable,
          "'. Required elements are `max_lag`, `argvar`, and `arglag`."
        )
      }

      max_lag_value <- current$max_lag

      if (!is.numeric(max_lag_value) ||
          length(max_lag_value) != 1L ||
          is.na(max_lag_value) ||
          !is.finite(max_lag_value) ||
          max_lag_value < 0 ||
          max_lag_value != as.integer(max_lag_value)) {
        .epix_stop(
          "Stored `max_lag` for variable '", variable,
          "' must be one non-negative finite integer."
        )
      }

      if (!is.list(current$argvar) || !is.list(current$arglag)) {
        .epix_stop(
          "Stored `argvar` and `arglag` must be lists for variable '",
          variable, "'."
        )
      }

      as.integer(max_lag_value)
    },
    integer(1)
  )
  names(fitted_max_lags) <- vars

  if (length(unique(fitted_max_lags)) != 1L) {
    .epix_stop(
      "All fitted exposure variables must use the same `max_lag` under the ",
      "EpiExposure exact-history contract. Stored values are: ",
      paste(
        paste0(names(fitted_max_lags), "=", fitted_max_lags),
        collapse = ", "
      ),
      ". Refit the model with a common lag window."
    )
  }

  common_max_lag <- unname(fitted_max_lags[1L])
  expected_history_length <- common_max_lag + 1L

  if (!.epix_is_scalar_number(stored_max_lag) ||
      stored_max_lag < 0 ||
      stored_max_lag != as.integer(stored_max_lag)) {
    .epix_stop(
      "Stored `epiexposure_max_lag` must be one non-negative integer."
    )
  }
  stored_max_lag <- as.integer(stored_max_lag)

  if (!identical(stored_max_lag, common_max_lag)) {
    .epix_stop(
      "Stored `epiexposure_max_lag` (", stored_max_lag,
      ") does not match the common `epiexposure_spec` max_lag (",
      common_max_lag, ")."
    )
  }

  if (!.epix_is_scalar_number(stored_history_length) ||
      stored_history_length < 1 ||
      stored_history_length != as.integer(stored_history_length)) {
    .epix_stop(
      "Stored `epiexposure_history_length` must be one positive integer."
    )
  }
  stored_history_length <- as.integer(stored_history_length)

  if (!identical(stored_history_length, expected_history_length)) {
    .epix_stop(
      "Stored `epiexposure_history_length` must equal `max_lag + 1`. ",
      "Expected ", expected_history_length, " but found ",
      stored_history_length, "."
    )
  }

  if (!identical(
    history_contract,
    "all_fitted_exposures_same_exact_max_lag_plus_one"
  )) {
    .epix_stop(
      "Invalid stored `epiexposure_history_contract`. Expected ",
      "'all_fitted_exposures_same_exact_max_lag_plus_one'."
    )
  }

  for (variable in vars) {
    spec[[variable]]$max_lag <- common_max_lag
  }

  if (!is.null(basis_objects)) {
    if (!is.list(basis_objects) || is.null(names(basis_objects)) ||
        anyNA(names(basis_objects)) || any(!nzchar(names(basis_objects))) ||
        anyDuplicated(names(basis_objects))) {
      .epix_stop(
        "Stored `epiexposure_basis_objects` must be NULL or a named list ",
        "with unique names."
      )
    }

    if (!setequal(names(basis_objects), vars)) {
      .epix_stop(
        "Stored basis-object names must match `epiexposure_vars` exactly."
      )
    }
    basis_objects <- basis_objects[vars]

    for (variable in vars) {
      basis_object <- basis_objects[[variable]]

      if (!inherits(basis_object, "crossbasis")) {
        .epix_stop(
          "Stored basis for variable '", variable,
          "' must inherit from `crossbasis` in EpiExposure v1."
        )
      }

      basis_lag <- attr(
        basis_object,
        "lag",
        exact = TRUE
      )

      if (is.null(basis_lag) ||
          !is.numeric(basis_lag) ||
          !length(basis_lag) ||
          anyNA(basis_lag) ||
          any(!is.finite(basis_lag))) {
        .epix_stop(
          "Stored basis for variable '", variable,
          "' has invalid lag metadata."
        )
      }

      basis_max_lag <- as.integer(max(basis_lag))

      if (!identical(basis_max_lag, common_max_lag)) {
        .epix_stop(
          "Stored basis/spec lag mismatch for variable '", variable,
          "': basis max lag = ", basis_max_lag,
          ", common fitted max lag = ", common_max_lag, "."
        )
      }
    }
  }

  if (is.null(cb_cols)) cb_cols <- character(0)
  if (!is.character(cb_cols) || anyNA(cb_cols) || anyDuplicated(cb_cols)) {
    .epix_stop("Stored `epiexposure_cb_cols` must be a unique character vector.")
  }

  cb_cols_source <- "metadata"

  # bdlnm can be fit from basis objects without an explicit cb_* data matrix.
  # In that one engine-specific case, create the package's canonical design
  # names from the *stored crossbasis objects*. The native coefficient mapping
  # remains name-based and is validated separately below.
  if (!length(cb_cols)) {
    if (!identical(engine, "bdlnm") || is.null(basis_objects)) {
      .epix_stop(
        "No canonical cross-basis columns are stored in `epiexposure_cb_cols`."
      )
    }

    cb_cols <- unlist(
      lapply(vars, function(variable) {
        paste0("cb_", variable, "_", seq_len(ncol(basis_objects[[variable]])))
      }),
      use.names = FALSE
    )
    cb_cols_source <- "stored_basis_objects"
  }

  for (variable in vars) {
    cols <- .epix_cb_cols_for_var(cb_cols, variable)
    if (!length(cols)) {
      .epix_stop(
        "No canonical cross-basis columns were found for variable '",
        variable, "'."
      )
    }

    if (!is.null(basis_objects) &&
        length(cols) != ncol(basis_objects[[variable]])) {
      .epix_stop(
        "Cross-basis dimension mismatch for variable '", variable,
        "': metadata contains ", length(cols), " column(s), but the stored ",
        "basis contains ", ncol(basis_objects[[variable]]), "."
      )
    }
  }

  if (!is.null(id_col) && !.epix_is_scalar_string(id_col)) {
    .epix_stop("Stored `epiexposure_id_col` must be NULL or one column name.")
  }

  if (!is.null(random_fit_col) && !.epix_is_scalar_string(random_fit_col)) {
    .epix_stop(
      "Stored `epiexposure_random_effect_fit_col` must be NULL or one column name."
    )
  }

  if (!.epix_is_scalar_string(random_structure) ||
      !random_structure %in% c("none", "random_intercept")) {
    .epix_stop(
      "Stored `epiexposure_random_structure` must be 'none' or ",
      "'random_intercept'."
    )
  }

  if (identical(random_structure, "random_intercept") && is.null(id_col)) {
    .epix_stop(
      "Random-intercept metadata are inconsistent: `epiexposure_id_col` is NULL."
    )
  }

  if (!identical(prediction_level, "population")) {
    .epix_stop(
      "This EpiExposure v1 prediction layer requires ",
      "`epiexposure_prediction_level = 'population'`."
    )
  }

  if (!identical(prediction_estimand, "expected_response")) {
    .epix_stop(
      "This EpiExposure v1 prediction layer requires ",
      "`epiexposure_prediction_estimand = 'expected_response'`."
    )
  }

  if (!identical(point_contract, "central_expected_response")) {
    .epix_stop(
      "Invalid stored point-prediction contract. Expected ",
      "'central_expected_response'."
    )
  }

  if (!identical(uncertainty_contract, "draw_by_draw_median_quantiles")) {
    .epix_stop(
      "Invalid stored uncertainty contract. Expected ",
      "'draw_by_draw_median_quantiles'."
    )
  }

  list(
    engine = engine,
    family_name = family_name,
    link = link,
    cb_cols = cb_cols,
    cb_cols_source = cb_cols_source,
    vars = vars,
    data_template = data_template[1, , drop = FALSE],
    id_col = id_col,
    random_effect_fit_col = random_fit_col,
    random_structure = random_structure,
    spec = spec,
    basis_objects = basis_objects,
    max_lag = common_max_lag,
    history_length = expected_history_length,
    history_contract = history_contract,
    prediction_level = prediction_level,
    prediction_estimand = prediction_estimand,
    point_prediction_contract = point_contract,
    uncertainty_contract = uncertainty_contract
  )
}


# ----------------------------------------------------------------------------
# Basis reconstruction for new observed data or explicit exposure profiles
# ----------------------------------------------------------------------------

.epix_basis_definition <- function(metadata, variable) {
  spec <- metadata$spec[[variable]]
  max_lag <- as.integer(spec$max_lag)

  if (!identical(max_lag, metadata$max_lag)) {
    .epix_stop(
      "Stored exposure specification for variable '", variable,
      "' does not match the common fitted `max_lag`."
    )
  }
  argvar <- spec$argvar
  arglag <- spec$arglag

  # `define_exposures()` stores the effective fitted exposure range in the
  # specification. Use it even when the complete training cross-basis object
  # was not retained by `fit_epidlnm()`. This keeps extrapolation checks active
  # without inferring the fitted range from prediction data.
  exposure_range <- NULL
  if (!is.null(spec$range) &&
      is.numeric(spec$range) &&
      length(spec$range) >= 2L &&
      !anyNA(spec$range) &&
      all(is.finite(spec$range))) {
    exposure_range <- range(as.numeric(spec$range))
  }

  expected_ncol <- length(.epix_cb_cols_for_var(metadata$cb_cols, variable))

  if (is.list(metadata$basis_objects) &&
      !is.null(names(metadata$basis_objects)) &&
      variable %in% names(metadata$basis_objects)) {
    template <- metadata$basis_objects[[variable]]

    if (!inherits(template, "crossbasis")) {
      .epix_stop(
        "Stored basis object for variable '", variable,
        "' must inherit from `crossbasis`."
      )
    }

    template_argvar <- attr(template, "argvar", exact = TRUE)
    template_arglag <- attr(template, "arglag", exact = TRUE)
    template_lag <- attr(template, "lag", exact = TRUE)
    template_range <- attr(template, "range", exact = TRUE)

    if (is.list(template_argvar)) argvar <- template_argvar
    if (is.list(template_arglag)) arglag <- template_arglag

    if (!is.null(template_lag)) {
      template_lag <- as.numeric(template_lag)
      if (length(template_lag) >= 1L && all(is.finite(template_lag))) {
        template_max_lag <- as.integer(max(template_lag))
        if (!identical(template_max_lag, max_lag)) {
          .epix_stop(
            "Stored basis/spec lag mismatch for variable '", variable,
            "': basis max lag = ", template_max_lag,
            ", spec max lag = ", max_lag, "."
          )
        }
      }
    }

    if (!is.null(template_range) && length(template_range) >= 2L &&
        all(is.finite(template_range))) {
      exposure_range <- range(as.numeric(template_range))
    }

    if (ncol(template) != expected_ncol) {
      .epix_stop(
        "Stored cross-basis dimension mismatch for variable '", variable, "'."
      )
    }
  }

  list(
    max_lag = max_lag,
    argvar = argvar,
    arglag = arglag,
    exposure_range = exposure_range,
    expected_ncol = expected_ncol
  )
}

.epix_handle_extrapolation <- function(values, exposure_range, variable,
                                       extrapolation) {
  if (is.null(exposure_range) || !length(values)) return(invisible(NULL))

  outside <- values < exposure_range[1] | values > exposure_range[2]
  if (!any(outside)) return(invisible(NULL))

  msg <- paste0(
    "Prediction values for variable '", variable,
    "' extend outside the exposure range used to construct the stored basis [",
    format(exposure_range[1]), ", ", format(exposure_range[2]), "]."
  )

  if (identical(extrapolation, "error")) {
    .epix_stop(msg)
  }

  if (identical(extrapolation, "warn")) {
    warning(msg, call. = FALSE)
  }

  invisible(NULL)
}

.epix_build_cb_row <- function(history, variable, metadata,
                               extrapolation = c("warn", "error", "allow")) {
  extrapolation <- match.arg(extrapolation)

  if (!requireNamespace("dlnm", quietly = TRUE)) {
    .epix_stop("Package 'dlnm' is required to construct prediction cross-bases.")
  }

  definition <- .epix_basis_definition(metadata, variable)

  if (!is.numeric(history) || !length(history) ||
      anyNA(history) || any(!is.finite(history))) {
    .epix_stop(
      "Exposure history for variable '", variable,
      "' must contain only finite numeric values."
    )
  }

  expected_length <- metadata$history_length

  if (length(history) != expected_length) {
    .epix_stop(
      "Exposure history for variable '", variable,
      "' must contain exactly ", expected_length,
      " observation(s) for fitted max_lag = ", metadata$max_lag,
      "; received ", length(history), ". Histories are not truncated, ",
      "padded, or silently realigned."
    )
  }

  .epix_handle_extrapolation(
    history,
    definition$exposure_range,
    variable,
    extrapolation
  )

  cb <- dlnm::crossbasis(
    history,
    lag = definition$max_lag,
    argvar = definition$argvar,
    arglag = definition$arglag
  )

  if (ncol(cb) != definition$expected_ncol) {
    .epix_stop(
      "Reconstructed cross-basis for variable '", variable,
      "' has ", ncol(cb), " column(s), but the fitted model expects ",
      definition$expected_ncol, ". The stored basis specification cannot be ",
      "safely applied to these prediction data."
    )
  }

  row <- as.numeric(cb[nrow(cb), , drop = TRUE])

  if (length(row) != definition$expected_ncol || any(!is.finite(row))) {
    .epix_stop(
      "Could not obtain a finite final cross-basis row for variable '",
      variable, "'."
    )
  }

  names(row) <- .epix_cb_cols_for_var(metadata$cb_cols, variable)
  row
}

.epix_validate_regular_time <- function(data, group_index, time_col) {
  group_ids <- unique(group_index)
  steps <- numeric(0)

  for (g in group_ids) {
    idx <- which(group_index == g)
    tt <- data[[time_col]][idx]

    if (!is.numeric(tt) || anyNA(tt) || any(!is.finite(tt))) {
      .epix_stop(
        "Prediction time column '", time_col,
        "' must contain only finite numeric values."
      )
    }

    tt <- sort(tt)

    if (anyDuplicated(tt)) {
      .epix_stop(
        "Duplicated time values were found within prediction group '",
        as.character(g), "'."
      )
    }

    if (length(tt) > 1L) {
      d <- diff(tt)
      step <- d[1]
      tol <- sqrt(.Machine$double.eps) * max(1, abs(step))

      if (step <= 0 || any(abs(d - step) > tol)) {
        .epix_stop(
          "Prediction histories must be complete and equally spaced within ",
          "each group. Irregular spacing was detected for group '",
          as.character(g), "'."
        )
      }

      steps <- c(steps, step)
    }
  }

  if (length(steps) > 1L) {
    ref_step <- steps[1]
    tol <- sqrt(.Machine$double.eps) * max(1, abs(ref_step))

    if (any(abs(steps - ref_step) > tol)) {
      .epix_stop(
        "All prediction groups must use the same temporal spacing because ",
        "DLNM lag indices must represent the same time unit across groups."
      )
    }
  }

  invisible(TRUE)
}

.epix_prediction_template <- function(metadata, n) {
  template <- metadata$data_template[1, , drop = FALSE]
  out <- template[rep(1L, n), , drop = FALSE]
  rownames(out) <- NULL

  # Ensure canonical cb_* columns exist even for an engine (notably bdlnm)
  # whose native fit may not have required them in the original data frame.
  for (column in metadata$cb_cols) {
    if (!column %in% names(out)) out[[column]] <- NA_real_
  }

  out
}

#' Build a prediction design using the fitted cross-basis definition
#'
#' Internal helper for both observed new data and explicit exposure profiles.
#'
#' Exactly one of `newdata` or `profiles` must be supplied.
#'
#' `newdata` is a long-format exposure history. One epidemic-level design row is
#' returned per `group`, using the final time point after chronological sorting.
#' The time series must be finite, unique, complete, equally spaced, and contain
#' exactly the common fitted `max_lag + 1` observations per group. If
#' `group = NULL`, the entire data frame is treated as one exposure history and
#' must have exactly that same length.
#'
#' `profiles` is a named list containing one chronological exposure-history
#' vector for every fitted exposure. Vectors run from the earliest observation
#' to the most recent observation. Every fitted exposure must use exactly the
#' same history length, and that length must equal the common fitted
#' `max_lag + 1`. Longer histories are not truncated and shorter histories are
#' not padded.
#'
#' The fitted basis is never re-estimated from prediction data. Stored crossbasis
#' `argvar`/`arglag` attributes are preferred; otherwise the strict stored
#' `epiexposure_spec` is used. The reconstructed basis dimension is checked
#' against the fitted canonical `cb_*` columns.
#'
#' @param fit Fitted EpiExposure model.
#' @param newdata Optional long-format prediction data.
#' @param profiles Optional named list of complete exposure histories.
#' @param group Grouping column in `newdata`, or `NULL` for one history.
#' @param time Chronological numeric time column in `newdata`.
#' @param extrapolation Behavior when exposure values extend beyond the range
#'   stored in a crossbasis object: `"warn"`, `"error"`, or `"allow"`.
#'
#' @return Epidemic-level prediction data frame with attributes
#'   `epiexposure_prediction_keys` and `epiexposure_prediction_source`.
#'
#' @noRd
.build_newdata_basis <- function(
    fit,
    newdata = NULL,
    profiles = NULL,
    group = "epi_id",
    time = "time",
    extrapolation = c("warn", "error", "allow")
) {
  extrapolation <- match.arg(extrapolation)
  metadata <- .get_epiexposure_metadata(fit)

  has_newdata <- !is.null(newdata)
  has_profiles <- !is.null(profiles)

  if (identical(has_newdata, has_profiles)) {
    .epix_stop("Supply exactly one of `newdata` or `profiles`.")
  }

  # --------------------------------------------------------------------------
  # Explicit profile
  # --------------------------------------------------------------------------
  if (has_profiles) {
    if (!is.list(profiles) || is.null(names(profiles)) ||
        anyNA(names(profiles)) || any(!nzchar(names(profiles))) ||
        anyDuplicated(names(profiles))) {
      .epix_stop("`profiles` must be a named list with unique exposure names.")
    }

    missing_vars <- setdiff(metadata$vars, names(profiles))
    extra_vars <- setdiff(names(profiles), metadata$vars)

    if (length(missing_vars) || length(extra_vars)) {
      details <- c(
        if (length(missing_vars)) {
          paste0("missing: ", paste(missing_vars, collapse = ", "))
        },
        if (length(extra_vars)) {
          paste0("unexpected: ", paste(extra_vars, collapse = ", "))
        }
      )
      .epix_stop(
        "`profiles` names must match the fitted exposures exactly (",
        paste(details, collapse = "; "), ")."
      )
    }

    profile_lengths <- vapply(
      profiles[metadata$vars],
      length,
      integer(1)
    )

    if (length(unique(profile_lengths)) != 1L) {
      .epix_stop(
        "All fitted exposure profiles must contain the same number of time ",
        "observations. Received: ",
        paste(
          paste0(names(profile_lengths), "=", profile_lengths),
          collapse = ", "
        ),
        "."
      )
    }

    if (any(profile_lengths != metadata$history_length)) {
      .epix_stop(
        "Every fitted exposure profile must contain exactly ",
        metadata$history_length, " observations for fitted max_lag = ",
        metadata$max_lag, ". Received: ",
        paste(
          paste0(names(profile_lengths), "=", profile_lengths),
          collapse = ", "
        ),
        ". Longer histories are not truncated and shorter histories are not ",
        "padded."
      )
    }

    design <- .epix_prediction_template(metadata, 1L)

    for (variable in metadata$vars) {
      row <- .epix_build_cb_row(
        history = profiles[[variable]],
        variable = variable,
        metadata = metadata,
        extrapolation = extrapolation
      )
      design[names(row)] <- as.list(row)
    }

    keys <- data.frame(profile_id = 1L, stringsAsFactors = FALSE)
    attr(design, "epiexposure_prediction_keys") <- keys
    attr(design, "epiexposure_prediction_source") <- "profiles"
    attr(design, "epiexposure_history_contract") <- metadata$history_contract
    attr(design, "epiexposure_max_lag") <- metadata$max_lag
    attr(design, "epiexposure_history_length") <- metadata$history_length
    attr(design, "epiexposure_metadata") <- metadata
    return(design)
  }

  # --------------------------------------------------------------------------
  # Long-format observed new data
  # --------------------------------------------------------------------------
  if (!is.data.frame(newdata) || !nrow(newdata)) {
    .epix_stop("`newdata` must be a non-empty data.frame.")
  }

  if (!.epix_is_scalar_string(time) || !time %in% names(newdata)) {
    .epix_stop("`time` must name a column present in `newdata`.")
  }

  missing_vars <- setdiff(metadata$vars, names(newdata))
  if (length(missing_vars)) {
    .epix_stop(
      "Prediction data are missing fitted exposure variable(s): ",
      paste(missing_vars, collapse = ", "), "."
    )
  }

  for (variable in metadata$vars) {
    if (!is.numeric(newdata[[variable]]) || anyNA(newdata[[variable]]) ||
        any(!is.finite(newdata[[variable]]))) {
      .epix_stop(
        "Prediction exposure variable '", variable,
        "' must contain only finite numeric values."
      )
    }
  }

  if (is.null(group)) {
    group_index <- rep(".epiexposure_single_history", nrow(newdata))
    group_values <- ".epiexposure_single_history"
    group_name <- NULL
  } else {
    if (!.epix_is_scalar_string(group) || !group %in% names(newdata)) {
      .epix_stop(
        "`group` must be NULL or name a grouping column present in `newdata`."
      )
    }
    if (anyNA(newdata[[group]])) {
      .epix_stop("Prediction grouping column cannot contain missing values.")
    }

    group_index <- as.character(newdata[[group]])
    group_values <- unique(group_index)
    group_name <- group
  }

  .epix_validate_regular_time(newdata, group_index, time)

  required_history_length <- metadata$history_length

  group_counts <- table(group_index)
  invalid_counts <- names(group_counts)[
    group_counts != required_history_length
  ]

  if (length(invalid_counts)) {
    examples <- paste0(
      invalid_counts,
      "=",
      as.integer(group_counts[invalid_counts])
    )

    .epix_stop(
      "Every prediction history must contain exactly ",
      required_history_length, " observations for fitted max_lag = ",
      metadata$max_lag, ". Non-matching group(s) include: ",
      paste(utils::head(examples, 5L), collapse = ", "),
      if (length(examples) > 5L) "; ..." else ".",
      " Histories are not truncated, padded, or silently realigned."
    )
  }

  design <- .epix_prediction_template(metadata, length(group_values))

  for (i in seq_along(group_values)) {
    g <- group_values[i]
    idx <- which(group_index == g)
    idx <- idx[order(newdata[[time]][idx])]

    for (variable in metadata$vars) {
      row <- .epix_build_cb_row(
        history = newdata[[variable]][idx],
        variable = variable,
        metadata = metadata,
        extrapolation = extrapolation
      )
      design[i, names(row)] <- as.list(row)
    }
  }

  if (is.null(group_name)) {
    keys <- data.frame(prediction_id = 1L, stringsAsFactors = FALSE)
  } else {
    keys <- data.frame(group_values, stringsAsFactors = FALSE)
    names(keys) <- group_name
  }

  attr(design, "epiexposure_prediction_keys") <- keys
  attr(design, "epiexposure_prediction_source") <- "newdata"
  attr(design, "epiexposure_history_contract") <- metadata$history_contract
  attr(design, "epiexposure_max_lag") <- metadata$max_lag
  attr(design, "epiexposure_history_length") <- metadata$history_length
  attr(design, "epiexposure_metadata") <- metadata
  design
}


# ----------------------------------------------------------------------------
# Canonical parameter-name mapping
# ----------------------------------------------------------------------------

.epix_expected_parameter_names <- function(metadata) {
  c("(Intercept)", metadata$cb_cols)
}

.epix_match_one_name <- function(candidates, raw_names, label) {
  candidates <- unique(candidates[nzchar(candidates)])
  matches <- intersect(candidates, raw_names)

  if (length(matches) == 1L) return(matches)

  if (length(matches) > 1L) {
    .epix_stop(
      "Ambiguous coefficient-name mapping for ", label, ": ",
      paste(matches, collapse = ", "), "."
    )
  }

  character(0)
}

.epix_bdlnm_parameter_map <- function(fit, metadata, raw_names) {
  if (is.null(metadata$basis_objects)) {
    .epix_stop(
      "bdlnm prediction requires stored `epiexposure_basis_objects` metadata."
    )
  }

  map <- setNames(character(0), character(0))

  intercept_matches <- raw_names[
    raw_names %in% c("(Intercept)", "Intercept") |
      grepl("Intercept", raw_names, fixed = TRUE)
  ]
  intercept_matches <- unique(intercept_matches)

  if (length(intercept_matches) != 1L) {
    .epix_stop(
      "Could not uniquely identify the bdlnm intercept by coefficient name. ",
      "Available fixed-effect names are: ",
      paste(raw_names, collapse = ", "), "."
    )
  }
  map["(Intercept)"] <- intercept_matches

  for (variable in metadata$vars) {
    basis <- metadata$basis_objects[[variable]]
    basis_names <- colnames(basis)
    standard_names <- .epix_cb_cols_for_var(metadata$cb_cols, variable)

    if (is.null(basis_names) || length(basis_names) != length(standard_names)) {
      .epix_stop(
        "Stored bdlnm crossbasis column names are unavailable or inconsistent ",
        "for variable '", variable, "'."
      )
    }

    for (j in seq_along(standard_names)) {
      standard <- standard_names[j]
      basis_name <- basis_names[j]

      exact_candidates <- unique(c(
        standard,
        paste0(variable, basis_name),
        basis_name
      ))

      matched <- .epix_match_one_name(
        exact_candidates,
        raw_names,
        paste0("bdlnm coefficient '", standard, "'")
      )

      if (!length(matched)) {
        suffix_candidates <- unique(c(
          paste0(variable, basis_name),
          basis_name
        ))
        suffix_matches <- raw_names[
          vapply(
            raw_names,
            function(x) any(vapply(suffix_candidates, function(s) endsWith(x, s), logical(1))),
            logical(1)
          )
        ]
        suffix_matches <- unique(suffix_matches)

        if (length(suffix_matches) == 1L) {
          matched <- suffix_matches
        } else if (length(suffix_matches) > 1L) {
          .epix_stop(
            "Ambiguous bdlnm coefficient mapping for canonical term '",
            standard, "': ", paste(suffix_matches, collapse = ", "), "."
          )
        }
      }

      if (!length(matched)) {
        .epix_stop(
          "Could not map canonical cross-basis coefficient '", standard,
          "' to a bdlnm fixed-effect coefficient using stored basis names."
        )
      }

      map[standard] <- matched
    }
  }

  expected <- .epix_expected_parameter_names(metadata)
  if (!identical(names(map), expected)) {
    map <- map[expected]
  }

  if (anyNA(map) || any(!nzchar(map)) || anyDuplicated(unname(map))) {
    .epix_stop("Invalid or duplicated bdlnm coefficient mapping.")
  }

  map
}

.epix_parameter_map <- function(fit, metadata, raw_names) {
  if (!is.character(raw_names) || !length(raw_names) || anyNA(raw_names) ||
      any(!nzchar(raw_names)) || anyDuplicated(raw_names)) {
    .epix_stop("Model coefficient names must be unique and non-empty.")
  }

  if (identical(metadata$engine, "bdlnm")) {
    return(.epix_bdlnm_parameter_map(fit, metadata, raw_names))
  }

  expected <- .epix_expected_parameter_names(metadata)
  map <- setNames(rep(NA_character_, length(expected)), expected)

  intercept <- .epix_match_one_name(
    c("(Intercept)", "Intercept"),
    raw_names,
    "the intercept"
  )
  if (!length(intercept)) {
    .epix_stop(
      "Could not identify the fitted intercept by name. Available names: ",
      paste(raw_names, collapse = ", "), "."
    )
  }
  map["(Intercept)"] <- intercept

  for (column in metadata$cb_cols) {
    if (!column %in% raw_names) {
      .epix_stop(
        "Fitted model is missing expected cross-basis coefficient '",
        column, "'. Coefficients are aligned by name, not by position."
      )
    }
    map[column] <- column
  }

  if (anyNA(map) || anyDuplicated(unname(map))) {
    .epix_stop("Invalid fixed-effect coefficient mapping.")
  }

  map
}

.epix_standardize_parameter_vector <- function(fit, metadata, raw_beta) {
  if (!is.numeric(raw_beta) || !length(raw_beta) || is.null(names(raw_beta))) {
    .epix_stop("Fixed-effect coefficients must be a named numeric vector.")
  }

  map <- .epix_parameter_map(fit, metadata, names(raw_beta))
  out <- as.numeric(raw_beta[unname(map)])
  names(out) <- names(map)

  if (any(!is.finite(out))) {
    .epix_stop("Extracted central fixed-effect coefficients are not all finite.")
  }

  out
}


# ----------------------------------------------------------------------------
# Central fixed/population parameters
# ----------------------------------------------------------------------------

#' Extract harmonized central fixed-effect parameters
#'
#' Frequentist models return fitted fixed-effect estimates. Bayesian models
#' return posterior means of fixed/population-level coefficients. The returned
#' vector is always named in the canonical EpiExposure design space:
#' `(Intercept)`, followed by `epiexposure_cb_cols`.
#'
#' @noRd
.extract_central_parameters <- function(fit) {
  metadata <- .get_epiexposure_metadata(fit)
  engine <- metadata$engine

  raw_beta <- switch(
    engine,

    glm = stats::coef(fit),

    glmmTMB = {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        .epix_stop("Package 'glmmTMB' is required for glmmTMB prediction.")
      }
      glmmTMB::fixef(fit)$cond
    },

    gam = stats::coef(fit),

    gamm = {
      if (!is.list(fit) || is.null(fit$gam) || !inherits(fit$gam, "gam")) {
        .epix_stop("Invalid stored `gamm` object: `$gam` component not found.")
      }
      stats::coef(fit$gam)
    },

    gls = stats::coef(fit),

    spamm = {
      if (!requireNamespace("spaMM", quietly = TRUE)) {
        .epix_stop("Package 'spaMM' is required for spaMM prediction.")
      }
      spaMM::fixef(fit)
    },

    brms = {
      if (!requireNamespace("brms", quietly = TRUE)) {
        .epix_stop("Package 'brms' is required for brms prediction.")
      }
      fixed <- brms::fixef(fit, summary = TRUE, robust = FALSE)
      beta <- fixed[, "Estimate"]
      names(beta) <- rownames(fixed)
      beta
    },

    inla = {
      if (is.null(fit$summary.fixed) || !nrow(fit$summary.fixed) ||
          is.null(rownames(fit$summary.fixed))) {
        .epix_stop("INLA fit does not contain named `summary.fixed` coefficients.")
      }
      beta <- fit$summary.fixed[, "mean"]
      names(beta) <- rownames(fit$summary.fixed)
      beta
    },

    bdlnm = {
      # Use the INLA marginal posterior means for deterministic prediction.
      # This avoids Monte Carlo noise from the finite coefficient sample matrix.
      if (is.null(fit$model) || is.null(fit$model$summary.fixed) ||
          !nrow(fit$model$summary.fixed) ||
          is.null(rownames(fit$model$summary.fixed))) {
        .epix_stop(
          "bdlnm fit does not contain named underlying INLA fixed-effect summaries."
        )
      }
      beta <- fit$model$summary.fixed[, "mean"]
      names(beta) <- rownames(fit$model$summary.fixed)
      beta
    },

    .epix_stop("Unsupported EpiExposure engine: ", engine, ".")
  )

  .epix_standardize_parameter_vector(fit, metadata, raw_beta)
}


# ----------------------------------------------------------------------------
# Frequentist fixed-effect covariance and MVN draws
# ----------------------------------------------------------------------------

.epix_extract_fixed_vcov <- function(fit, metadata) {
  engine <- metadata$engine

  raw_vcov <- switch(
    engine,

    glm = stats::vcov(fit),

    glmmTMB = {
      vc <- stats::vcov(fit)
      if (is.list(vc) && !is.null(vc$cond)) vc$cond else vc
    },

    gam = {
      vc <- tryCatch(
        stats::vcov(fit, unconditional = TRUE),
        error = function(e) NULL
      )
      if (is.null(vc)) stats::vcov(fit) else vc
    },

    gamm = {
      if (!is.list(fit) || is.null(fit$gam) || !inherits(fit$gam, "gam")) {
        .epix_stop("Invalid stored `gamm` object: `$gam` component not found.")
      }
      vc <- tryCatch(
        stats::vcov(fit$gam, unconditional = TRUE),
        error = function(e) NULL
      )
      if (is.null(vc)) stats::vcov(fit$gam) else vc
    },

    gls = stats::vcov(fit),

    spamm = {
      if (!requireNamespace("spaMM", quietly = TRUE)) {
        .epix_stop("Package 'spaMM' is required for spaMM uncertainty.")
      }
      stats::vcov(fit)
    },

    .epix_stop(
      "`.epix_extract_fixed_vcov()` is only used for frequentist engines."
    )
  )

  raw_vcov <- as.matrix(raw_vcov)

  if (nrow(raw_vcov) != ncol(raw_vcov) ||
      is.null(rownames(raw_vcov)) || is.null(colnames(raw_vcov))) {
    .epix_stop("Fixed-effect covariance matrix must be square and named.")
  }

  if (!identical(rownames(raw_vcov), colnames(raw_vcov))) {
    if (!setequal(rownames(raw_vcov), colnames(raw_vcov))) {
      .epix_stop("Fixed-effect covariance row and column names are inconsistent.")
    }
    raw_vcov <- raw_vcov[rownames(raw_vcov), rownames(raw_vcov), drop = FALSE]
  }

  map <- .epix_parameter_map(fit, metadata, rownames(raw_vcov))
  selected <- unname(map)
  vcov_standard <- raw_vcov[selected, selected, drop = FALSE]
  rownames(vcov_standard) <- names(map)
  colnames(vcov_standard) <- names(map)
  storage.mode(vcov_standard) <- "double"

  if (any(!is.finite(vcov_standard))) {
    .epix_stop("Fixed-effect covariance matrix contains non-finite values.")
  }

  # Numerical symmetry guard.
  vcov_standard <- (vcov_standard + t(vcov_standard)) / 2
  vcov_standard
}

.epix_mvn_draws <- function(mu, sigma, n_samples) {
  n_samples <- .epix_validate_n_samples(n_samples)

  if (!is.numeric(mu) || !length(mu) || is.null(names(mu)) ||
      any(!is.finite(mu))) {
    .epix_stop("MVN mean must be a finite named numeric vector.")
  }

  sigma <- as.matrix(sigma)
  if (nrow(sigma) != length(mu) || ncol(sigma) != length(mu) ||
      any(!is.finite(sigma))) {
    .epix_stop("MVN covariance matrix is incompatible with the coefficient vector.")
  }

  sigma <- (sigma + t(sigma)) / 2
  eig <- eigen(sigma, symmetric = TRUE)
  scale <- max(1, max(abs(eig$values)))
  tolerance <- 1e-8 * scale

  if (any(eig$values < -tolerance)) {
    .epix_stop(
      "Fixed-effect covariance matrix is not positive semidefinite; ",
      "parameter draws cannot be generated safely."
    )
  }

  values <- pmax(eig$values, 0)
  p <- length(mu)

  if (p == 1L) {
    draws <- matrix(
      stats::rnorm(n_samples, mean = mu, sd = sqrt(values[1])),
      ncol = 1L
    )
  } else {
    root <- eig$vectors %*%
      diag(sqrt(values), nrow = p, ncol = p)
    z <- matrix(stats::rnorm(n_samples * p), nrow = n_samples, ncol = p)
    draws <- z %*% t(root)
    draws <- sweep(draws, 2L, mu, FUN = "+")
  }

  colnames(draws) <- names(mu)
  rownames(draws) <- paste0("sample", seq_len(nrow(draws)))
  draws
}


# ----------------------------------------------------------------------------
# Bayesian fixed/population parameter draws
# ----------------------------------------------------------------------------

.epix_subsample_rows <- function(draws, n_samples) {
  n_samples <- .epix_validate_n_samples(n_samples)
  draws <- as.matrix(draws)

  if (nrow(draws) < n_samples) {
    .epix_stop(
      "Requested ", n_samples, " parameter draws, but only ", nrow(draws),
      " stored posterior draws are available. Reduce `n_samples` or refit the ",
      "model with more posterior draws."
    )
  }

  if (nrow(draws) == n_samples) return(draws)

  draws[sample.int(nrow(draws), n_samples, replace = FALSE), , drop = FALSE]
}

.epix_inla_fixed_draws <- function(inla_fit, n_samples) {
  n_samples <- .epix_validate_n_samples(n_samples)

  if (!requireNamespace("INLA", quietly = TRUE)) {
    .epix_stop("Package 'INLA' is required for INLA posterior draws.")
  }

  if (is.null(inla_fit$summary.fixed) || !nrow(inla_fit$summary.fixed) ||
      is.null(rownames(inla_fit$summary.fixed))) {
    .epix_stop("INLA fit does not contain named fixed-effect summaries.")
  }

  fixed_names <- rownames(inla_fit$summary.fixed)
  selection <- as.list(rep(1L, length(fixed_names)))
  names(selection) <- fixed_names

  posterior <- INLA::inla.posterior.sample(
    n = n_samples,
    result = inla_fit,
    selection = selection,
    add.names = TRUE
  )

  if (!is.list(posterior) || length(posterior) != n_samples) {
    .epix_stop("INLA returned an unexpected posterior-sample object.")
  }

  rows <- lapply(posterior, function(sample) {
    latent <- sample$latent
    if (is.null(latent)) {
      .epix_stop("An INLA posterior sample does not contain `$latent`.")
    }

    values <- as.numeric(latent)
    nm <- names(latent)
    if (is.null(nm)) nm <- rownames(latent)
    if (is.null(nm)) {
      .epix_stop("INLA posterior latent values do not have names.")
    }

    nm <- sub(":1$", "", nm)
    if (anyDuplicated(nm)) {
      .epix_stop("INLA posterior fixed-effect names are duplicated after normalization.")
    }

    names(values) <- nm

    missing <- setdiff(fixed_names, nm)
    if (length(missing)) {
      .epix_stop(
        "INLA posterior sample is missing fixed effect(s): ",
        paste(missing, collapse = ", "), "."
      )
    }

    values[fixed_names]
  })

  out <- do.call(rbind, rows)
  colnames(out) <- fixed_names
  rownames(out) <- paste0("sample", seq_len(nrow(out)))
  storage.mode(out) <- "double"

  if (any(!is.finite(out))) {
    .epix_stop("INLA fixed-effect posterior draws contain non-finite values.")
  }

  out
}

.epix_standardize_draw_matrix <- function(fit, metadata, raw_draws) {
  raw_draws <- as.matrix(raw_draws)

  if (!nrow(raw_draws) || !ncol(raw_draws) || is.null(colnames(raw_draws))) {
    .epix_stop("Parameter draws must be a matrix with named columns.")
  }

  map <- .epix_parameter_map(fit, metadata, colnames(raw_draws))
  out <- raw_draws[, unname(map), drop = FALSE]
  colnames(out) <- names(map)
  storage.mode(out) <- "double"

  if (any(!is.finite(out))) {
    .epix_stop("Parameter draws contain non-finite values.")
  }

  rownames(out) <- paste0("sample", seq_len(nrow(out)))
  out
}

#' Extract harmonized fixed/population parameter draws
#'
#' Frequentist engines use an asymptotic multivariate-normal distribution based
#' on the fitted fixed-effect coefficient vector and its covariance matrix.
#' Bayesian engines use their joint posterior fixed/population-effect draws.
#' Group-specific random effects and residual/process noise are excluded.
#'
#' A returned row is one coherent parameter draw and must be used jointly across
#' all prediction rows when constructing a curve, scenario set, baseline, or
#' contrast. This preserves covariance across predictions.
#'
#' @noRd
.extract_parameter_draws <- function(fit, n_samples = 1000L) {
  n_samples <- .epix_validate_n_samples(n_samples)
  metadata <- .get_epiexposure_metadata(fit)
  engine <- metadata$engine

  frequentist <- c("glm", "glmmTMB", "gam", "gamm", "gls", "spamm")

  if (engine %in% frequentist) {
    mu <- .extract_central_parameters(fit)
    sigma <- .epix_extract_fixed_vcov(fit, metadata)
    return(.epix_mvn_draws(mu, sigma, n_samples))
  }

  if (identical(engine, "brms")) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      .epix_stop("Package 'brms' is required for brms posterior draws.")
    }

    raw <- brms::fixef(fit, summary = FALSE)
    raw <- .epix_subsample_rows(raw, n_samples)
    return(.epix_standardize_draw_matrix(fit, metadata, raw))
  }

  if (identical(engine, "inla")) {
    raw <- .epix_inla_fixed_draws(fit, n_samples)
    return(.epix_standardize_draw_matrix(fit, metadata, raw))
  }

  if (identical(engine, "bdlnm")) {
    if (!is.null(fit$coefficients) && is.matrix(fit$coefficients) &&
        ncol(fit$coefficients) >= n_samples &&
        !is.null(rownames(fit$coefficients))) {
      # bdlnm stores coefficients as parameters x samples.
      raw <- t(fit$coefficients)
      raw <- .epix_subsample_rows(raw, n_samples)
      return(.epix_standardize_draw_matrix(fit, metadata, raw))
    }

    if (is.null(fit$model)) {
      .epix_stop(
        "bdlnm does not contain enough stored coefficient draws and its ",
        "underlying INLA model is unavailable."
      )
    }

    raw <- .epix_inla_fixed_draws(fit$model, n_samples)
    return(.epix_standardize_draw_matrix(fit, metadata, raw))
  }

  .epix_stop("Unsupported EpiExposure engine: ", engine, ".")
}


# ----------------------------------------------------------------------------
# Standard population-level fixed design
# ----------------------------------------------------------------------------

.epix_standard_fixed_design <- function(newdata, metadata) {
  if (!is.data.frame(newdata) || !nrow(newdata)) {
    .epix_stop("Prediction design must be a non-empty data.frame.")
  }

  missing <- setdiff(metadata$cb_cols, names(newdata))
  if (length(missing)) {
    .epix_stop(
      "Prediction design is missing fitted cross-basis column(s): ",
      paste(missing, collapse = ", "), "."
    )
  }

  x_cb <- newdata[, metadata$cb_cols, drop = FALSE]

  for (column in metadata$cb_cols) {
    if (!is.numeric(x_cb[[column]]) || anyNA(x_cb[[column]]) ||
        any(!is.finite(x_cb[[column]]))) {
      .epix_stop(
        "Prediction cross-basis column '", column,
        "' must contain only finite numeric values."
      )
    }
  }

  X <- cbind(
    "(Intercept)" = rep(1, nrow(newdata)),
    as.matrix(x_cb)
  )
  storage.mode(X) <- "double"

  expected <- .epix_expected_parameter_names(metadata)
  if (!identical(colnames(X), expected)) {
    .epix_stop("Internal prediction-design column ordering is inconsistent.")
  }

  X
}

.epix_manual_point_population <- function(fit, newdata, type) {
  metadata <- .get_epiexposure_metadata(fit)
  beta <- .extract_central_parameters(fit)
  X <- .epix_standard_fixed_design(newdata, metadata)

  if (!identical(colnames(X), names(beta))) {
    .epix_stop(
      "Prediction design and central parameter names do not align exactly."
    )
  }

  eta <- as.numeric(X %*% beta)

  if (identical(type, "link")) return(eta)

  link <- .epix_link_object(metadata$link)
  out <- as.numeric(link$linkinv(eta))

  if (any(!is.finite(out))) {
    .epix_stop("Inverse-link transformation produced non-finite predictions.")
  }

  out
}

.epix_gam_random_smooth_labels <- function(gam_fit, random_fit_col) {
  if (is.null(random_fit_col)) return(character(0))
  if (is.null(gam_fit$smooth) || !length(gam_fit$smooth)) return(character(0))

  labels <- vapply(
    gam_fit$smooth,
    function(sm) {
      term <- sm$term
      is_re <- inherits(sm, "random.effect")
      matches_term <- !is.null(term) && random_fit_col %in% term

      if (is_re && matches_term && !is.null(sm$label)) {
        as.character(sm$label)[1]
      } else {
        NA_character_
      }
    },
    character(1)
  )

  unique(labels[!is.na(labels) & nzchar(labels)])
}

.epix_native_point_population <- function(fit, newdata, type, metadata) {
  engine <- metadata$engine

  if (identical(engine, "glm")) {
    return(as.numeric(stats::predict(fit, newdata = newdata, type = type)))
  }

  if (identical(engine, "glmmTMB")) {
    return(as.numeric(stats::predict(
      fit,
      newdata = newdata,
      type = type,
      re.form = NA,
      allow.new.levels = TRUE
    )))
  }

  if (identical(engine, "gam")) {
    exclude <- character(0)

    if (identical(metadata$random_structure, "random_intercept")) {
      exclude <- .epix_gam_random_smooth_labels(
        fit,
        metadata$random_effect_fit_col
      )

      if (!length(exclude)) {
        .epix_stop(
          "Could not identify the fitted GAM random-effect smooth to exclude ",
          "for population-level prediction."
        )
      }
    }

    return(as.numeric(stats::predict(
      fit,
      newdata = newdata,
      type = type,
      exclude = exclude
    )))
  }

  if (identical(engine, "gamm")) {
    if (!is.list(fit) || is.null(fit$gam) || !inherits(fit$gam, "gam")) {
      .epix_stop("Invalid `gamm` object: `$gam` component not found.")
    }

    # mgcv documents prediction from the returned $gam component as having the
    # gamm random effects set to zero, which is exactly the v1 population target.
    return(as.numeric(stats::predict(
      fit$gam,
      newdata = newdata,
      type = type
    )))
  }

  if (identical(engine, "gls")) {
    # EpiExposure restricts GLS to Gaussian identity-link models.
    return(as.numeric(stats::predict(fit, newdata = newdata)))
  }

  if (identical(engine, "spamm")) {
    if (!requireNamespace("spaMM", quietly = TRUE)) {
      .epix_stop("Package 'spaMM' is required for spaMM prediction.")
    }

    return(as.numeric(stats::predict(
      fit,
      newdata = newdata,
      type = type,
      re.form = NA
    )))
  }

  .epix_stop(
    "No native deterministic point-prediction route is used for engine '",
    engine, "'."
  )
}

#' Predict the harmonized population-level point estimate
#'
#' Uses native engine prediction for frequentist engines when it represents the
#' same estimand (expected response or linear predictor with random effects set
#' to zero). If that native route fails, an algebraically equivalent fixed-effect
#' `X %*% beta` route is used. Bayesian engines deliberately use the central
#' fixed/population parameters and `X %*% beta` rather than summarizing
#' posterior expected predictions, so `uncertainty = FALSE` has the same
#' parameter-estimate meaning across engines.
#'
#' @param fit Fitted EpiExposure model.
#' @param newdata Epidemic-level design containing canonical `cb_*` columns,
#'   generally returned by `.build_newdata_basis()`.
#' @param type `"response"` or `"link"`.
#' @param warn_fallback Logical; warn if a known frequentist native prediction
#'   fails and the equivalent manual fixed-effect route is used.
#'
#' @noRd
.predict_point_population <- function(
    fit,
    newdata,
    type = c("response", "link"),
    warn_fallback = FALSE
) {
  type <- match.arg(type)
  metadata <- .get_epiexposure_metadata(fit)

  if (!is.data.frame(newdata) || !nrow(newdata)) {
    .epix_stop("`newdata` must be a non-empty epidemic-level prediction design.")
  }

  # Validate canonical design even when native prediction will be used.
  invisible(.epix_standard_fixed_design(newdata, metadata))

  frequentist <- c("glm", "glmmTMB", "gam", "gamm", "gls", "spamm")

  if (metadata$engine %in% frequentist) {
    native_error <- NULL
    native <- tryCatch(
      .epix_native_point_population(fit, newdata, type, metadata),
      error = function(e) {
        native_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(native) && length(native) == nrow(newdata) &&
        all(is.finite(native))) {
      attr(native, "epiexposure_prediction_method") <- "native_population"
      attr(native, "epiexposure_prediction_type") <- type
      return(native)
    }

    if (isTRUE(warn_fallback)) {
      warning(
        "Native population-level prediction for engine '", metadata$engine,
        "' failed; using the equivalent fixed-effect X-beta route. ",
        if (!is.null(native_error)) paste0("Reason: ", native_error) else "",
        call. = FALSE
      )
    }
  }

  out <- .epix_manual_point_population(fit, newdata, type)
  attr(out, "epiexposure_prediction_method") <- "central_parameters_xbeta"
  attr(out, "epiexposure_prediction_type") <- type
  out
}


# ----------------------------------------------------------------------------
# Draw-by-draw population-level prediction
# ----------------------------------------------------------------------------

#' Predict population-level expected responses draw by draw
#'
#' Parameter draws are applied jointly to all rows in `newdata`. The same draw
#' index therefore represents the same coefficient vector across the entire
#' prediction grid, preserving covariance among baseline, target, scenario, and
#' lag/exposure predictions.
#'
#' No group-specific random-effect draws and no residual/process noise are
#' included. For `type = "response"`, the inverse link is applied separately to
#' every parameter draw before any summary is computed.
#'
#' @noRd
.predict_draws_population <- function(
    fit,
    newdata,
    n_samples = 1000L,
    type = c("response", "link"),
    parameter_draws = NULL
) {
  type <- match.arg(type)
  n_samples <- .epix_validate_n_samples(n_samples)
  metadata <- .get_epiexposure_metadata(fit)
  X <- .epix_standard_fixed_design(newdata, metadata)

  if (is.null(parameter_draws)) {
    parameter_draws <- .extract_parameter_draws(fit, n_samples = n_samples)
  } else {
    parameter_draws <- as.matrix(parameter_draws)

    if (nrow(parameter_draws) != n_samples) {
      .epix_stop(
        "Supplied `parameter_draws` has ", nrow(parameter_draws),
        " row(s), but `n_samples = ", n_samples, "`."
      )
    }

    expected <- .epix_expected_parameter_names(metadata)
    if (is.null(colnames(parameter_draws)) ||
        !identical(colnames(parameter_draws), expected)) {
      .epix_stop(
        "Supplied `parameter_draws` must have canonical columns in this exact ",
        "order: ", paste(expected, collapse = ", "), "."
      )
    }

    storage.mode(parameter_draws) <- "double"
    if (any(!is.finite(parameter_draws))) {
      .epix_stop("Supplied `parameter_draws` contain non-finite values.")
    }
  }

  if (!identical(colnames(parameter_draws), colnames(X))) {
    .epix_stop("Parameter draws and prediction design do not align exactly.")
  }

  eta <- parameter_draws %*% t(X)
  storage.mode(eta) <- "double"

  if (identical(type, "link")) {
    out <- eta
  } else {
    link <- .epix_link_object(metadata$link)
    transformed <- link$linkinv(as.vector(eta))
    out <- matrix(
      as.numeric(transformed),
      nrow = nrow(eta),
      ncol = ncol(eta),
      byrow = FALSE
    )
  }

  if (any(!is.finite(out))) {
    .epix_stop(
      "Draw-by-draw prediction produced non-finite values. Check the fitted ",
      "model, link, and extrapolation range."
    )
  }

  rownames(out) <- paste0("sample", seq_len(nrow(out)))
  colnames(out) <- paste0("prediction", seq_len(ncol(out)))
  attr(out, "epiexposure_prediction_type") <- type
  attr(out, "epiexposure_prediction_level") <- "population"
  attr(out, "epiexposure_prediction_estimand") <- "expected_response"
  attr(out, "epiexposure_parameter_draws") <- parameter_draws
  out
}


# ----------------------------------------------------------------------------
# Shared empirical summary of draw-based expected predictions
# ----------------------------------------------------------------------------

#' Summarize draw-by-draw EpiExposure predictions
#'
#' Summaries are computed independently for each prediction column *after* all
#' draw-specific transformations have been performed. The central uncertainty
#' summary is the median, with standard deviation and empirical lower/upper
#' quantiles. This helper does not create new draws and does not add residual
#' variability.
#'
#' @noRd
.summarise_prediction_draws <- function(
    draws,
    probs = c(0.025, 0.975),
    value_name = "prediction"
) {
  probs <- .epix_validate_probs(probs)

  if (!.epix_is_scalar_string(value_name)) {
    .epix_stop("`value_name` must be one non-empty character value.")
  }

  draws <- as.matrix(draws)
  storage.mode(draws) <- "double"

  if (!nrow(draws) || !ncol(draws) || any(!is.finite(draws))) {
    .epix_stop("`draws` must be a non-empty finite numeric matrix.")
  }

  med <- apply(draws, 2L, stats::median)
  sdv <- if (nrow(draws) > 1L) {
    apply(draws, 2L, stats::sd)
  } else {
    rep(NA_real_, ncol(draws))
  }
  lower <- apply(draws, 2L, stats::quantile, probs = probs[1], names = FALSE)
  upper <- apply(draws, 2L, stats::quantile, probs = probs[2], names = FALSE)

  out <- data.frame(
    row_id = seq_len(ncol(draws)),
    stringsAsFactors = FALSE
  )
  out[[value_name]] <- as.numeric(med)
  out[[paste0(value_name, "_sd")]] <- as.numeric(sdv)
  out[[paste0(value_name, "_lower")]] <- as.numeric(lower)
  out[[paste0(value_name, "_upper")]] <- as.numeric(upper)

  attr(out, "epiexposure_probs") <- probs
  attr(out, "epiexposure_summary_center") <- "median"
  attr(out, "epiexposure_n_samples") <- nrow(draws)
  out
}


# ----------------------------------------------------------------------------
# Optional convenience accessor for prediction keys
# ----------------------------------------------------------------------------

.epiexposure_prediction_keys <- function(prediction_design) {
  keys <- attr(
    prediction_design,
    "epiexposure_prediction_keys",
    exact = TRUE
  )

  if (is.null(keys)) {
    data.frame(row_id = seq_len(nrow(prediction_design)))
  } else {
    keys
  }
}
