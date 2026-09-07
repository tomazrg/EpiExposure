#' Build an epidemic-level DLNM design matrix
#'
#' Converts complete long-format exposure histories into one epidemic-level
#' design row per `epi_id`. For each exposure, the function reconstructs a new
#' DLNM cross-basis from that epidemic's observed exposure history using the
#' **effective basis parameterization stored in the supplied training template**,
#' then retains the final cross-basis row after chronological ordering.
#'
#' `build_design()` never re-estimates spline knots, boundary knots, degrees of
#' freedom, or lag-basis parameters from the input `data`. The supplied
#' `cb_templates` are authoritative for basis transport.
#'
#' @param data Non-empty long-format data frame containing `epi_id`, `time`, and
#'   every exposure represented in `cb_templates`. Rows are ordered internally
#'   by `epi_id` and `time`. Under the EpiExposure exact-history contract, every
#'   epidemic must contain exactly the common fitted `max_lag + 1` time points.
#'   Because all exposure variables are finite columns of those same rows, all
#'   fitted exposures necessarily use the same temporal support and history
#'   length.
#'
#'   If a column `y` is present and `include_response = TRUE`, `y` must contain
#'   exactly one finite numeric outcome value per epidemic.
#' @param cb_templates Named list of `dlnm::crossbasis` objects, normally
#'   returned by `define_exposures()`. Every element must contain valid effective
#'   `argvar`, `arglag`, and `lag` attributes.
#'
#'   The cross-basis objects themselves are not reused numerically for new
#'   epidemics. Their stored effective parameterization is reused to transform
#'   each new exposure history.
#' @param max_lag Optional lag validation argument retained for backward
#'   compatibility. `cb_templates` are the authoritative source of the fitted
#'   lag definition and `max_lag` never overrides them.
#'
#'   All templates must share one common maximum lag. If `NULL`, that common lag
#'   is read directly from the templates. If supplied, `max_lag` may be one
#'   non-negative integer scalar or the legacy two-element form `c(0, L)`;
#'   after normalization its maximum must equal the common template lag.
#'   Exposure-specific/named lag windows are not supported by the current
#'   EpiExposure exact-history contract.
#' @param include_response Logical. If `TRUE` and `data` contains `y`, append
#'   one response value per epidemic. If `FALSE`, `y` is ignored.
#'
#' @return A data frame with one row per epidemic containing `epi_id`,
#'   canonical cross-basis columns named `cb_<exposure>_<column_index>`, and,
#'   optionally, `y`.
#'
#'   The returned data frame stores the following attributes:
#'
#'   - `"cb_templates"`: the original supplied cross-basis templates;
#'   - `"epiexposure_spec"`: effective per-exposure basis specification
#'     synchronized to the template attributes;
#'   - `"epiexposure_cb_cols"`: canonical cross-basis column order;
#'   - `"epiexposure_vars"`: exposure names;
#'   - `"epiexposure_max_lag"`: common fitted maximum lag;
#'   - `"epiexposure_history_length"`: exact required history length,
#'     `max_lag + 1`;
#'   - `"epiexposure_history_contract"`:
#'     `"all_fitted_exposures_same_exact_max_lag_plus_one"`;
#'   - `"epiexposure_time_step"`: common time-series spacing, or `NA` when
#'     `max_lag = 0` because a one-point history has no estimable spacing;
#'   - `"epiexposure_design_contract"`:
#'     `"final_crossbasis_row_per_group"`.
#'
#' @details
#' ## Temporal requirements
#'
#' A vector supplied to `dlnm::crossbasis()` is interpreted as one complete,
#' equally spaced exposure history. Accordingly, `build_design()` requires
#' finite numeric time, unique times within epidemics, constant spacing within
#' each epidemic when more than one time point is present, and the same spacing
#' across epidemics.
#'
#' If the common fitted maximum lag is `L`, every epidemic must contain exactly
#'
#' \deqn{
#'   L + 1
#' }
#'
#' observations. Histories with fewer or more observations are rejected.
#' `build_design()` never truncates an older history, selects a trailing window,
#' pads a shorter history, or silently aligns exposures with different lag
#' windows. Since every exposure is evaluated on the same rows of `data`, all
#' fitted exposures have identical temporal length and support.
#'
#' For `max_lag = 0`, the exact history length is one observation. In that
#' special case no within-history time interval exists, so
#' `epiexposure_time_step` is stored as `NA`.
#'
#' ## Training-template transport
#'
#' For each exposure `x`, reconstruction is conceptually:
#'
#' ```
#' dlnm::crossbasis(
#'   x,
#'   lag    = attr(training_template, "lag"),
#'   argvar = attr(training_template, "argvar"),
#'   arglag = attr(training_template, "arglag")
#' )
#' ```
#'
#' The numerical cross-basis values change with the new exposure history while
#' the training basis definition remains fixed. Reconstructed bases must retain
#' the training template's number of columns and native column names/order
#' before canonical EpiExposure names are assigned.
#'
#' ## Template storage
#'
#' `attr(result, "cb_templates")` stores the **original supplied templates**,
#' not a cross-basis reconstructed from the first epidemic. This preserves the
#' original training basis definition for later prediction and validation.
#'
#' @export
build_design <- function(
    data,
    cb_templates,
    max_lag = NULL,
    include_response = TRUE
) {

  valid_flag <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  is_integerish <- function(x) {
    is.numeric(x) &&
      length(x) >= 1L &&
      !anyNA(x) &&
      all(is.finite(x)) &&
      all(x >= 0) &&
      all(x == as.integer(x))
  }

  same_numeric_scalar <- function(a, b, tolerance = 1e-10) {
    is.numeric(a) &&
      length(a) == 1L &&
      is.finite(a) &&
      is.numeric(b) &&
      length(b) == 1L &&
      is.finite(b) &&
      abs(a - b) <= tolerance * max(1, abs(a), abs(b))
  }

  # ==========================================================================
  # BASIC INPUT VALIDATION
  # ==========================================================================

  if (!is.data.frame(data) || !nrow(data)) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }

  if (!valid_flag(include_response)) {
    stop("`include_response` must be TRUE or FALSE.", call. = FALSE)
  }

  missing_basic <- setdiff(c("epi_id", "time"), names(data))
  if (length(missing_basic)) {
    stop(
      "`data` is missing required column(s): ",
      paste(missing_basic, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (anyNA(data$epi_id)) {
    stop("`data$epi_id` cannot contain missing values.", call. = FALSE)
  }

  if (!is.numeric(data$time) ||
      anyNA(data$time) ||
      any(!is.finite(data$time))) {
    stop(
      "`data$time` must contain only finite numeric values.",
      call. = FALSE
    )
  }

  if (!is.list(cb_templates) ||
      !length(cb_templates) ||
      is.null(names(cb_templates)) ||
      anyNA(names(cb_templates)) ||
      any(!nzchar(names(cb_templates))) ||
      anyDuplicated(names(cb_templates))) {
    stop(
      "`cb_templates` must be a non-empty named list with unique, non-empty ",
      "exposure names.",
      call. = FALSE
    )
  }

  vars <- names(cb_templates)

  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop(
      "`data` is missing exposure variable(s) required by `cb_templates`: ",
      paste(missing_vars, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  for (variable in vars) {
    x <- data[[variable]]
    if (!is.numeric(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "Exposure variable '",
        variable,
        "' must contain only finite numeric values.",
        call. = FALSE
      )
    }
  }

  if (!requireNamespace("dlnm", quietly = TRUE)) {
    stop(
      "Package 'dlnm' is required by `build_design()`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # VALIDATE TRAINING TEMPLATE DEFINITIONS
  # ==========================================================================

  stored_spec <- attr(cb_templates, "spec", exact = TRUE)

  if (!is.null(stored_spec)) {
    if (!is.list(stored_spec) ||
        is.null(names(stored_spec)) ||
        anyNA(names(stored_spec)) ||
        any(!nzchar(names(stored_spec))) ||
        anyDuplicated(names(stored_spec)) ||
        !setequal(names(stored_spec), vars)) {
      stop(
        "`attr(cb_templates, 'spec')` is present but inconsistent with the ",
        "template names.",
        call. = FALSE
      )
    }
    stored_spec <- stored_spec[vars]
  }

  template_info <- stats::setNames(vector("list", length(vars)), vars)
  effective_spec <- if (!is.null(stored_spec)) {
    stored_spec
  } else {
    stats::setNames(vector("list", length(vars)), vars)
  }

  cb_cols <- character(0)
  template_max_lag <- stats::setNames(integer(length(vars)), vars)

  for (variable in vars) {
    template <- cb_templates[[variable]]

    if (!inherits(template, "crossbasis")) {
      stop(
        "`cb_templates[['",
        variable,
        "']]` must inherit from `crossbasis`.",
        call. = FALSE
      )
    }

    argvar <- attr(template, "argvar", exact = TRUE)
    arglag <- attr(template, "arglag", exact = TRUE)
    lag_attr <- attr(template, "lag", exact = TRUE)

    if (!is.list(argvar) || !length(argvar)) {
      stop(
        "Training cross-basis for exposure '",
        variable,
        "' is missing a valid effective `argvar` attribute.",
        call. = FALSE
      )
    }

    if (!is.list(arglag) || !length(arglag)) {
      stop(
        "Training cross-basis for exposure '",
        variable,
        "' is missing a valid effective `arglag` attribute.",
        call. = FALSE
      )
    }

    if (!is.numeric(lag_attr) ||
        !length(lag_attr) ||
        anyNA(lag_attr) ||
        any(!is.finite(lag_attr)) ||
        any(lag_attr < 0) ||
        any(lag_attr != as.integer(lag_attr))) {
      stop(
        "Training cross-basis for exposure '",
        variable,
        "' is missing a valid non-negative integer `lag` attribute.",
        call. = FALSE
      )
    }

    lag_attr <- as.integer(lag_attr)

    if (min(lag_attr) != 0L) {
      stop(
        "EpiExposure v1 requires lag histories to start at lag 0. ",
        "Template for exposure '",
        variable,
        "' starts at lag ",
        min(lag_attr),
        ".",
        call. = FALSE
      )
    }

    current_max_lag <- as.integer(max(lag_attr))

    if (ncol(template) < 1L) {
      stop(
        "Training cross-basis for exposure '",
        variable,
        "' contains no basis columns.",
        call. = FALSE
      )
    }

    native_names <- colnames(template)

    if (!is.null(native_names) &&
        (length(native_names) != ncol(template) ||
         anyNA(native_names) ||
         any(!nzchar(native_names)) ||
         anyDuplicated(native_names))) {
      stop(
        "Training cross-basis for exposure '",
        variable,
        "' has invalid native column names.",
        call. = FALSE
      )
    }

    canonical_names <- paste0(
      "cb_",
      variable,
      "_",
      seq_len(ncol(template))
    )

    if (!is.null(stored_spec)) {
      current_spec <- stored_spec[[variable]]

      if (!is.list(current_spec)) {
        stop(
          "Stored exposure specification for '",
          variable,
          "' must be a list.",
          call. = FALSE
        )
      }

      if (!is.null(current_spec$max_lag)) {
        spec_lag <- current_spec$max_lag

        if (!is.numeric(spec_lag) ||
            length(spec_lag) != 1L ||
            is.na(spec_lag) ||
            !is.finite(spec_lag) ||
            spec_lag < 0 ||
            spec_lag != as.integer(spec_lag)) {
          stop(
            "Stored `max_lag` specification for exposure '",
            variable,
            "' must be one non-negative finite integer.",
            call. = FALSE
          )
        }

        if (!identical(as.integer(spec_lag), current_max_lag)) {
          stop(
            "Stored exposure specification and cross-basis template disagree ",
            "on `max_lag` for exposure '",
            variable,
            "'.",
            call. = FALSE
          )
        }
      }
    }

    effective_spec[[variable]]$max_lag <- current_max_lag
    effective_spec[[variable]]$argvar <- argvar
    effective_spec[[variable]]$arglag <- arglag

    template_max_lag[[variable]] <- current_max_lag

    template_info[[variable]] <- list(
      lag = lag_attr,
      max_lag = current_max_lag,
      argvar = argvar,
      arglag = arglag,
      ncol = ncol(template),
      native_names = native_names,
      canonical_names = canonical_names
    )

    cb_cols <- c(cb_cols, canonical_names)
  }

  if (anyDuplicated(cb_cols)) {
    stop(
      "Canonical cross-basis column names are duplicated across exposures.",
      call. = FALSE
    )
  }

  if (length(unique(template_max_lag)) != 1L) {
    stop(
      "All exposure templates must use the same `max_lag` under the ",
      "EpiExposure exact-history contract. Template values are: ",
      paste(
        paste0(names(template_max_lag), "=", template_max_lag),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  common_max_lag <- unname(template_max_lag[[1L]])
  required_history_length <- common_max_lag + 1L
  history_contract <- "all_fitted_exposures_same_exact_max_lag_plus_one"

  # ==========================================================================
  # OPTIONAL max_lag VALIDATION
  # ==========================================================================

  if (!is.null(max_lag)) {
    if (!is_integerish(max_lag)) {
      stop(
        "`max_lag` must be NULL, one non-negative integer, or the legacy ",
        "two-element form `c(0, L)`.",
        call. = FALSE
      )
    }

    if (!is.null(names(max_lag)) &&
        any(nzchar(names(max_lag)))) {
      stop(
        "Named/exposure-specific `max_lag` values are not supported. All ",
        "exposure templates must share one common lag window.",
        call. = FALSE
      )
    }

    if (length(max_lag) == 1L) {
      supplied_common_lag <- as.integer(max_lag)

    } else if (length(max_lag) == 2L &&
               min(max_lag) == 0L) {
      supplied_common_lag <- as.integer(max(max_lag))

    } else {
      stop(
        "`max_lag` must be one non-negative integer or the legacy two-element ",
        "form `c(0, L)`.",
        call. = FALSE
      )
    }

    if (!identical(supplied_common_lag, common_max_lag)) {
      stop(
        "`max_lag` contradicts the common training-template lag: supplied ",
        supplied_common_lag,
        ", template maximum lag ",
        common_max_lag,
        ".",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # TEMPORAL REGULARITY AND COVERAGE
  # ==========================================================================

  data_ordered <- data[
    order(data$epi_id, data$time),
    ,
    drop = FALSE
  ]
  rownames(data_ordered) <- NULL

  group_key <- as.character(data_ordered$epi_id)
  group_rows <- split(seq_len(nrow(data_ordered)), group_key, drop = TRUE)

  first_rows <- vapply(group_rows, function(idx) idx[1L], integer(1))
  group_rows <- group_rows[order(first_rows)]
  first_rows <- vapply(group_rows, function(idx) idx[1L], integer(1))
  group_ids <- data_ordered$epi_id[first_rows]

  common_time_step <- NA_real_
  tolerance <- sqrt(.Machine$double.eps)

  for (i in seq_along(group_rows)) {
    idx <- group_rows[[i]]
    current_id <- data_ordered$epi_id[idx[1L]]
    tt <- data_ordered$time[idx]

    if (anyDuplicated(tt)) {
      stop(
        "Duplicated `time` values were detected within epi_id = ",
        as.character(current_id),
        ".",
        call. = FALSE
      )
    }

    if (length(tt) != required_history_length) {
      stop(
        "epi_id = ",
        as.character(current_id),
        " must contain exactly ",
        required_history_length,
        " time observation(s) for max_lag = ",
        common_max_lag,
        "; received ",
        length(tt),
        ". Histories are not truncated, padded, or silently realigned.",
        call. = FALSE
      )
    }

    if (required_history_length > 1L) {
      dt <- diff(tt)

      if (any(dt <= 0) || any(!is.finite(dt))) {
        stop(
          "Time must increase strictly within epi_id = ",
          as.character(current_id),
          ".",
          call. = FALSE
        )
      }

      current_step <- dt[1L]

      if (any(
        abs(dt - current_step) >
        tolerance * pmax(1, abs(dt), abs(current_step))
      )) {
        stop(
          "Irregular time spacing was detected within epi_id = ",
          as.character(current_id),
          ". DLNM vector histories must be complete and equally spaced.",
          call. = FALSE
        )
      }

      if (is.na(common_time_step)) {
        common_time_step <- as.numeric(current_step)
      } else if (!same_numeric_scalar(
        current_step,
        common_time_step,
        tolerance = tolerance
      )) {
        stop(
          "Time spacing differs across epidemics. epi_id = ",
          as.character(current_id),
          " uses step ",
          format(current_step),
          ", while the common step is ",
          format(common_time_step),
          ".",
          call. = FALSE
        )
      }
    }
  }

  # ==========================================================================
  # BUILD ONE DESIGN ROW PER EPIDEMIC
  # ==========================================================================

  out <- data.frame(
    epi_id = group_ids,
    check.names = FALSE
  )

  for (variable in vars) {
    detail <- template_info[[variable]]

    variable_matrix <- matrix(
      NA_real_,
      nrow = length(group_rows),
      ncol = detail$ncol
    )
    colnames(variable_matrix) <- detail$canonical_names

    for (i in seq_along(group_rows)) {
      idx <- group_rows[[i]]
      x <- data_ordered[[variable]][idx]

      if (length(x) != required_history_length) {
        stop(
          "Exposure history for '",
          variable,
          "' and epi_id = ",
          as.character(data_ordered$epi_id[idx[1L]]),
          " must contain exactly ",
          required_history_length,
          " observations for max_lag = ",
          common_max_lag,
          "; received ",
          length(x),
          ".",
          call. = FALSE
        )
      }

      cb_new <- dlnm::crossbasis(
        x,
        lag = detail$lag,
        argvar = detail$argvar,
        arglag = detail$arglag
      )

      if (!inherits(cb_new, "crossbasis")) {
        stop(
          "Reconstructed basis for exposure '",
          variable,
          "' did not return a `crossbasis` object.",
          call. = FALSE
        )
      }

      if (ncol(cb_new) != detail$ncol) {
        stop(
          "Reconstructed cross-basis dimension mismatch for exposure '",
          variable,
          "': training template has ",
          detail$ncol,
          " column(s), reconstructed basis has ",
          ncol(cb_new),
          ".",
          call. = FALSE
        )
      }

      reconstructed_native_names <- colnames(cb_new)

      if (!is.null(detail$native_names) &&
          (is.null(reconstructed_native_names) ||
           !identical(reconstructed_native_names, detail$native_names))) {
        stop(
          "Reconstructed cross-basis native column names/order do not match ",
          "the training template for exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      reconstructed_lag <- attr(cb_new, "lag", exact = TRUE)

      if (is.null(reconstructed_lag) ||
          !identical(
            as.integer(reconstructed_lag),
            as.integer(detail$lag)
          )) {
        stop(
          "Reconstructed lag definition does not match the training template ",
          "for exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      final_row <- as.numeric(
        cb_new[nrow(cb_new), , drop = TRUE]
      )

      if (length(final_row) != detail$ncol ||
          anyNA(final_row) ||
          any(!is.finite(final_row))) {
        stop(
          "The final cross-basis row is not finite for exposure '",
          variable,
          "' and epi_id = ",
          as.character(data_ordered$epi_id[idx[1L]]),
          ".",
          call. = FALSE
        )
      }

      variable_matrix[i, ] <- final_row
    }

    for (column in detail$canonical_names) {
      out[[column]] <- variable_matrix[, column]
    }
  }

  actual_cb_cols <- grep("^cb_", names(out), value = TRUE)

  if (!identical(actual_cb_cols, cb_cols)) {
    stop(
      "Internal cross-basis column order does not match the template contract.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # OPTIONAL RESPONSE
  # ==========================================================================

  if (isTRUE(include_response) && "y" %in% names(data_ordered)) {
    if (!is.numeric(data_ordered$y) ||
        anyNA(data_ordered$y) ||
        any(!is.finite(data_ordered$y))) {
      stop(
        "When included, `data$y` must contain only finite numeric values.",
        call. = FALSE
      )
    }

    y_out <- numeric(length(group_rows))

    for (i in seq_along(group_rows)) {
      idx <- group_rows[[i]]
      yy <- unique(data_ordered$y[idx])

      if (length(yy) != 1L) {
        stop(
          "Multiple distinct `y` values were found within epi_id = ",
          as.character(data_ordered$epi_id[idx[1L]]),
          ". Exactly one epidemic-level outcome is required.",
          call. = FALSE
        )
      }

      y_out[i] <- yy
    }

    out$y <- y_out
  }

  # ==========================================================================
  # OUTPUT METADATA
  # ==========================================================================

  attr(out, "cb_templates") <- cb_templates
  attr(out, "epiexposure_spec") <- effective_spec
  attr(out, "epiexposure_cb_cols") <- cb_cols
  attr(out, "epiexposure_vars") <- vars
  attr(out, "epiexposure_max_lag") <- common_max_lag
  attr(out, "epiexposure_history_length") <- required_history_length
  attr(out, "epiexposure_history_contract") <- history_contract
  attr(out, "epiexposure_time_step") <- common_time_step
  attr(out, "epiexposure_design_contract") <-
    "final_crossbasis_row_per_group"
  attr(out, "epiexposure_profile_order") <- "chronological"

  out
}
