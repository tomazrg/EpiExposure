#' Generate exposure-value scenarios across DLNM lag periods
#'
#' Creates structured scenario definitions for `simulate_scenarios()` by
#' generating combinations of user-specified exposure values and copying each
#' resulting combination unchanged to every supplied lag period.
#'
#' This function generates **scenario definitions only**. It does not fit a
#' model, construct a cross-basis, or predict an outcome. Conversion from
#' retrospective lag-period definitions to complete chronological exposure
#' histories is performed later by `simulate_scenarios()`, which then calls
#' `predict_outcomes()`.
#'
#' @param periods A non-empty data frame, typically returned by
#'   `define_periods()`, containing `period`, `lag_start`, and `lag_end`.
#'
#'   Periods are expressed on the retrospective lag scale. For example,
#'   `lag_start = 0` and `lag_end = 20` represent lags 0 through 20, the portion
#'   of the exposure history closest to outcome assessment.
#'
#'   Periods must not overlap. Gaps are allowed: any exposure-lag positions not
#'   covered later by a scenario remain available for background filling in
#'   `simulate_scenarios()`.
#' @param vary Named list of exposure variables and finite numeric values to
#'   evaluate. At least one varying variable is required.
#'
#'   With `scenario_type = "grid"`, the Cartesian product across the variables
#'   in `vary` is generated **once**. Each resulting combination is then copied
#'   to every supplied lag period. No additional Cartesian product is generated
#'   across periods.
#'
#'   With `scenario_type = "paired"`, values are matched by position across
#'   variables. All vectors in `vary` must therefore have the same length.
#'
#'   Values within each varying variable must be unique. Repeated values would
#'   create duplicated scenario combinations and are rejected explicitly.
#' @param fixed Named list of exposure variables held fixed in every generated
#'   scenario and every supplied period. Each fixed variable must contain
#'   exactly one finite numeric value. Use `list()` when none are required.
#' @param scenario_type Character. `"grid"` or `"paired"`.
#' @param scenario_names Optional character vector containing exactly one unique,
#'   non-empty name for each generated scenario.
#' @param scenario_var Optional name of one variable in `vary`. It is used only
#'   to construct the automatic scenario-name prefix and is retained as
#'   metadata. It does **not** give that exposure a special statistical role.
#'   If `NULL`, the first variable in `vary` is used.
#'
#' @return A structured list containing:
#'
#'   \describe{
#'     \item{`scenarios`}{Named scenario definitions. Each generated exposure
#'       combination is copied to every supplied period.}
#'     \item{`periods`}{The validated lag-period definitions.}
#'     \item{`info`}{One metadata row per generated scenario. It contains the
#'       scenario name and the values of both varying and fixed exposures.}
#'   }
#'
#'   Additional attributes record `scenario_type`, `scenario_var`, and simple
#'   period-coverage diagnostics.
#'
#' @details
#' ## How values are repeated across periods
#'
#' `simulate_ranges()` deliberately does **not** generate separate exposure
#' values for each period. Suppose four periods cover lags 0--85 and one grid
#' combination is:
#'
#' ```
#' tmean   = 20
#' rain    = 8
#' wetness = 2
#' ```
#'
#' The generated scenario is equivalent to:
#'
#' ```
#' P1 = list(tmean = 20, rain = 8, wetness = 2)
#' P2 = list(tmean = 20, rain = 8, wetness = 2)
#' P3 = list(tmean = 20, rain = 8, wetness = 2)
#' P4 = list(tmean = 20, rain = 8, wetness = 2)
#' ```
#'
#' `simulate_scenarios()` then expands each scalar across all lags belonging to
#' its period. If P1--P4 jointly cover the entire fitted lag window and the
#' variables in `vary` plus `fixed` include every exposure in the fitted model,
#' the resulting exposure history is complete and no background/reference value
#' is needed for those positions.
#'
#' If periods leave gaps, do not extend to the fitted maximum lag, or omit one
#' or more fitted exposures, `simulate_scenarios()` fills only the unassigned
#' positions from `ref_vals` or from exposure medians derived from `data`.
#'
#' ## Grid versus paired scenarios
#'
#' With:
#'
#' ```
#' vary = list(
#'   tmean = c(20, 34),
#'   rain = 0:15,
#'   wetness = c(2, 10)
#' )
#' ```
#'
#' `scenario_type = "grid"` generates
#'
#' \deqn{2 \times 16 \times 2 = 64}
#'
#' scenarios. If four periods are supplied, there are still 64 scenarios, not
#' \eqn{64^4}: each of the 64 combinations is simply copied to all four
#' periods.
#'
#' `scenario_type = "paired"` instead combines the first value of each varying
#' exposure, then the second value of each, and so forth.
#'
#' @export
simulate_ranges <- function(
    periods,
    vary = list(),
    fixed = list(),
    scenario_type = c("grid", "paired"),
    scenario_names = NULL,
    scenario_var = NULL
) {

  scenario_type <- match.arg(scenario_type)

  # ==========================================================================
  # PERIODS
  # ==========================================================================

  if (!is.data.frame(periods) || !nrow(periods)) {
    stop("`periods` must be a non-empty data.frame.", call. = FALSE)
  }

  required_period_columns <- c("period", "lag_start", "lag_end")

  if (!all(required_period_columns %in% names(periods))) {
    stop(
      "`periods` must contain `period`, `lag_start`, and `lag_end`.",
      call. = FALSE
    )
  }

  periods <- periods

  if (anyNA(periods$period) ||
      any(!nzchar(as.character(periods$period))) ||
      anyDuplicated(as.character(periods$period))) {
    stop(
      "`periods$period` must contain unique, non-missing, non-empty labels.",
      call. = FALSE
    )
  }

  periods$period <- as.character(periods$period)

  for (column in c("lag_start", "lag_end")) {
    x <- periods[[column]]

    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
        any(x < 0) || any(x != as.integer(x))) {
      stop(
        "`periods$", column,
        "` must contain non-negative finite integers.",
        call. = FALSE
      )
    }

    periods[[column]] <- as.integer(x)
  }

  if (any(periods$lag_start > periods$lag_end)) {
    stop(
      "Every period must satisfy `lag_start <= lag_end`.",
      call. = FALSE
    )
  }

  # Because every generated exposure combination is assigned to every period,
  # overlapping periods would assign the same exposure-lag position more than
  # once. Reject this here rather than relying on simulate_scenarios() to detect
  # the conflict later.
  period_lags <- lapply(
    seq_len(nrow(periods)),
    function(i) seq.int(periods$lag_start[i], periods$lag_end[i])
  )

  if (length(period_lags) > 1L) {
    for (i in seq_len(length(period_lags) - 1L)) {
      for (j in seq.int(i + 1L, length(period_lags))) {
        overlap_lags <- intersect(period_lags[[i]], period_lags[[j]])

        if (length(overlap_lags)) {
          stop(
            "Lag periods '", periods$period[i], "' and '",
            periods$period[j], "' overlap at lag(s): ",
            paste(overlap_lags, collapse = ", "),
            ". `simulate_ranges()` assigns every generated exposure ",
            "combination to every period, so periods must not overlap.",
            call. = FALSE
          )
        }
      }
    }
  }

  covered_lags <- sort(unique(unlist(period_lags, use.names = FALSE)))
  lag_min <- min(covered_lags)
  lag_max_defined <- max(covered_lags)

  expected_within_defined_range <- seq.int(lag_min, lag_max_defined)
  gap_lags <- setdiff(expected_within_defined_range, covered_lags)

  # ==========================================================================
  # VARYING EXPOSURES
  # ==========================================================================

  if (!is.list(vary) || !length(vary)) {
    stop(
      "`vary` must be a non-empty named list.",
      call. = FALSE
    )
  }

  var_names <- names(vary)

  if (is.null(var_names) || anyNA(var_names) ||
      any(!nzchar(var_names)) || anyDuplicated(var_names)) {
    stop(
      "`vary` must have unique, non-empty variable names.",
      call. = FALSE
    )
  }

  for (nm in var_names) {
    values <- vary[[nm]]

    if (!is.numeric(values) || !length(values) ||
        !is.null(dim(values)) || anyNA(values) ||
        any(!is.finite(values))) {
      stop(
        "`vary[['", nm,
        "']]` must contain at least one finite numeric value.",
        call. = FALSE
      )
    }

    values <- as.numeric(values)

    if (anyDuplicated(values)) {
      stop(
        "`vary[['", nm,
        "']]` contains duplicated values. Use unique exposure values to avoid ",
        "duplicated scenario combinations.",
        call. = FALSE
      )
    }

    vary[[nm]] <- values
  }

  # ==========================================================================
  # FIXED EXPOSURES
  # ==========================================================================

  if (!is.list(fixed)) {
    stop("`fixed` must be a named list.", call. = FALSE)
  }

  if (length(fixed)) {
    fixed_names <- names(fixed)

    if (is.null(fixed_names) || anyNA(fixed_names) ||
        any(!nzchar(fixed_names)) || anyDuplicated(fixed_names)) {
      stop(
        "`fixed` must have unique, non-empty variable names.",
        call. = FALSE
      )
    }

    for (nm in fixed_names) {
      value <- fixed[[nm]]

      if (!is.numeric(value) || length(value) != 1L ||
          is.na(value) || !is.finite(value)) {
        stop(
          "`fixed[['", nm,
          "']]` must contain exactly one finite numeric value.",
          call. = FALSE
        )
      }

      fixed[[nm]] <- as.numeric(value)
    }
  }

  overlap_variables <- intersect(names(vary), names(fixed))

  if (length(overlap_variables)) {
    stop(
      "Variables cannot appear in both `vary` and `fixed`: ",
      paste(overlap_variables, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # ==========================================================================
  # EXPOSURE COMBINATIONS
  # ==========================================================================

  if (identical(scenario_type, "grid")) {
    grid <- expand.grid(
      vary,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    lengths <- vapply(vary, length, integer(1))

    if (length(unique(lengths)) != 1L) {
      stop(
        "For `scenario_type = 'paired'`, all vectors in `vary` must have ",
        "the same length.",
        call. = FALSE
      )
    }

    grid <- as.data.frame(
      vary,
      stringsAsFactors = FALSE,
      optional = TRUE
    )
  }

  n_scenarios <- nrow(grid)

  if (n_scenarios < 1L) {
    stop(
      "No scenarios could be generated from `vary`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(grid)) {
    stop(
      "The requested exposure settings generated duplicated scenario ",
      "combinations.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # SCENARIO VARIABLE AND NAMES
  # ==========================================================================

  if (is.null(scenario_var)) {
    scenario_var <- var_names[1L]
  }

  if (!is.character(scenario_var) || length(scenario_var) != 1L ||
      is.na(scenario_var) || !nzchar(scenario_var)) {
    stop(
      "`scenario_var` must be NULL or one non-empty character value.",
      call. = FALSE
    )
  }

  if (!scenario_var %in% var_names) {
    stop(
      "`scenario_var` must be one of the variable names in `vary`.",
      call. = FALSE
    )
  }

  if (is.null(scenario_names)) {
    prefix <- toupper(scenario_var)
    scenario_names <- paste0(
      prefix,
      "_s",
      seq_len(n_scenarios)
    )
  } else {
    if (!is.character(scenario_names) ||
        length(scenario_names) != n_scenarios ||
        anyNA(scenario_names) ||
        any(!nzchar(scenario_names))) {
      stop(
        "`scenario_names` must contain exactly one non-empty character name ",
        "for each generated scenario.",
        call. = FALSE
      )
    }

    if (anyDuplicated(scenario_names)) {
      stop(
        "`scenario_names` must be unique.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # SCENARIO DEFINITIONS
  # ==========================================================================

  period_labels <- periods$period
  scenarios_out <- vector("list", n_scenarios)

  for (i in seq_len(n_scenarios)) {
    scenario_values <- as.list(
      grid[i, , drop = FALSE]
    )

    # Convert one-row data-frame cells explicitly to numeric scalars.
    scenario_values <- lapply(
      scenario_values,
      function(x) as.numeric(x[[1L]])
    )

    complete_period_block <- c(
      scenario_values,
      fixed
    )

    scenario_definition <- rep(
      list(complete_period_block),
      length(period_labels)
    )
    names(scenario_definition) <- period_labels

    scenarios_out[[i]] <- scenario_definition
  }

  names(scenarios_out) <- scenario_names

  # ==========================================================================
  # SCENARIO METADATA
  # ==========================================================================

  info_df <- grid

  if (length(fixed)) {
    for (nm in names(fixed)) {
      info_df[[nm]] <- rep(
        as.numeric(fixed[[nm]]),
        n_scenarios
      )
    }
  }

  info_df$scenario <- scenario_names

  info_df <- info_df[
    ,
    c(
      "scenario",
      setdiff(names(info_df), "scenario")
    ),
    drop = FALSE
  ]

  rownames(info_df) <- NULL

  # ==========================================================================
  # RETURN
  # ==========================================================================

  out <- list(
    scenarios = scenarios_out,
    periods = periods,
    info = info_df
  )

  attr(out, "scenario_type") <- scenario_type
  attr(out, "scenario_var") <- scenario_var
  attr(out, "period_lag_min") <- lag_min
  attr(out, "period_lag_max") <- lag_max_defined
  attr(out, "period_starts_at_zero") <- identical(lag_min, 0L)
  attr(out, "period_has_internal_gaps") <- length(gap_lags) > 0L
  attr(out, "period_gap_lags") <- as.integer(gap_lags)

  class(out) <- c(
    "epiexposure_ranges",
    "list"
  )

  out
}
