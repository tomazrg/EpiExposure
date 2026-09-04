#' Generate DLNM simulation scenarios across lag periods
#'
#' Creates structured scenario input for `simulate_scenarios()` by combining
#' user-defined exposure levels and assigning each resulting combination as a
#' constant condition within every supplied lag period.
#'
#' `periods` is interpreted on the retrospective lag scale. For example,
#' `lag_start = 0` and `lag_end = 20` represent lags 0 through 20, i.e. the
#' portion of the exposure history closest to the outcome assessment.
#'
#' The function generates scenario definitions only. Conversion from lag-period
#' definitions to complete chronological exposure histories is performed later
#' by `simulate_scenarios()` before prediction with `predict_outcome()`.
#'
#' @param periods Output from `define_periods()`. Must be a data.frame containing
#'   `period`, `lag_start`, and `lag_end`.
#' @param vary Named list of exposure variables and the finite numeric values to
#'   evaluate. Under `scenario_type = "grid"`, all combinations are generated.
#'   Under `scenario_type = "paired"`, values are matched by position and all
#'   vectors must have the same length.
#' @param fixed Named list of exposure variables to keep fixed. Each fixed
#'   variable must contain exactly one finite numeric value. Use `list()` when
#'   no fixed variables are required.
#' @param scenario_type Character. Either `"grid"` or `"paired"`.
#' @param scenario_names Optional character vector with one unique name for each
#'   generated scenario.
#' @param scenario_var Optional name of a variable in `vary` used to construct
#'   automatic scenario names. If `NULL`, the first variable in `vary` is used.
#'
#' @return A list containing:
#'   - `scenarios`: named scenario definitions indexed by lag period;
#'   - `periods`: the supplied lag-period definitions;
#'   - `info`: a data.frame describing the exposure combination in each scenario.
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

  # ---------------------------------------------------------------------------
  # Validate lag periods
  # ---------------------------------------------------------------------------
  if (!is.data.frame(periods) || nrow(periods) == 0L) {
    stop("`periods` must be a non-empty data.frame.")
  }

  required_period_columns <- c("period", "lag_start", "lag_end")
  if (!all(required_period_columns %in% names(periods))) {
    stop("`periods` must contain columns: period, lag_start, lag_end.")
  }

  if (anyNA(periods$period) || any(as.character(periods$period) == "")) {
    stop("`periods$period` must contain non-missing, non-empty labels.")
  }

  if (anyDuplicated(periods$period)) {
    stop("`periods$period` must contain unique labels.")
  }

  for (column in c("lag_start", "lag_end")) {
    x <- periods[[column]]
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
        any(x < 0) || any(x != as.integer(x))) {
      stop("`periods$", column, "` must contain non-negative finite integers.")
    }
  }

  if (any(periods$lag_start > periods$lag_end)) {
    stop("Every period must satisfy `lag_start <= lag_end`.")
  }

  # ---------------------------------------------------------------------------
  # Validate varying exposure values
  # ---------------------------------------------------------------------------
  if (!is.list(vary) || length(vary) == 0L) {
    stop("`vary` must be a non-empty named list.")
  }

  var_names <- names(vary)
  if (is.null(var_names) || anyNA(var_names) || any(var_names == "") ||
      anyDuplicated(var_names)) {
    stop("`vary` must have unique, non-empty variable names.")
  }

  for (nm in var_names) {
    values <- vary[[nm]]
    if (!is.numeric(values) || length(values) == 0L || anyNA(values) ||
        any(!is.finite(values))) {
      stop(
        "`vary[['", nm,
        "']]` must contain at least one finite numeric value."
      )
    }
    vary[[nm]] <- as.numeric(values)
  }

  # ---------------------------------------------------------------------------
  # Validate fixed exposure values
  # ---------------------------------------------------------------------------
  if (!is.list(fixed)) {
    stop("`fixed` must be a named list.")
  }

  if (length(fixed) > 0L) {
    fixed_names <- names(fixed)
    if (is.null(fixed_names) || anyNA(fixed_names) || any(fixed_names == "") ||
        anyDuplicated(fixed_names)) {
      stop("`fixed` must have unique, non-empty variable names.")
    }

    for (nm in fixed_names) {
      value <- fixed[[nm]]
      if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
          !is.finite(value)) {
        stop(
          "`fixed[['", nm,
          "']]` must contain exactly one finite numeric value."
        )
      }
      fixed[[nm]] <- as.numeric(value)
    }
  }

  overlap <- intersect(names(vary), names(fixed))
  if (length(overlap)) {
    stop(
      "Variables cannot appear in both `vary` and `fixed`: ",
      paste(overlap, collapse = ", "), "."
    )
  }

  # ---------------------------------------------------------------------------
  # Generate exposure combinations
  # ---------------------------------------------------------------------------
  if (scenario_type == "grid") {
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
        "the same length."
      )
    }

    grid <- as.data.frame(vary, stringsAsFactors = FALSE)
  }

  n_scenarios <- nrow(grid)
  if (n_scenarios < 1L) {
    stop("No scenarios could be generated from `vary`.")
  }

  # ---------------------------------------------------------------------------
  # Scenario names
  # ---------------------------------------------------------------------------
  if (is.null(scenario_var)) {
    scenario_var <- var_names[1L]
  }

  if (!is.character(scenario_var) || length(scenario_var) != 1L ||
      is.na(scenario_var) || !nzchar(scenario_var)) {
    stop("`scenario_var` must be NULL or one non-empty character value.")
  }

  if (!scenario_var %in% var_names) {
    stop("`scenario_var` must be one of the names in `vary`.")
  }

  if (is.null(scenario_names)) {
    prefix <- toupper(scenario_var)
    scenario_names <- paste0(prefix, "_s", seq_len(n_scenarios))
  } else {
    if (!is.character(scenario_names) || length(scenario_names) != n_scenarios ||
        anyNA(scenario_names) || any(!nzchar(scenario_names))) {
      stop(
        "`scenario_names` must be a character vector with one non-empty name ",
        "for each generated scenario."
      )
    }
  }

  if (anyDuplicated(scenario_names)) {
    stop("`scenario_names` must be unique.")
  }

  # ---------------------------------------------------------------------------
  # Build piecewise-constant scenario definitions
  # ---------------------------------------------------------------------------
  period_labels <- as.character(periods$period)
  scenarios_out <- vector("list", n_scenarios)

  for (i in seq_len(n_scenarios)) {
    scenario_values <- as.list(grid[i, , drop = FALSE])
    scenario_definition <- vector("list", length(period_labels))
    names(scenario_definition) <- period_labels

    for (j in seq_along(period_labels)) {
      scenario_definition[[j]] <- c(scenario_values, fixed)
    }

    scenarios_out[[i]] <- scenario_definition
  }

  names(scenarios_out) <- scenario_names

  # ---------------------------------------------------------------------------
  # Scenario metadata
  # ---------------------------------------------------------------------------
  info_df <- grid
  info_df$scenario <- scenario_names
  info_df <- info_df[
    , c("scenario", setdiff(names(info_df), "scenario")),
    drop = FALSE
  ]
  rownames(info_df) <- NULL

  list(
    scenarios = scenarios_out,
    periods = periods,
    info = info_df
  )
}
