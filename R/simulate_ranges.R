#' Generate DLNM simulation scenarios (scenario or profile mode)
#'
#' Creates structured input for `simulate_scenarios()`, supporting both
#' discrete scenarios and chronological exposure profiles.
#'
#' In `mode = "profile"`, supplied vectors must be chronological exposure
#' profiles ordered from the earliest observation to the most recent
#' observation. Internal lag mapping is handled by downstream functions.
#'
#' @param periods Output from define_periods(). Must contain columns
#'   `period`, `lag_start`, and `lag_end`.
#' @param vary Named list of variables to vary.
#' @param fixed Named list of variables to keep fixed.
#' @param scenario_type Character. `"grid"` or `"paired"`.
#' @param mode Character. `"scenario"` or `"profile"`.
#' @param scenario_names Optional scenario names.
#' @param scenario_var Optional variable used to define automatic scenario names.
#'
#' @return A list containing:
#'   - scenarios
#'   - periods
#'   - info
#'
#' @export
simulate_ranges <- function(
    periods,
    vary = list(),
    fixed = list(),
    scenario_type = c("grid", "paired"),
    mode = c("scenario", "profile"),
    scenario_names = NULL,
    scenario_var = NULL
) {

  scenario_type <- match.arg(scenario_type)
  mode <- match.arg(mode)

  if (!is.data.frame(periods)) {
    stop("`periods` must be a data.frame.")
  }

  req_cols <- c("period", "lag_start", "lag_end")
  if (!all(req_cols %in% names(periods))) {
    stop("`periods` must contain columns: period, lag_start, lag_end.")
  }

  if (anyDuplicated(periods$period)) {
    stop("`periods$period` must contain unique labels.")
  }

  if (!is.list(vary) || length(vary) == 0) {
    stop("`vary` must be a non-empty named list.")
  }

  var_names <- names(vary)
  if (is.null(var_names) || anyNA(var_names) || any(var_names == "")) {
    stop("`vary` must be a named list.")
  }

  for (nm in var_names) {
    if (length(vary[[nm]]) == 0) {
      stop("Variable '", nm, "' contains no values.")
    }
  }

  if (!is.list(fixed)) {
    stop("`fixed` must be a named list.")
  }

  overlap <- intersect(names(vary), names(fixed))
  if (length(overlap)) {
    stop(
      "Variables cannot appear in both `vary` and `fixed`: ",
      paste(overlap, collapse = ", "), "."
    )
  }

  if (mode == "scenario") {

    if (scenario_type == "grid") {
      grid <- expand.grid(
        vary,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      lens <- vapply(vary, length, integer(1))
      if (length(unique(lens)) != 1L) {
        stop(
          "For scenario_type = 'paired', all vectors in `vary` ",
          "must have the same length."
        )
      }
      grid <- as.data.frame(vary, stringsAsFactors = FALSE)
    }

    n_scen <- nrow(grid)

    if (is.null(scenario_var)) {
      scenario_var <- var_names[1]
    }

    if (!scenario_var %in% var_names) {
      stop("`scenario_var` must be one of the names in `vary`.")
    }

    if (is.null(scenario_names)) {
      prefix <- toupper(scenario_var)
      scenario_names <- paste0(prefix, "_s", seq_len(n_scen))
    }

    if (length(scenario_names) != n_scen) {
      stop("Length of `scenario_names` must equal number of scenarios.")
    }

    if (anyNA(scenario_names) || any(scenario_names == "")) {
      stop("`scenario_names` must not contain NA or empty values.")
    }

    if (anyDuplicated(scenario_names)) {
      stop("`scenario_names` must be unique.")
    }

    scenarios_out <- vector("list", n_scen)

    for (i in seq_len(n_scen)) {
      scen_row <- as.list(grid[i, , drop = FALSE])
      scen_i <- list()

      for (p in periods$period) {
        scen_period <- c(scen_row, fixed)
        scen_i[[p]] <- scen_period
      }

      scenarios_out[[i]] <- scen_i
    }

    names(scenarios_out) <- scenario_names

    info_df <- grid
    info_df$scenario <- scenario_names
    info_df <- info_df[, c("scenario", setdiff(names(info_df), "scenario")), drop = FALSE]

    return(list(
      scenarios = scenarios_out,
      periods = periods,
      info = info_df
    ))
  }

  if (mode == "profile") {

    if (length(vary) != 1L) {
      stop("For mode = 'profile', `vary` must contain exactly one variable.")
    }

    var <- var_names[1]
    x_vals <- vary[[1]]

    if (!is.numeric(x_vals) || anyNA(x_vals) || any(!is.finite(x_vals))) {
      stop("Profile values must be finite numeric values.")
    }

    if (is.null(scenario_var)) {
      scenario_var <- var
    }

    if (is.null(scenario_names)) {
      scenario_name <- paste0(toupper(scenario_var), "_profile")
    } else {
      scenario_name <- scenario_names[1]
    }

    if (is.na(scenario_name) || scenario_name == "") {
      stop("Invalid `scenario_names` for profile mode.")
    }

    scen <- list()

    for (p in periods$period) {
      scen_period <- c(
        setNames(list(as.numeric(x_vals)), var),
        fixed
      )
      scen[[p]] <- scen_period
    }

    scenarios_out <- list(scen)
    names(scenarios_out) <- scenario_name

    info_df <- data.frame(
      scenario = scenario_name,
      variable = var,
      mode = "profile",
      profile_order = "chronological",
      stringsAsFactors = FALSE
    )

    return(list(
      scenarios = scenarios_out,
      periods = periods,
      info = info_df
    ))
  }
}
