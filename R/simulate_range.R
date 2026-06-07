#' Generate DLNM simulation scenarios (scenario or profile mode)
#'
#' Creates structured input for simulate_scenarios(), supporting both
#' discrete scenarios and continuous profiles (ranges).
#'
#' @param lag_windows Output from define_lag_windows()
#' @param vary Named list of variables to vary (vectors of values)
#' @param fixed Named list of fixed values
#' @param scenario_type "grid" (all combinations) or "paired"
#' @param mode "scenario" (default) or "profile"
#' @param scenario_names Optional custom names
#' @param scenario_var Optional variable used to define scenario naming (recommended)
#'
#' @return A structured list with:
#'   - scenarios: named list of scenarios
#'   - lag_windows: lag window table used to define timing
#'   - info: data.frame describing scenario values (NULL for profile mode if not needed)
#'
#' @export
simulate_range <- function(
    lag_windows,
    vary = list(),
    fixed = list(),
    scenario_type = c("grid", "paired"),
    mode = c("scenario", "profile"),
    scenario_names = NULL,
    scenario_var = NULL
) {

  scenario_type <- match.arg(scenario_type)
  mode <- match.arg(mode)

  # ------------------------------------------------------------
  # validations
  # ------------------------------------------------------------
  if (!is.data.frame(lag_windows)) {
    stop("`lag_windows` must be a data.frame.")
  }

  req_cols <- c("window_id", "lag_start", "lag_end")
  if (!all(req_cols %in% names(lag_windows))) {
    stop("`lag_windows` must contain columns: window_id, lag_start, lag_end.")
  }

  if (!is.list(vary) || length(vary) == 0) {
    stop("`vary` must be a non-empty named list.")
  }

  var_names <- names(vary)
  if (is.null(var_names) || any(var_names == "")) {
    stop("`vary` must be a named list.")
  }

  if (!is.list(fixed)) {
    stop("`fixed` must be a named list.")
  }

  # ============================================================
  # ✅ SCENARIO MODE
  # ============================================================
  if (mode == "scenario") {

    # ---- build combinations ----
    if (scenario_type == "grid") {

      grid <- expand.grid(
        vary,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

    } else {

      lens <- sapply(vary, length)
      if (length(unique(lens)) != 1) {
        stop("For scenario_type = 'paired', all vectors in `vary` must have the same length.")
      }

      grid <- as.data.frame(vary, stringsAsFactors = FALSE)
    }

    n_scen <- nrow(grid)

    # ---- define naming variable ----
    if (is.null(scenario_var)) {
      scenario_var <- var_names[1]
    }

    if (!scenario_var %in% var_names) {
      stop("`scenario_var` must be one of the names in `vary`.")
    }

    # ---- scenario names ----
    if (is.null(scenario_names)) {
      prefix <- toupper(scenario_var)
      scenario_names <- paste0(prefix, "_s", seq_len(n_scen))
    }

    if (length(scenario_names) != n_scen) {
      stop("Length of `scenario_names` must equal number of scenarios.")
    }

    # ---- build scenarios ----
    scenarios_out <- vector("list", n_scen)

    for (i in seq_len(n_scen)) {

      scen_row <- as.list(grid[i, , drop = FALSE])
      scen_i <- list()

      for (w in lag_windows$window_id) {

        scen_window <- scen_row

        if (length(fixed) > 0) {
          scen_window <- c(scen_window, fixed)
        }

        scen_i[[w]] <- scen_window
      }

      scenarios_out[[i]] <- scen_i
    }

    names(scenarios_out) <- scenario_names

    info_df <- grid
    info_df$scenario <- scenario_names
    info_df <- info_df[, c("scenario", setdiff(names(info_df), "scenario")), drop = FALSE]

    return(list(
      scenarios = scenarios_out,
      lag_windows = lag_windows,
      info = info_df
    ))
  }

  # ============================================================
  # ✅ PROFILE MODE
  # ============================================================
  if (mode == "profile") {

    if (length(vary) != 1) {
      stop("For mode = 'profile', `vary` must contain exactly one variable.")
    }

    var <- var_names[1]
    x_vals <- vary[[1]]

    if (is.null(scenario_var)) {
      scenario_var <- var
    }

    if (is.null(scenario_names)) {
      scenario_name <- paste0(toupper(scenario_var), "_profile")
    } else {
      scenario_name <- scenario_names[1]
    }

    scen <- list()

    for (w in lag_windows$window_id) {

      scen_window <- c(
        setNames(list(x_vals), var),
        fixed
      )

      scen[[w]] <- scen_window
    }

    scenarios_out <- list(scen)
    names(scenarios_out) <- scenario_name

    info_df <- data.frame(
      scenario = scenario_name,
      variable = var,
      mode = "profile",
      stringsAsFactors = FALSE
    )

    return(list(
      scenarios = scenarios_out,
      lag_windows = lag_windows,
      info = info_df
    ))
  }
}
