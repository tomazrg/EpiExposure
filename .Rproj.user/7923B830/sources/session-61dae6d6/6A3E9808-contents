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
#' @return Named list of scenarios
#'
#' @export
simulate_range <- function(lag_windows,
                           vary = list(),
                           fixed = list(),
                           scenario_type = c("grid", "paired"),
                           mode = c("scenario", "profile"),
                           scenario_names = NULL,
                           scenario_var = NULL) {
  
  scenario_type <- match.arg(scenario_type)
  mode <- match.arg(mode)
  
  stopifnot(length(vary) > 0)
  
  var_names <- names(vary)
  stopifnot(!is.null(var_names))
  
  # ============================================================
  # ✅ SCENARIO MODE (default)
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
      stopifnot(length(unique(lens)) == 1)
      
      grid <- as.data.frame(vary)
    }
    
    n_scen <- nrow(grid)
    
    # ---- define naming variable ----
    if (is.null(scenario_var)) {
      scenario_var <- var_names[1]
    }
    
    # ---- scenario names ----
    if (is.null(scenario_names)) {
      prefix <- toupper(scenario_var)
      scenario_names <- paste0(prefix, "_s", seq_len(n_scen))
    }
    
    stopifnot(length(scenario_names) == n_scen)
    
    # ---- build scenarios ----
    scenarios <- vector("list", n_scen)
    
    for (i in seq_len(n_scen)) {
      
      scen_row <- as.list(grid[i, ])
      scen_i <- list()
      
      for (w in lag_windows$window_id) {
        
        scen_window <- scen_row
        
        if (length(fixed) > 0) {
          scen_window <- c(scen_window, fixed)
        }
        
        scen_i[[w]] <- scen_window
      }
      
      scenarios[[i]] <- scen_i
    }
    
    names(scenarios) <- scenario_names
    return(scenarios)
  }
  
  # ============================================================
  # ✅ PROFILE MODE (continuous curve like article)
  # ============================================================
  if (mode == "profile") {
    
    # ---- only ONE variable can vary ----
    stopifnot(length(vary) == 1)
    
    var <- var_names[1]
    x_vals <- vary[[1]]
    
    # ---- scenario naming ----
    if (is.null(scenario_var)) {
      scenario_var <- var
    }
    
    if (is.null(scenario_names)) {
      scenario_name <- paste0(toupper(scenario_var), "_profile")
    } else {
      scenario_name <- scenario_names[1]
    }
    
    # ---- build single profile scenario ----
    scen <- list()
    
    for (w in lag_windows$window_id) {
      
      scen_window <- c(
        setNames(list(x_vals), var),
        fixed
      )
      
      scen[[w]] <- scen_window
    }
    
    scenarios <- list(scen)
    names(scenarios) <- scenario_name
    
    return(scenarios)
  }
}