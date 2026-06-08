#' Simulate epidemiological DLNM scenarios
#'
#' Generates predicted outcomes for one or more user-defined epidemiological
#' scenarios, using exposure profiles assembled across lag periods.
#'
#' Scenarios can be supplied either as:
#' - a structured object returned by `simulate_range()`, or
#' - a named list of scenario definitions.
#'
#' When `uncertainty = TRUE`, uncertainty is propagated through
#' `predict_outcome()`. If `output = "summary"`, the central estimate is
#' computed as the median of simulated predictions, and interval limits are
#' obtained from empirical quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** in the summary output, the column `prediction` represents the
#' central estimate, computed as the median when uncertainty is propagated.
#'
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param scenarios Output from `simulate_range()` or a named list of scenarios.
#'   Recommended structured object:
#'   `list(scenarios = ..., info = ..., periods = ...)`.
#' @param wx_long Optional long-format weather data. Used only as fallback
#'   to derive default reference values when `ref_vals` are not supplied
#'   and the fitted model does not store centering information.
#' @param periods Optional period table. If `NULL`, the function will
#'   try to read it from `scenarios$periods`.
#' @param lag_max Optional maximum lag. Used only if `fit` does not store lag_max.
#' @param df_var Optional degrees of freedom (exposure). Ignored if `fit` contains
#'   `epiexposure_spec`. Kept only for fallback compatibility.
#' @param df_lag Optional degrees of freedom (lag). Ignored if `fit` contains
#'   `epiexposure_spec`. Kept only for fallback compatibility.
#' @param fun_var Optional basis function for exposure. Ignored if `fit` contains
#'   `epiexposure_spec`. Kept only for fallback compatibility.
#' @param fun_lag Optional basis function for lag. Ignored if `fit` contains
#'   `epiexposure_spec`. Kept only for fallback compatibility.
#' @param ref_vals Optional named list of reference values for each variable.
#'   If `NULL`, the function tries in order:
#'   1. median from `wx_long`
#'   2. centering value stored in `fit` spec (`argvar$cen`)
#' @param pop_level Logical. If `TRUE`, predictions exclude random effects
#'   where supported.
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Integer. Number of samples used for uncertainty quantification.
#'
#' @return A data.frame.
#'
#' If `uncertainty = FALSE`, returns one row per scenario point with column:
#' - `prediction`
#'
#' If `uncertainty = TRUE` and `output = "summary"`, returns one row per scenario
#' point with columns:
#' - `prediction` (median-based central estimate)
#' - `sd`
#' - `lower`
#' - `upper`
#'
#' If `uncertainty = TRUE` and `output = "samples"`, returns one row per sample
#' with columns:
#' - `sample`
#' - `prediction`
#'
#' Additional metadata columns from `scenarios$info` are preserved.
#'
#' @details
#' Uncertainty is propagated through `predict_outcome()`, which uses
#' model-consistent sampling:
#' - Bayesian models use posterior draws
#' - Frequentist models use simulation from the asymptotic coefficient distribution
#'
#' For summary outputs under uncertainty, the median is used instead of the mean
#' to provide a more robust central estimate under asymmetric predictive
#' distributions.
#'
#' @export
simulate_scenarios <- function(
    fit,
    scenarios,
    wx_long = NULL,
    periods = NULL,
    lag_max = NULL,
    df_var = NULL,
    df_lag = NULL,
    fun_var = NULL,
    fun_lag = NULL,
    ref_vals = NULL,
    pop_level = TRUE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  output <- match.arg(output)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  # ------------------------------------------------------------
  # validations
  # ------------------------------------------------------------
  if (is.null(fit)) stop("`fit` cannot be NULL.")
  if (!is.logical(pop_level) || length(pop_level) != 1L) {
    stop("`pop_level` must be TRUE or FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  fit_spec     <- attr(fit, "epiexposure_spec")
  fit_vars     <- attr(fit, "epiexposure_vars")
  dat_template <- attr(fit, "epiexposure_dat_template")

  if (is.null(dat_template) || !is.data.frame(dat_template) || nrow(dat_template) < 1) {
    stop("The fitted model does not contain `epiexposure_dat_template`. Refit with fit_epidlnm().")
  }

  if (is.null(fit_spec) || !is.list(fit_spec)) {
    stop("The fitted model does not contain `epiexposure_spec`. Refit with fit_epidlnm().")
  }

  # ------------------------------------------------------------
  # Handle structured scenarios
  # ------------------------------------------------------------
  scenario_info <- NULL
  periods_from_scenarios <- NULL

  if (is.list(scenarios) && "scenarios" %in% names(scenarios)) {
    scenario_info <- scenarios$info %||% NULL
    periods_from_scenarios <- scenarios$periods %||% NULL
    scenarios <- scenarios$scenarios
  }

  if (!is.list(scenarios) || is.null(names(scenarios))) {
    stop("`scenarios` must be a named list or a structured object containing `$scenarios`.")
  }

  # ------------------------------------------------------------
  # Resolve periods
  # ------------------------------------------------------------
  periods <- periods %||% periods_from_scenarios

  if (is.null(periods)) {
    stop("`periods` is missing. Provide it explicitly or include it inside the `scenarios` object.")
  }

  req_cols <- c("period", "lag_start", "lag_end")
  if (!all(req_cols %in% names(periods))) {
    stop("`periods` must contain columns: period, lag_start, lag_end.")
  }

  if (anyDuplicated(periods$period)) {
    stop("`periods$period` must contain unique labels.")
  }

  # ------------------------------------------------------------
  # Resolve variables
  # ------------------------------------------------------------
  if (!is.null(fit_vars) && length(fit_vars) > 0) {
    vars <- fit_vars
  } else {
    vars <- unique(unlist(lapply(scenarios, function(scen) {
      unique(unlist(lapply(scen, names)))
    })))
  }

  if (length(vars) == 0) {
    stop("No variables could be detected from fit metadata or scenarios.")
  }

  # ------------------------------------------------------------
  # Resolve per-variable lag length
  # ------------------------------------------------------------
  get_var_lagmax <- function(v) {
    if (!is.null(fit_spec[[v]]) && !is.null(fit_spec[[v]]$lag_max)) {
      return(as.integer(fit_spec[[v]]$lag_max))
    }
    if (!is.null(lag_max)) {
      return(as.integer(lag_max))
    }
    stop(
      "Could not determine lag_max for variable: ", v,
      ". The model must contain `epiexposure_spec[[var]]$lag_max`, or you must supply `lag_max`."
    )
  }

  var_lagmax <- setNames(lapply(vars, get_var_lagmax), vars)

  # ------------------------------------------------------------
  # Resolve reference values
  # ------------------------------------------------------------
  if (is.null(ref_vals)) {
    ref_vals <- vector("list", length(vars))
    names(ref_vals) <- vars

    for (v in vars) {

      if (!is.null(wx_long) && is.data.frame(wx_long) && v %in% names(wx_long)) {
        ref_vals[[v]] <- as.numeric(stats::median(wx_long[[v]], na.rm = TRUE))
        next
      }

      if (!is.null(fit_spec[[v]]) &&
          !is.null(fit_spec[[v]]$argvar) &&
          !is.null(fit_spec[[v]]$argvar$cen)) {
        ref_vals[[v]] <- as.numeric(fit_spec[[v]]$argvar$cen)
        next
      }

      stop(
        "Could not determine default reference value for variable '", v, "'. ",
        "Provide `ref_vals` explicitly or `wx_long`."
      )
    }
  }

  # ------------------------------------------------------------
  # Lag utility
  # ------------------------------------------------------------
  lag_to_idx <- function(lags, N) {
    idx <- N - lags
    idx[idx >= 1 & idx <= N]
  }

  # ------------------------------------------------------------
  # Build exposure profile for ONE variable under ONE scenario point
  # ------------------------------------------------------------
  build_profile <- function(var, scen_values) {

    N_var <- var_lagmax[[var]] + 1L
    x <- rep(ref_vals[[var]], N_var)

    for (i in seq_len(nrow(periods))) {

      p_id <- periods$period[i]

      if (!is.null(scen_values[[p_id]]) &&
          !is.null(scen_values[[p_id]][[var]])) {

        lags <- seq(periods$lag_start[i], periods$lag_end[i])
        idx <- lag_to_idx(lags, N_var)

        x[idx] <- scen_values[[p_id]][[var]]
      }
    }

    x
  }

  # ------------------------------------------------------------
  # Expand continuous scenario points
  # ------------------------------------------------------------
  expand_scenario_points <- function(scen_values) {

    lens <- unlist(lapply(scen_values, function(block) {
      sapply(block, length)
    }), use.names = FALSE)

    lens_gt1 <- unique(lens[lens > 1])

    if (length(lens_gt1) == 0) {
      return(list(scen_values))
    }

    if (length(lens_gt1) > 1) {
      stop("All varying scenario vectors must have the same length.")
    }

    n_pts <- lens_gt1[1]
    pts <- vector("list", n_pts)

    for (i in seq_len(n_pts)) {
      pts[[i]] <- lapply(scen_values, function(block) {
        lapply(block, function(val) {
          if (length(val) > 1) val[i] else val
        })
      })
    }

    pts
  }

  # ------------------------------------------------------------
  # Attach scenario metadata
  # ------------------------------------------------------------
  attach_scenario_info <- function(name, n_rows) {

    if (is.null(scenario_info)) {
      return(data.frame(
        scenario = rep(name, n_rows),
        stringsAsFactors = FALSE
      ))
    }

    df_info <- scenario_info[scenario_info$scenario == name, , drop = FALSE]

    if (nrow(df_info) == 0) {
      return(data.frame(
        scenario = rep(name, n_rows),
        stringsAsFactors = FALSE
      ))
    }

    if (nrow(df_info) == 1 && n_rows > 1) {
      df_info <- df_info[rep(1, n_rows), , drop = FALSE]
      rownames(df_info) <- NULL
      return(df_info)
    }

    if (nrow(df_info) != n_rows) {
      stop("scenario_info rows for scenario '", name, "' do not match the number of simulated points.")
    }

    df_info
  }

  # ------------------------------------------------------------
  # Apply scenarios using predict_outcome()
  # ------------------------------------------------------------
  re_mode <- if (isTRUE(pop_level)) "population" else "conditional"

  pred_output <- if (uncertainty) "samples" else output

  out_list <- purrr::imap(scenarios, function(scen, name) {

    scen_pts <- expand_scenario_points(scen)
    n_pts <- length(scen_pts)

    pred_list <- vector("list", n_pts)

    for (j in seq_len(n_pts)) {

      profiles_j <- lapply(vars, function(v) build_profile(v, scen_pts[[j]]))
      names(profiles_j) <- vars

      pred_j <- predict_outcome(
        fit = fit,
        profiles = profiles_j,
        re = re_mode,
        allow_new_levels = TRUE,
        type = "response",
        uncertainty = uncertainty,
        output = pred_output,
        n_samples = n_samples
      )

      meta_j <- attach_scenario_info(name, 1)

      if (!uncertainty) {

        meta_j$prediction <- pred_j$prediction[1]
        pred_list[[j]] <- meta_j

      } else if (output == "summary") {

        meta_j$prediction <- stats::median(pred_j$prediction, na.rm = TRUE)
        meta_j$sd <- stats::sd(pred_j$prediction, na.rm = TRUE)
        meta_j$lower <- safe_quantile(pred_j$prediction)[1]
        meta_j$upper <- safe_quantile(pred_j$prediction)[2]
        pred_list[[j]] <- meta_j

      } else {

        tmp <- meta_j[rep(1, nrow(pred_j)), , drop = FALSE]
        tmp$sample <- pred_j$sample
        tmp$prediction <- pred_j$prediction
        pred_list[[j]] <- tmp
      }
    }

    do.call(rbind, pred_list)
  })

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
