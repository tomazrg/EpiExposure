#' Simulate epidemiological exposure-profile scenarios
#'
#' Builds complete chronological exposure profiles from user-defined lag-period
#' scenarios and predicts their population-level expected outcomes with
#' `predict_outcomes()`.
#'
#' `simulate_scenarios()` is intentionally a thin scenario-construction wrapper.
#' It does not contain engine-specific prediction code, does not refit or
#' redefine DLNM bases, and does not summarize parameter draws independently of
#' `predict_outcomes()`. All model-specific prediction, link handling, and
#' uncertainty propagation are delegated to the harmonized prediction layer.
#'
#' @param fit Fitted model returned by the current `fit_epidlnm()`.
#' @param scenarios Scenario definitions. Two input forms are supported:
#'
#'   \enumerate{
#'     \item A named list in which each top-level element is one scenario.
#'       Within a scenario, elements are named lag periods, and each period is a
#'       named list of fitted exposure variables.
#'     \item A structured object containing a `$scenarios` component and,
#'       optionally, `$periods` and `$info` components, such as the object
#'       produced by EpiExposure scenario-generation utilities.
#'   }
#'
#'   For a period/exposure entry, a scalar applies one value uniformly to every
#'   lag in that period. A numeric vector of length greater than one represents
#'   multiple matched **scenario points**, not a lag-varying within-period
#'   trajectory. All varying vectors within the same scenario must have the same
#'   length; scalar entries are recycled across those points.
#'
#'   If a genuinely lag-varying chronological exposure profile is required,
#'   supply that profile directly to `predict_outcomes(profiles = ...)` instead
#'   of encoding it as a period vector here.
#'
#' @param data Optional long-format exposure data used only when a scenario
#'   leaves one or more fitted exposure-lag positions unspecified and no
#'   corresponding value was supplied in `ref_vals`. In that case, the median
#'   of the required fitted exposure in `data` is used as the background value.
#'   If every scenario completely specifies every fitted exposure over its full
#'   lag window, `data` is not required. `data` is never used to reconstruct,
#'   refit, or redefine a cross-basis.
#' @param periods Optional data frame containing `period`, `lag_start`, and
#'   `lag_end`. If `scenarios` is a structured object with an embedded
#'   `$periods`, `periods` may be omitted. If both are supplied, their required
#'   columns must agree exactly; conflicting definitions are rejected.
#' @param ref_vals Optional named list of finite numeric background values.
#'   Values are used **only** for fitted exposure-lag positions left unspecified
#'   by a scenario. The list may contain all fitted exposures or only the
#'   exposures for which background filling is needed.
#'
#'   If an unassigned position requires a background value and that exposure is
#'   absent from `ref_vals`, its median is derived from `data`. If neither source
#'   is available, the function stops with an explicit error naming the exposure
#'   that still requires a background value.
#'
#'   Complete scenarios therefore do not need `ref_vals` or `data`. For example,
#'   when `simulate_ranges()` copies every fitted exposure combination across
#'   non-overlapping periods that jointly cover the full fitted lag window,
#'   every exposure-lag position is already specified.
#'
#'   EpiExposure does not use `argvar$cen` or another spline-centering attribute
#'   as an implicit scenario background because DLNM centering and a complete
#'   epidemiological background profile are distinct concepts.
#' @param uncertainty Logical. If `FALSE`, predict from the harmonized central
#'   fixed/population parameter estimate. If `TRUE`, propagate joint
#'   fixed/population parameter uncertainty draw by draw through
#'   `predict_outcomes()`.
#' @param output Character. `"summary"` returns deterministic predictions when
#'   `uncertainty = FALSE`, or median, SD, and empirical uncertainty intervals
#'   when `uncertainty = TRUE`. `"samples"` returns one row per parameter draw
#'   and scenario point and requires `uncertainty = TRUE`.
#' @param n_samples Positive integer number of parameter draws used when
#'   `uncertainty = TRUE`. At least two are required.
#' @param probs Numeric vector of length two defining the empirical uncertainty
#'   interval when `uncertainty = TRUE` and `output = "summary"`. The default
#'   `c(0.025, 0.975)` gives a 95 percent interval.
#' @param seed Optional finite integer for reproducible parameter sampling.
#'   Random-number handling is delegated to `predict_outcomes()`, which restores
#'   the caller's global random-number state when it exits.
#' @param extrapolation Character. Behavior when a complete scenario profile
#'   contains exposure values outside the fitted cross-basis exposure range:
#'   `"warn"` (default), `"error"`, or `"allow"`. The fitted basis is never
#'   re-estimated from scenario values.
#'
#' @return A data frame containing scenario identifiers, optional scenario
#'   metadata, and population-level expected-response predictions.
#'
#'   `scenario` identifies the named scenario and `scenario_point` identifies
#'   matched points when a scenario contains varying vectors.
#'
#'   With `uncertainty = FALSE`, the result contains `prediction`.
#'
#'   With `uncertainty = TRUE` and `output = "summary"`, it contains
#'   `prediction`, `prediction_sd`, `prediction_lower`, and
#'   `prediction_upper`.
#'
#'   With `uncertainty = TRUE` and `output = "samples"`, it contains `sample`
#'   and `prediction`, with the same `sample` index referring to the same joint
#'   fixed/population parameter draw across all scenario points.
#'
#' @details
#' ## Scenario interpretation
#'
#' Each scenario is converted into one complete chronological exposure profile
#' for every fitted exposure. Period-specific assignments are written first.
#' Only positions that remain unassigned are filled from `ref_vals` or, when
#' needed, exposure medians derived from `data`.
#'
#' This means a complete scenario is predicted exactly as supplied. Background
#' values do not modify or recenter positions that the scenario already defines.
#'
#' For a fitted maximum lag \eqn{L}, a profile has \eqn{L + 1} elements ordered
#' chronologically:
#'
#' \deqn{(x_L, x_{L-1}, \ldots, x_1, x_0),}
#'
#' so the final profile element corresponds to lag 0.
#'
#' A scenario prediction therefore describes the expected outcome under the
#' **entire assembled exposure profile**, not the effect of a single lag in
#' isolation.
#'
#' ## Population-level expected response
#'
#' EpiExposure v1 reports scenario predictions from the fixed/population
#' component of the fitted model:
#'
#' \deqn{\eta = X\beta,}
#'
#' with fitted random effects set to zero. Random effects may have contributed
#' to model estimation, but no fitted group-specific random effect is inherited
#' by a scenario. With response-scale output, the target is
#'
#' \deqn{E(Y \mid X) = g^{-1}(X\beta).}
#'
#' This is not integration over the random-effect distribution and is not a
#' simulation of a future observed outcome.
#'
#' ## Uncertainty and scenario comparisons
#'
#' When `uncertainty = TRUE`, all scenario profiles are passed to
#' `predict_outcomes()` in a **single prediction call**. Consequently, parameter
#' draw \eqn{s} is applied jointly to every scenario:
#'
#' \deqn{\mu_j^{(s)} = g^{-1}(X_j\beta^{(s)}),}
#'
#' where \eqn{j} indexes scenario points. This preserves covariance among
#' scenario predictions and allows sample-wise differences between scenarios to
#' be calculated correctly downstream.
#'
#' The uncertainty distribution describes uncertainty in expected responses.
#' Residual, observation, process, dispersion, posterior-predictive, and
#' group-specific random-effect noise are not added.
#'
#' ## Overlapping periods
#'
#' `periods` may contain overlapping intervals, but within a single scenario the
#' same exposure cannot be assigned by two overlapping period blocks. Such a
#' conflict is rejected rather than allowing one block to silently overwrite the
#' other.
#'
#' @export
simulate_scenarios <- function(
    fit,
    scenarios,
    data = NULL,
    periods = NULL,
    ref_vals = NULL,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    probs = c(0.025, 0.975),
    seed = NULL,
    extrapolation = c("warn", "error", "allow")
) {

  # ARGUMENTS AND STRICT MODEL CONTRACT

  output <- match.arg(output)
  extrapolation <- match.arg(extrapolation)

  if (is.null(fit)) {
    stop("`fit` cannot be NULL.", call. = FALSE)
  }

  metadata <- .get_epiexposure_metadata(fit)

  if (!identical(metadata$prediction_level, "population")) {
    stop(
      "`simulate_scenarios()` supports population-level prediction only.",
      call. = FALSE
    )
  }

  if (!identical(metadata$prediction_estimand, "expected_response")) {
    stop(
      "`simulate_scenarios()` requires the fitted prediction estimand to be ",
      "'expected_response'.",
      call. = FALSE
    )
  }

  if (!is.logical(uncertainty) || length(uncertainty) != 1L ||
      is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!uncertainty && identical(output, "samples")) {
    stop(
      "`output = 'samples'` requires `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  n_samples <- .epix_validate_n_samples(n_samples)

  if (uncertainty && n_samples < 2L) {
    stop(
      "`n_samples` must be at least 2 when `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  if (uncertainty && identical(output, "summary")) {
    probs <- .epix_validate_probs(probs)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed != as.integer(seed)) {
      stop("`seed` must be NULL or one finite integer.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # NORMALIZE STRUCTURED SCENARIO INPUT

  scenario_info <- NULL
  embedded_periods <- NULL

  if (is.list(scenarios) && "scenarios" %in% names(scenarios)) {
    scenario_info <- scenarios$info %||% NULL
    embedded_periods <- scenarios$periods %||% NULL
    scenarios <- scenarios$scenarios
  }

  if (!is.list(scenarios) || !length(scenarios) ||
      is.null(names(scenarios)) || anyNA(names(scenarios)) ||
      any(!nzchar(names(scenarios))) || anyDuplicated(names(scenarios))) {
    stop(
      "`scenarios` must be a non-empty named list with unique scenario names.",
      call. = FALSE
    )
  }

  scenario_names <- names(scenarios)

  # PERIOD DEFINITIONS

  validate_periods <- function(x, label) {
    required <- c("period", "lag_start", "lag_end")

    if (!is.data.frame(x) || !nrow(x) ||
        !all(required %in% names(x))) {
      stop(
        "`", label, "` must be a non-empty data.frame containing `period`, ",
        "`lag_start`, and `lag_end`.",
        call. = FALSE
      )
    }

    out <- x

    if (anyNA(out$period) ||
        any(!nzchar(as.character(out$period))) ||
        anyDuplicated(as.character(out$period))) {
      stop(
        "`", label, "$period` must contain unique non-empty labels.",
        call. = FALSE
      )
    }

    out$period <- as.character(out$period)

    for (nm in c("lag_start", "lag_end")) {
      value <- out[[nm]]

      if (!is.numeric(value) || anyNA(value) ||
          any(!is.finite(value)) || any(value < 0) ||
          any(value != as.integer(value))) {
        stop(
          "`", label, "$", nm,
          "` must contain non-negative finite integers.",
          call. = FALSE
        )
      }

      out[[nm]] <- as.integer(value)
    }

    if (any(out$lag_start > out$lag_end)) {
      stop(
        "Every row of `", label,
        "` must satisfy `lag_start <= lag_end`.",
        call. = FALSE
      )
    }

    out
  }

  explicit_periods <- if (!is.null(periods)) {
    validate_periods(periods, "periods")
  } else {
    NULL
  }

  embedded_periods <- if (!is.null(embedded_periods)) {
    validate_periods(embedded_periods, "scenarios$periods")
  } else {
    NULL
  }

  if (is.null(explicit_periods) && is.null(embedded_periods)) {
    stop(
      "`periods` is required unless the structured `scenarios` object contains ",
      "a valid `$periods` component.",
      call. = FALSE
    )
  }

  if (!is.null(explicit_periods) && !is.null(embedded_periods)) {
    required <- c("period", "lag_start", "lag_end")

    explicit_required <- explicit_periods[, required, drop = FALSE]
    embedded_required <- embedded_periods[, required, drop = FALSE]

    same_periods <- identical(
      explicit_required,
      embedded_required
    )

    if (!same_periods) {
      stop(
        "`periods` conflicts with the `$periods` definition embedded in ",
        "`scenarios`. Supply only one definition or make the required columns ",
        "identical.",
        call. = FALSE
      )
    }
  }

  periods_use <- explicit_periods %||% embedded_periods

  # Fast lookup by period label while preserving user order elsewhere.
  period_row <- stats::setNames(
    seq_len(nrow(periods_use)),
    periods_use$period
  )

  # OPTIONAL SCENARIO METADATA


  reserved_info_names <- c(
    "profile",
    "scenario_point",
    "sample",
    "prediction",
    "prediction_sd",
    "prediction_lower",
    "prediction_upper"
  )

  if (!is.null(scenario_info)) {
    if (!is.data.frame(scenario_info)) {
      stop(
        "`scenarios$info` must be NULL or a data.frame.",
        call. = FALSE
      )
    }

    if (!"scenario" %in% names(scenario_info)) {
      stop(
        "`scenarios$info` must contain a `scenario` column.",
        call. = FALSE
      )
    }

    if (anyNA(scenario_info$scenario) ||
        any(!nzchar(as.character(scenario_info$scenario)))) {
      stop(
        "`scenarios$info$scenario` cannot contain missing or empty values.",
        call. = FALSE
      )
    }

    unknown_info <- setdiff(
      unique(as.character(scenario_info$scenario)),
      scenario_names
    )

    if (length(unknown_info)) {
      stop(
        "`scenarios$info` contains scenario name(s) absent from `scenarios`: ",
        paste(unknown_info, collapse = ", "), ".",
        call. = FALSE
      )
    }

    conflicts <- intersect(names(scenario_info), reserved_info_names)
    if (length(conflicts)) {
      stop(
        "`scenarios$info` uses reserved output column name(s): ",
        paste(conflicts, collapse = ", "), ".",
        call. = FALSE
      )
    }

    if ("profile_order" %in% names(scenario_info)) {
      profile_order <- as.character(scenario_info$profile_order)

      if (anyNA(profile_order) ||
          any(profile_order != "chronological")) {
        stop(
          "Every non-missing `profile_order` entry in `scenarios$info` must ",
          "equal 'chronological'.",
          call. = FALSE
        )
      }
    }
  }

  # OPTIONAL BACKGROUND SOURCES

  vars <- metadata$vars

  # `ref_vals` is validated here but is not required. Background values are
  # resolved lazily only if a scenario leaves exposure-lag positions unassigned.
  if (!is.null(ref_vals)) {
    if (!is.list(ref_vals) || !length(ref_vals) ||
        is.null(names(ref_vals)) || anyNA(names(ref_vals)) ||
        any(!nzchar(names(ref_vals))) ||
        anyDuplicated(names(ref_vals))) {
      stop(
        "`ref_vals` must be NULL or a named list with unique non-empty ",
        "exposure names.",
        call. = FALSE
      )
    }

    extra_ref <- setdiff(names(ref_vals), vars)

    if (length(extra_ref)) {
      stop(
        "`ref_vals` contains exposure variable(s) not fitted in the model: ",
        paste(extra_ref, collapse = ", "), ".",
        call. = FALSE
      )
    }

    valid_ref <- vapply(
      ref_vals,
      function(x) {
        is.numeric(x) &&
          length(x) == 1L &&
          !is.na(x) &&
          is.finite(x)
      },
      logical(1)
    )

    if (any(!valid_ref)) {
      stop(
        "Every supplied `ref_vals` element must be one finite numeric value.",
        call. = FALSE
      )
    }

    ref_vals <- lapply(ref_vals, as.numeric)
  }

  if (!is.null(data)) {
    if (!is.data.frame(data) || !nrow(data)) {
      stop(
        "`data` must be NULL or a non-empty data.frame.",
        call. = FALSE
      )
    }
  }

  background_cache <- stats::setNames(
    vector("list", length(vars)),
    vars
  )
  background_source <- stats::setNames(
    rep(NA_character_, length(vars)),
    vars
  )

  get_background_value <- function(variable) {

    if (!is.null(background_cache[[variable]])) {
      return(background_cache[[variable]])
    }

    if (!is.null(ref_vals) && variable %in% names(ref_vals)) {
      value <- as.numeric(ref_vals[[variable]])
      background_cache[[variable]] <<- value
      background_source[[variable]] <<- "ref_vals"
      return(value)
    }

    if (!is.null(data)) {
      if (!variable %in% names(data)) {
        stop(
          "Scenario construction requires a background value for fitted ",
          "exposure '", variable, "', but that variable is absent from `data` ",
          "and was not supplied in `ref_vals`.",
          call. = FALSE
        )
      }

      x <- data[[variable]]

      if (!is.numeric(x) || all(is.na(x)) ||
          any(!is.finite(x[!is.na(x)]))) {
        stop(
          "Scenario construction requires a background value for fitted ",
          "exposure '", variable, "', but `data[['", variable,
          "']]` does not contain usable finite numeric values.",
          call. = FALSE
        )
      }

      value <- as.numeric(stats::median(x, na.rm = TRUE))
      background_cache[[variable]] <<- value
      background_source[[variable]] <<- "data_median"
      return(value)
    }

    stop(
      "Scenario construction left one or more lag positions unspecified for ",
      "fitted exposure '", variable, "'. Supply `ref_vals[['", variable,
      "']]`, provide `data` so its median can be used, or define that exposure ",
      "over the complete fitted lag window.",
      call. = FALSE
    )
  }

  # SCENARIO STRUCTURE VALIDATION

  validate_scenario <- function(scenario, scenario_name) {

    # Empty scenario = complete reference/background profile.
    if (is.null(scenario)) {
      scenario <- list()
    }

    if (!is.list(scenario)) {
      stop(
        "Scenario '", scenario_name, "' must be a list of period blocks.",
        call. = FALSE
      )
    }

    if (!length(scenario)) {
      return(list(
        scenario = scenario,
        n_points = 1L
      ))
    }

    if (is.null(names(scenario)) || anyNA(names(scenario)) ||
        any(!nzchar(names(scenario))) ||
        anyDuplicated(names(scenario))) {
      stop(
        "Scenario '", scenario_name,
        "' must use unique non-empty period names.",
        call. = FALSE
      )
    }

    unknown_periods <- setdiff(names(scenario), periods_use$period)
    if (length(unknown_periods)) {
      stop(
        "Scenario '", scenario_name,
        "' refers to period(s) not found in `periods`: ",
        paste(unknown_periods, collapse = ", "), ".",
        call. = FALSE
      )
    }

    varying_lengths <- integer(0)

    for (period_name in names(scenario)) {
      block <- scenario[[period_name]]

      if (is.null(block)) {
        block <- list()
        scenario[[period_name]] <- block
      }

      if (!is.list(block)) {
        stop(
          "Scenario '", scenario_name, "', period '", period_name,
          "' must be a named list of fitted exposure values.",
          call. = FALSE
        )
      }

      if (!length(block)) {
        next
      }

      if (is.null(names(block)) || anyNA(names(block)) ||
          any(!nzchar(names(block))) ||
          anyDuplicated(names(block))) {
        stop(
          "Scenario '", scenario_name, "', period '", period_name,
          "' must use unique non-empty exposure names.",
          call. = FALSE
        )
      }

      unknown_vars <- setdiff(names(block), vars)
      if (length(unknown_vars)) {
        stop(
          "Scenario '", scenario_name, "', period '", period_name,
          "' contains exposure variable(s) not fitted in the model: ",
          paste(unknown_vars, collapse = ", "), ".",
          call. = FALSE
        )
      }

      for (variable in names(block)) {
        value <- block[[variable]]

        if (!is.numeric(value) || !length(value) ||
            !is.null(dim(value)) || anyNA(value) ||
            any(!is.finite(value))) {
          stop(
            "Scenario '", scenario_name, "', period '", period_name,
            "', variable '", variable,
            "' must be a finite numeric scalar or vector.",
            call. = FALSE
          )
        }

        if (length(value) > 1L) {
          varying_lengths <- c(varying_lengths, length(value))
        }
      }
    }

    unique_lengths <- unique(varying_lengths)

    if (length(unique_lengths) > 1L) {
      stop(
        "All varying vectors within scenario '", scenario_name,
        "' must have the same length. Received varying lengths: ",
        paste(sort(unique_lengths), collapse = ", "), ".",
        call. = FALSE
      )
    }

    n_points <- if (length(unique_lengths)) {
      as.integer(unique_lengths[1L])
    } else {
      1L
    }

    list(
      scenario = scenario,
      n_points = n_points
    )
  }

  validated_scenarios <- lapply(
    scenario_names,
    function(name) validate_scenario(scenarios[[name]], name)
  )
  names(validated_scenarios) <- scenario_names

  # SCENARIO INFO MAPPING

  scenario_info_for_point <- function(
    scenario_name,
    point_index,
    n_points
  ) {
    if (is.null(scenario_info)) {
      return(data.frame(row.names = 1L))
    }

    idx <- which(
      as.character(scenario_info$scenario) == scenario_name
    )

    if (!length(idx)) {
      return(data.frame(row.names = 1L))
    }

    current <- scenario_info[idx, , drop = FALSE]

    if (nrow(current) == 1L) {
      selected <- current[1L, , drop = FALSE]
    } else if (nrow(current) == n_points) {
      selected <- current[point_index, , drop = FALSE]
    } else {
      stop(
        "`scenarios$info` contains ", nrow(current),
        " row(s) for scenario '", scenario_name,
        "', but that scenario expands to ", n_points,
        " point(s). Supply either one metadata row or exactly one row per ",
        "scenario point.",
        call. = FALSE
      )
    }

    selected$scenario <- NULL
    rownames(selected) <- NULL
    selected
  }

  # BUILD COMPLETE CHRONOLOGICAL PROFILES

  all_profiles <- lapply(vars, function(x) list())
  names(all_profiles) <- vars

  key_rows <- list()
  profile_id <- 0L

  lag_to_index <- function(lags, n_profile) {
    # Chronological profile:
    # position n_profile = lag 0
    # position 1         = max_lag
    n_profile - as.integer(lags)
  }

  build_one_profile_set <- function(
    scenario,
    scenario_name,
    point_index,
    n_points
  ) {

    # Start undefined. Scenario blocks are written first. Background values are
    # requested only for positions that remain unassigned afterwards.
    profiles <- lapply(
      vars,
      function(variable) {
        rep(
          NA_real_,
          metadata$spec[[variable]]$max_lag + 1L
        )
      }
    )
    names(profiles) <- vars

    # Track assignments separately for each exposure so overlapping scenario
    # period blocks cannot silently overwrite one another.
    assigned <- lapply(
      vars,
      function(variable) {
        rep(FALSE, metadata$spec[[variable]]$max_lag + 1L)
      }
    )
    names(assigned) <- vars

    if (length(scenario)) {
      for (period_name in names(scenario)) {
        block <- scenario[[period_name]]

        if (!length(block)) {
          next
        }

        period_idx <- period_row[[period_name]]
        lag_start <- periods_use$lag_start[period_idx]
        lag_end <- periods_use$lag_end[period_idx]
        lags <- seq.int(lag_start, lag_end)

        for (variable in names(block)) {

          max_lag <- metadata$spec[[variable]]$max_lag

          if (lag_end > max_lag) {
            stop(
              "Scenario '", scenario_name, "', period '", period_name,
              "' assigns variable '", variable, "' through lag ", lag_end,
              ", but that exposure was fitted only through max_lag = ",
              max_lag, ".",
              call. = FALSE
            )
          }

          raw_value <- block[[variable]]

          value <- if (length(raw_value) == 1L) {
            as.numeric(raw_value)
          } else {
            as.numeric(raw_value[point_index])
          }

          n_profile <- max_lag + 1L
          idx <- lag_to_index(lags, n_profile)

          if (any(idx < 1L | idx > n_profile)) {
            stop(
              "Internal lag-to-profile mapping failed for scenario '",
              scenario_name, "', period '", period_name,
              "', variable '", variable, "'.",
              call. = FALSE
            )
          }

          if (any(assigned[[variable]][idx])) {
            previous_lags <- lags[assigned[[variable]][idx]]

            stop(
              "Scenario '", scenario_name, "' assigns variable '", variable,
              "' more than once at overlapping lag(s): ",
              paste(previous_lags, collapse = ", "),
              ". Revise the overlapping period blocks rather than relying on ",
              "silent overwriting.",
              call. = FALSE
            )
          }

          profiles[[variable]][idx] <- value
          assigned[[variable]][idx] <- TRUE
        }
      }
    }

    # Fill only positions that were not defined by the scenario.
    for (variable in vars) {
      missing_idx <- which(!assigned[[variable]])

      if (length(missing_idx)) {
        profiles[[variable]][missing_idx] <-
          get_background_value(variable)
      }

      if (any(!is.finite(profiles[[variable]]))) {
        stop(
          "Internal scenario construction produced non-finite values for ",
          "exposure '", variable, "'.",
          call. = FALSE
        )
      }
    }

    profiles
  }

  for (scenario_name in scenario_names) {
    current <- validated_scenarios[[scenario_name]]
    scenario <- current$scenario
    n_points <- current$n_points

    for (point_index in seq_len(n_points)) {
      profile_id <- profile_id + 1L

      current_profiles <- build_one_profile_set(
        scenario = scenario,
        scenario_name = scenario_name,
        point_index = point_index,
        n_points = n_points
      )

      for (variable in vars) {
        all_profiles[[variable]][[profile_id]] <-
          current_profiles[[variable]]
      }

      info_row <- scenario_info_for_point(
        scenario_name = scenario_name,
        point_index = point_index,
        n_points = n_points
      )

      key <- data.frame(
        profile = profile_id,
        scenario = scenario_name,
        scenario_point = point_index,
        stringsAsFactors = FALSE
      )

      if (ncol(info_row)) {
        key <- cbind(key, info_row)
      }

      key_rows[[profile_id]] <- key
    }
  }

  scenario_keys <- do.call(rbind, key_rows)
  rownames(scenario_keys) <- NULL

  if (profile_id < 1L) {
    stop(
      "No scenario profiles could be constructed.",
      call. = FALSE
    )
  }

  # ONE HARMONIZED PREDICTION CALL FOR ALL SCENARIOS

  prediction <- predict_outcomes(
    fit = fit,
    profiles = all_profiles,
    type = "response",
    uncertainty = uncertainty,
    output = output,
    n_samples = n_samples,
    probs = probs,
    seed = seed,
    extrapolation = extrapolation
  )

  # `predict_outcomes()` omits `profile` when only one prediction exists.
  if (!"profile" %in% names(prediction)) {
    if (profile_id != 1L || nrow(prediction) < 1L) {
      stop(
        "Internal profile-key mismatch after scenario prediction.",
        call. = FALSE
      )
    }
    prediction$profile <- 1L
  }

  if (anyNA(prediction$profile) ||
      any(!prediction$profile %in% scenario_keys$profile)) {
    stop(
      "Scenario prediction returned an invalid internal profile identifier.",
      call. = FALSE
    )
  }

  key_index <- match(
    prediction$profile,
    scenario_keys$profile
  )

  if (anyNA(key_index)) {
    stop(
      "Could not align scenario predictions with scenario metadata.",
      call. = FALSE
    )
  }

  output_keys <- scenario_keys[
    key_index,
    setdiff(names(scenario_keys), "profile"),
    drop = FALSE
  ]
  rownames(output_keys) <- NULL

  prediction_values <- prediction[
    ,
    setdiff(names(prediction), "profile"),
    drop = FALSE
  ]
  rownames(prediction_values) <- NULL

  result <- cbind(
    output_keys,
    prediction_values
  )
  rownames(result) <- NULL

  # OUTPUT CONTRACT METADATA

  background_values_used <- lapply(
    vars,
    function(variable) {
      value <- background_cache[[variable]]
      if (is.null(value)) NA_real_ else as.numeric(value)
    }
  )
  names(background_values_used) <- vars

  attr(result, "epiexposure_background_values") <- background_values_used
  attr(result, "epiexposure_background_sources") <- background_source
  attr(result, "epiexposure_background_used") <- any(
    !is.na(background_source)
  )
  attr(result, "epiexposure_periods") <- periods_use
  attr(result, "epiexposure_prediction_level") <- "population"
  attr(result, "epiexposure_prediction_estimand") <- "expected_response"
  attr(result, "epiexposure_prediction_type") <- "response"
  attr(result, "epiexposure_uncertainty") <- uncertainty
  attr(result, "epiexposure_scenario_count") <- length(scenario_names)
  attr(result, "epiexposure_scenario_point_count") <- profile_id

  if (uncertainty) {
    attr(result, "epiexposure_n_samples") <- n_samples

    if (identical(output, "summary")) {
      attr(result, "epiexposure_probs") <- probs
      attr(result, "epiexposure_summary_center") <- "median"
    } else {
      attr(result, "epiexposure_summary_center") <- "samples"
    }
  }

  result
}
