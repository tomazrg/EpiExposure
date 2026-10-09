#' Simulate chronological exposure profiles
#'
#' Generates one or more complete exposure profiles in chronological order,
#' from the earliest/oldest observation to the most recent observation.
#'
#' The first profile value corresponds internally to the maximum retrospective
#' lag and the final profile value corresponds internally to lag 0:
#'
#' ```
#' time_0         = earliest/oldest observation = lag max_lag
#' time_max_lag   = most recent observation     = lag 0
#' ```
#'
#' `simulate_exposures()` simulates exposure profiles only. It does not draw
#' model coefficients and therefore has no `"summary"`/`"samples"` output
#' switch. The number of simulated exposure profiles is controlled exclusively
#' by `n`.
#'
#' @param max_lag Non-negative integer maximum retrospective lag. Every
#'   generated profile has exactly `max_lag + 1` chronological positions.
#' @param n Positive integer number of exposure profiles to simulate. Default
#'   is 1.
#' @param mode Character. `"profile"` generates the background process only;
#'   `"pattern"` generates the background and then applies one explicit
#'   chronological pattern.
#' @param pattern Character pattern used when `mode = "pattern"`:
#'
#'   - `"random_times"`: assign `fixed_value` at `n_times` distinct
#'     chronological positions sampled without replacement from `time_range`;
#'   - `"alternating_times"`: recycle `alternating_values` through
#'     `time_range` in chronological order;
#'   - `"block_times"`: assign `fixed_value` to the complete inclusive interval
#'     defined by `block_times`.
#' @param fixed_value Finite numeric scalar used by `"random_times"` and
#'   `"block_times"`.
#' @param n_times Positive integer number of selected chronological positions for
#'   `"random_times"`.
#' @param time_range Optional unique integer vector of chronological times used
#'   by `"random_times"` or `"alternating_times"`. Values must be between 0 and
#'   `max_lag`, inclusive. If `NULL`, every chronological position is eligible.
#'
#'   `time_range` is not used by `"block_times"`; supplying it with that pattern
#'   is an error rather than being silently ignored.
#' @param alternating_values Numeric vector containing at least two finite
#'   values. Values are recycled across the selected `time_range` in
#'   chronological order when `pattern = "alternating_times"`.
#' @param block_times Unique integer vector of length two defining the first and
#'   final chronological positions of the inclusive block when
#'   `pattern = "block_times"`.
#' @param background Named list describing the background process. Supported
#'   distributions are:
#'
#'   \describe{
#'     \item{`"normal"`}{`list(dist = "normal", mean = ..., sd = ...)`}
#'     \item{`"empirical"`}{`list(dist = "empirical", values = ...)`}
#'     \item{`"ar1"`}{`list(dist = "ar1", mean = ..., sd = ..., phi = ...)`}
#'     \item{`"fixed"`}{`list(dist = "fixed", value = ...)`}
#'   }
#'
#'   For the stationary Gaussian AR(1) generator, `sd` is the stationary
#'   marginal standard deviation and `phi` must lie strictly inside `(-1, 1)`.
#'
#'   Unknown background fields are rejected explicitly to catch misspelled
#'   parameter names.
#' @param bounds Optional finite numeric vector `c(lower, upper)`. Background
#'   values are clipped to these bounds before the pattern and values are clipped
#'   again after the pattern.
#'
#'   If `cumulative = TRUE`, `bounds` constrain the exposure **increments**
#'   before `cumsum()`; the final cumulative trajectory is not clipped again.
#' @param seed Optional finite integer seed controlling the complete set of `n`
#'   exposure simulations. The caller's random-number state is restored after
#'   the function returns.
#' @param cumulative Logical. If `TRUE`, applies `cumsum()` in chronological
#'   order after background generation, pattern assignment, and bounds.
#'
#' @return An object of class `"epiexposure_simulated_exposures"` with:
#'
#'   \describe{
#'     \item{`profiles`}{Canonical list containing exactly `n` simulated numeric
#'       chronological profiles, named `"simulation_1"`, ..., `"simulation_n"`.}
#'     \item{`profile`}{Backward-compatible numeric alias when `n = 1`;
#'       otherwise `NULL`.}
#'     \item{`simulation_data`}{Long data frame containing `simulation`,
#'       `position`, `time`, `lag`, and `profile`. This is a tabular view of the
#'       same simulated exposure profiles, not model-parameter samples.}
#'     \item{`meta`}{Common simulation metadata, including `max_lag`,
#'       chronological order, `n`, background specification, and RNG contract.}
#'     \item{`simulation_meta`}{List with simulation-specific pattern metadata,
#'       such as selected chronological positions and corresponding
#'       retrospective lags.}
#'   }
#'
#' @details
#' ## One simulation contract
#'
#' `n` is the only argument controlling how many exposure trajectories are
#' generated:
#'
#' ```
#' x1 <- simulate_exposures(..., n = 1)
#' length(x1$profiles)
#' # 1
#'
#' x100 <- simulate_exposures(..., n = 100)
#' length(x100$profiles)
#' # 100
#' ```
#'
#' There is no statistical summarization inside this function. If `n = 100`,
#' all 100 simulated profiles are retained.
#'
#' This is deliberately different from `predict_outcomes()` and
#' `compare_predictions()`, where `"samples"` refers to uncertainty draws of
#' fitted model parameters and `"summary"` summarizes those parameter-draw
#' predictions.
#'
#' ## Direct use in `compare_predictions()`
#'
#' The whole returned object can be used directly as one exposure component:
#'
#' ```
#' rain_A <- simulate_exposures(..., n = 100)
#'
#' scenario_A <- list(
#'   tmean = 25,
#'   rain = rain_A,
#'   wetness = 10
#' )
#' ```
#'
#' `compare_predictions()` recognizes `$profiles` as the canonical collection.
#' If other exposure variables in the same scenario also contain 100 profiles,
#' they are paired by simulation index:
#'
#' ```
#' tmean profile 1 + rain profile 1 + wetness profile 1
#' tmean profile 2 + rain profile 2 + wetness profile 2
#' ...
#' ```
#'
#' A variable with exactly one profile may be recycled across the common
#' scenario profile count. No Cartesian product is constructed.
#'
#' ## Chronological versus retrospective indexing
#'
#' All user-facing pattern arguments use chronological time:
#'
#' ```
#' time = 0        -> earliest/oldest observation
#' time = max_lag  -> most recent observation
#' ```
#'
#' The corresponding retrospective lag is:
#'
#' \deqn{lag = max\_lag - time.}
#'
#' No profile reversal is performed after simulation.
#'
#' ## Background generators
#'
#' For `background$dist = "normal"`, positions are independent Gaussian draws.
#'
#' For `"empirical"`, positions are sampled with replacement from the supplied
#' finite empirical values.
#'
#' For `"ar1"`, the process is initialized from its stationary marginal
#' distribution and uses innovation SD
#'
#' \deqn{
#'   sd_{innovation} = sd\sqrt{1-\phi^2}.
#' }
#'
#' Thus `background$sd` is the stationary marginal SD rather than the innovation
#' SD.
#'
#' For `"fixed"`, every background position equals `background$value`.
#'
#' ## Seeds across exposure variables
#'
#' `seed` controls one call to `simulate_exposures()`. Reusing exactly the same
#' seed in separate stochastic calls for different exposure variables can reuse
#' the same underlying pseudo-random sequence and unintentionally induce
#' dependence among simulated exposures. Use different seeds for independently
#' simulated variables unless synchronized random structure is intentional.
#'
#' ## Pattern-specific arguments
#'
#' Pattern arguments that do not belong to the selected mode/pattern are
#' rejected rather than silently ignored.
#'
#' @export
simulate_exposures <- function(
    max_lag,
    n = 1,
    mode = c("profile", "pattern"),
    pattern = c("random_times", "alternating_times", "block_times"),
    fixed_value = NULL,
    n_times = NULL,
    time_range = NULL,
    alternating_values = NULL,
    block_times = NULL,
    background = list(
      dist = "normal",
      mean = 0,
      sd = 1
    ),
    bounds = NULL,
    seed = NULL,
    cumulative = FALSE
) {

  # HELPERS

  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  is_whole_scalar <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x) &&
      x == as.integer(x)
  }

  validate_finite_scalar <- function(x, argument) {
    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x)) {
      stop(
        "`", argument, "` must be one finite numeric value.",
        call. = FALSE
      )
    }
    as.numeric(x)
  }

  validate_time_vector <- function(
    x,
    argument,
    max_lag_i,
    exact_length = NULL
  ) {
    if (!is.numeric(x) ||
        !length(x) ||
        anyNA(x) ||
        any(!is.finite(x)) ||
        any(x != as.integer(x))) {
      stop(
        "`", argument,
        "` must contain finite integer chronological times.",
        call. = FALSE
      )
    }

    x <- as.integer(x)

    if (!is.null(exact_length) &&
        length(x) != exact_length) {
      stop(
        "`", argument, "` must contain exactly ",
        exact_length, " chronological time value(s).",
        call. = FALSE
      )
    }

    if (any(x < 0L | x > max_lag_i)) {
      stop(
        "`", argument,
        "` must contain chronological times between 0 and `max_lag`.",
        call. = FALSE
      )
    }

    if (anyDuplicated(x)) {
      stop(
        "`", argument,
        "` must contain unique chronological times.",
        call. = FALSE
      )
    }

    x
  }

  reject_non_null <- function(argument_name, value, context) {
    if (!is.null(value)) {
      stop(
        "`", argument_name, "` is not used for ", context,
        ". Remove it rather than relying on an ignored argument.",
        call. = FALSE
      )
    }
  }

  # BASIC VALIDATION

  if (!is_whole_scalar(max_lag) || max_lag < 0) {
    stop(
      "`max_lag` must be one non-negative integer.",
      call. = FALSE
    )
  }
  max_lag_i <- as.integer(max_lag)

  if (!is_whole_scalar(n) || n <= 0) {
    stop(
      "`n` must be one positive integer.",
      call. = FALSE
    )
  }
  n <- as.integer(n)

  if (!is.logical(cumulative) ||
      length(cumulative) != 1L ||
      is.na(cumulative)) {
    stop(
      "`cumulative` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  mode_i <- match.arg(mode)

  # RNG CONTRACT

  seed_i <- NULL

  if (!is.null(seed)) {
    if (!is_whole_scalar(seed)) {
      stop(
        "`seed` must be NULL or one finite integer.",
        call. = FALSE
      )
    }

    seed_i <- as.integer(seed)

    had_random_seed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )

    if (had_random_seed) {
      old_random_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }

    on.exit(
      {
        if (had_random_seed) {
          assign(
            ".Random.seed",
            old_random_seed,
            envir = .GlobalEnv
          )
        } else if (exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )) {
          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )

    set.seed(seed_i)
  }

  # CHRONOLOGICAL INDEX

  profile_length <- max_lag_i + 1L

  chronological_time <- seq.int(
    0L,
    max_lag_i
  )

  retrospective_lag <- max_lag_i -
    chronological_time

  chronological_position <- seq_len(
    profile_length
  )

  profile_names <- paste0(
    "time_",
    chronological_time
  )

  # PATTERN CONTRACT

  pattern_i <- NA_character_

  if (identical(mode_i, "profile")) {
    if (!missing(pattern)) {
      stop(
        "`pattern` must not be supplied when `mode = 'profile'`.",
        call. = FALSE
      )
    }

    reject_non_null(
      "fixed_value",
      fixed_value,
      "`mode = 'profile'`"
    )
    reject_non_null(
      "n_times",
      n_times,
      "`mode = 'profile'`"
    )
    reject_non_null(
      "time_range",
      time_range,
      "`mode = 'profile'`"
    )
    reject_non_null(
      "alternating_values",
      alternating_values,
      "`mode = 'profile'`"
    )
    reject_non_null(
      "block_times",
      block_times,
      "`mode = 'profile'`"
    )

    time_range_i <- chronological_time

  } else {
    pattern_i <- match.arg(pattern)

    if (identical(pattern_i, "random_times")) {
      fixed_value <- validate_finite_scalar(
        fixed_value,
        "fixed_value"
      )

      if (!is_whole_scalar(n_times) ||
          n_times <= 0) {
        stop(
          "`n_times` must be one positive integer for ",
          "`pattern = 'random_times'`.",
          call. = FALSE
        )
      }
      n_times <- as.integer(n_times)

      reject_non_null(
        "alternating_values",
        alternating_values,
        "`pattern = 'random_times'`"
      )
      reject_non_null(
        "block_times",
        block_times,
        "`pattern = 'random_times'`"
      )

      if (is.null(time_range)) {
        time_range_i <- chronological_time
      } else {
        time_range_i <- sort(
          validate_time_vector(
            time_range,
            "time_range",
            max_lag_i
          )
        )
      }

      if (n_times > length(time_range_i)) {
        stop(
          "`n_times` cannot exceed the number of eligible chronological ",
          "positions in `time_range`.",
          call. = FALSE
        )
      }

    } else if (identical(
      pattern_i,
      "alternating_times"
    )) {
      if (!is.numeric(alternating_values) ||
          length(alternating_values) < 2L ||
          anyNA(alternating_values) ||
          any(!is.finite(alternating_values))) {
        stop(
          "`alternating_values` must contain at least two finite numeric ",
          "values for `pattern = 'alternating_times'`.",
          call. = FALSE
        )
      }

      alternating_values <- as.numeric(
        alternating_values
      )

      reject_non_null(
        "fixed_value",
        fixed_value,
        "`pattern = 'alternating_times'`"
      )
      reject_non_null(
        "n_times",
        n_times,
        "`pattern = 'alternating_times'`"
      )
      reject_non_null(
        "block_times",
        block_times,
        "`pattern = 'alternating_times'`"
      )

      if (is.null(time_range)) {
        time_range_i <- chronological_time
      } else {
        time_range_i <- sort(
          validate_time_vector(
            time_range,
            "time_range",
            max_lag_i
          )
        )
      }

      if (length(time_range_i) < 2L) {
        stop(
          "`pattern = 'alternating_times'` requires at least two selected ",
          "chronological positions.",
          call. = FALSE
        )
      }

    } else if (identical(
      pattern_i,
      "block_times"
    )) {
      fixed_value <- validate_finite_scalar(
        fixed_value,
        "fixed_value"
      )

      block_times <- sort(
        validate_time_vector(
          block_times,
          "block_times",
          max_lag_i,
          exact_length = 2L
        )
      )

      reject_non_null(
        "n_times",
        n_times,
        "`pattern = 'block_times'`"
      )
      reject_non_null(
        "alternating_values",
        alternating_values,
        "`pattern = 'block_times'`"
      )
      reject_non_null(
        "time_range",
        time_range,
        "`pattern = 'block_times'`; use `block_times` to define the interval"
      )

      time_range_i <- seq.int(
        block_times[1L],
        block_times[2L]
      )
    }
  }

  # BACKGROUND VALIDATION

  if (!is.list(background) ||
      !length(background)) {
    stop(
      "`background` must be a non-empty named list.",
      call. = FALSE
    )
  }

  background_names <- names(background)

  if (is.null(background_names) ||
      anyNA(background_names) ||
      any(!nzchar(background_names)) ||
      anyDuplicated(background_names)) {
    stop(
      "`background` must have unique, non-empty field names.",
      call. = FALSE
    )
  }

  dist <- background$dist %||% "normal"

  if (!is.character(dist) ||
      length(dist) != 1L ||
      is.na(dist) ||
      !nzchar(dist)) {
    stop(
      "`background$dist` must be one non-empty character value.",
      call. = FALSE
    )
  }

  dist <- match.arg(
    tolower(dist),
    choices = c(
      "normal",
      "empirical",
      "ar1",
      "fixed"
    )
  )

  allowed_background_fields <- switch(
    dist,
    normal = c("dist", "mean", "sd"),
    empirical = c("dist", "values"),
    ar1 = c("dist", "mean", "sd", "phi"),
    fixed = c("dist", "value")
  )

  unknown_background_fields <- setdiff(
    background_names,
    allowed_background_fields
  )

  if (length(unknown_background_fields)) {
    stop(
      "Unknown field(s) for `background$dist = '",
      dist, "'`: ",
      paste(
        unknown_background_fields,
        collapse = ", "
      ),
      ". Allowed fields are: ",
      paste(
        allowed_background_fields,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  background_parameters <- list(
    dist = dist
  )

  if (identical(dist, "normal")) {
    mean_value <- validate_finite_scalar(
      background$mean %||% 0,
      "background$mean"
    )

    sd_value <- validate_finite_scalar(
      background$sd %||% 1,
      "background$sd"
    )

    if (sd_value < 0) {
      stop(
        "`background$sd` must be non-negative.",
        call. = FALSE
      )
    }

    background_parameters$mean <- mean_value
    background_parameters$sd <- sd_value

  } else if (identical(dist, "empirical")) {
    values <- background$values

    if (!is.numeric(values) ||
        !length(values) ||
        anyNA(values) ||
        any(!is.finite(values))) {
      stop(
        "For `background$dist = 'empirical'`, `background$values` must be a ",
        "non-empty vector containing only finite numeric values.",
        call. = FALSE
      )
    }

    background_parameters$values <- as.numeric(
      values
    )

  } else if (identical(dist, "ar1")) {
    mean_value <- validate_finite_scalar(
      background$mean %||% 0,
      "background$mean"
    )

    sd_value <- validate_finite_scalar(
      background$sd %||% 1,
      "background$sd"
    )

    phi <- validate_finite_scalar(
      background$phi %||% 0.7,
      "background$phi"
    )

    if (sd_value < 0) {
      stop(
        "`background$sd` must be non-negative.",
        call. = FALSE
      )
    }

    if (abs(phi) >= 1) {
      stop(
        "For `background$dist = 'ar1'`, `background$phi` must lie strictly ",
        "inside (-1, 1).",
        call. = FALSE
      )
    }

    background_parameters$mean <- mean_value
    background_parameters$sd <- sd_value
    background_parameters$phi <- phi

  } else if (identical(dist, "fixed")) {
    fixed_background <- validate_finite_scalar(
      background$value,
      "background$value"
    )

    background_parameters$value <- fixed_background
  }

  # BOUNDS

  if (!is.null(bounds)) {
    if (!is.numeric(bounds) ||
        length(bounds) != 2L ||
        anyNA(bounds) ||
        any(!is.finite(bounds))) {
      stop(
        "`bounds` must be NULL or a finite numeric vector of length two.",
        call. = FALSE
      )
    }

    bounds <- as.numeric(bounds)

    if (bounds[1L] > bounds[2L]) {
      stop(
        "`bounds[1]` must be less than or equal to `bounds[2]`.",
        call. = FALSE
      )
    }

    if (identical(mode_i, "pattern") &&
        pattern_i %in% c(
          "random_times",
          "block_times"
        ) &&
        (
          fixed_value < bounds[1L] ||
          fixed_value > bounds[2L]
        )) {
      warning(
        "`fixed_value` lies outside `bounds` and will be clipped.",
        call. = FALSE
      )
    }

    if (identical(mode_i, "pattern") &&
        identical(pattern_i, "alternating_times") &&
        any(
          alternating_values < bounds[1L] |
          alternating_values > bounds[2L]
        )) {
      warning(
        "Some `alternating_values` lie outside `bounds` and will be clipped.",
        call. = FALSE
      )
    }
  }

  clip_bounds <- function(x) {
    if (is.null(bounds)) return(x)

    pmin(
      pmax(x, bounds[1L]),
      bounds[2L]
    )
  }

  # BACKGROUND GENERATOR

  generate_background <- function() {
    if (identical(dist, "normal")) {
      return(
        stats::rnorm(
          profile_length,
          mean = background_parameters$mean,
          sd = background_parameters$sd
        )
      )
    }

    if (identical(dist, "empirical")) {
      return(
        sample(
          background_parameters$values,
          size = profile_length,
          replace = TRUE
        )
      )
    }

    if (identical(dist, "ar1")) {
      mean_value <- background_parameters$mean
      sd_value <- background_parameters$sd
      phi <- background_parameters$phi

      x <- numeric(profile_length)

      x[1L] <- stats::rnorm(
        1L,
        mean = mean_value,
        sd = sd_value
      )

      innovation_sd <- sd_value *
        sqrt(1 - phi^2)

      if (profile_length > 1L) {
        for (i in seq.int(
          2L,
          profile_length
        )) {
          x[i] <-
            mean_value +
            phi * (
              x[i - 1L] -
                mean_value
            ) +
            stats::rnorm(
              1L,
              mean = 0,
              sd = innovation_sd
            )
        }
      }

      return(x)
    }

    rep(
      background_parameters$value,
      profile_length
    )
  }

  # ONE SIMULATION

  simulate_one <- function(simulation_id) {
    x <- generate_background()

    if (length(x) != profile_length ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "Background generation produced an invalid exposure profile.",
        call. = FALSE
      )
    }

    x <- clip_bounds(x)

    selected_times <- integer(0)

    if (identical(mode_i, "pattern")) {
      if (identical(pattern_i, "random_times")) {
        selected_times <- sort(
          sample(
            time_range_i,
            size = n_times,
            replace = FALSE
          )
        )

        x[selected_times + 1L] <- fixed_value

      } else if (identical(
        pattern_i,
        "alternating_times"
      )) {
        selected_times <- sort(
          time_range_i
        )

        x[selected_times + 1L] <- rep(
          alternating_values,
          length.out = length(
            selected_times
          )
        )

      } else if (identical(
        pattern_i,
        "block_times"
      )) {
        selected_times <- seq.int(
          block_times[1L],
          block_times[2L]
        )

        x[selected_times + 1L] <- fixed_value
      }
    }

    x <- clip_bounds(x)

    if (cumulative) {
      x <- cumsum(x)
    }

    if (length(x) != profile_length ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "Exposure simulation produced a non-finite final profile.",
        call. = FALSE
      )
    }

    names(x) <- profile_names

    list(
      profile = x,
      simulation_meta = list(
        simulation = simulation_id,
        selected_times = selected_times,
        selected_positions = selected_times + 1L,
        selected_internal_lags =
          max_lag_i - selected_times
      )
    )
  }

  # RUN SIMULATIONS

  simulations <- lapply(
    seq_len(n),
    simulate_one
  )

  profiles <- lapply(
    simulations,
    `[[`,
    "profile"
  )

  names(profiles) <- paste0(
    "simulation_",
    seq_len(n)
  )

  simulation_meta <- lapply(
    simulations,
    `[[`,
    "simulation_meta"
  )

  names(simulation_meta) <- names(
    profiles
  )

  simulation_matrix <- do.call(
    rbind,
    profiles
  )

  if (!is.matrix(simulation_matrix) ||
      nrow(simulation_matrix) != n ||
      ncol(simulation_matrix) != profile_length ||
      anyNA(simulation_matrix) ||
      any(!is.finite(simulation_matrix))) {
    stop(
      "Internal simulation assembly failed.",
      call. = FALSE
    )
  }

  colnames(simulation_matrix) <- profile_names

  # LONG TABULAR REPRESENTATION

  simulation_data <- data.frame(
    simulation = rep(
      seq_len(n),
      each = profile_length
    ),
    position = rep(
      chronological_position,
      times = n
    ),
    time = rep(
      chronological_time,
      times = n
    ),
    lag = rep(
      retrospective_lag,
      times = n
    ),
    profile = as.vector(
      t(simulation_matrix)
    ),
    stringsAsFactors = FALSE
  )

  # COMMON METADATA

  meta <- list(
    max_lag = max_lag_i,
    profile_length = profile_length,
    profile_order = "chronological",
    chronological_time = chronological_time,
    retrospective_lag = retrospective_lag,
    mode = mode_i,
    pattern = pattern_i,
    background = background_parameters,
    bounds = bounds,
    seed = seed_i,
    time_range = time_range_i,
    cumulative = cumulative,
    bounds_apply_to = if (cumulative) {
      "increments_before_cumulative_sum"
    } else {
      "profile_values"
    },
    n = n,
    n_profiles = n,
    profile_contract =
      "canonical_profiles_list_chronological",
    lag_contract =
      "lag_0_most_recent",
    rng_contract =
      "local_seed_restores_caller_rng",
    simulation_contract =
      "n_exposure_profiles_no_model_parameter_uncertainty"
  )

  if (n == 1L) {
    meta$selected_times <-
      simulation_meta[[1L]]$selected_times
    meta$selected_positions <-
      simulation_meta[[1L]]$selected_positions
    meta$selected_internal_lags <-
      simulation_meta[[1L]]$selected_internal_lags
  }

  # FINAL OBJECT

  out <- list(
    profiles = profiles,
    profile = if (n == 1L) {
      profiles[[1L]]
    } else {
      NULL
    },
    simulation_data = simulation_data,
    meta = meta,
    simulation_meta = simulation_meta
  )

  class(out) <- c(
    "epiexposure_simulated_exposures",
    "list"
  )

  attr(
    out,
    "epiexposure_profile_order"
  ) <- "chronological"

  attr(
    out,
    "epiexposure_lag_scale"
  ) <- "lag_0_most_recent"

  attr(
    out,
    "epiexposure_max_lag"
  ) <- max_lag_i

  attr(
    out,
    "epiexposure_n_profiles"
  ) <- n

  attr(
    out,
    "epiexposure_profile_contract"
  ) <- "canonical_profiles_list_chronological"

  attr(
    out,
    "epiexposure_simulation_contract"
  ) <- "n_exposure_profiles_no_model_parameter_uncertainty"

  out
}
