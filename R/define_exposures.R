#' Define DLNM exposure templates
#'
#' Defines the training cross-basis parameterization for one or more exposure
#' variables from long-format epidemic histories.
#'
#' `define_exposures()` is the first step in the standard EpiExposure workflow:
#'
#' ```
#' templates <- define_exposures(...)
#' design    <- build_design(data, templates)
#' fit       <- fit_epidlnm(...)
#' ```
#'
#' The function estimates **basis definitions**, not disease effects. For each
#' exposure it uses all training epidemics jointly to determine any
#' data-dependent exposure-basis features (for example spline knots and boundary
#' knots), while `dlnm::crossbasis(..., group = ...)` keeps the individual
#' epidemic histories as independent time series so lags never cross from one
#' epidemic into another.
#'
#' @param data Non-empty long-format data frame containing `epi_id`, `time`, and
#'   every exposure named in `vars`.
#'
#'   Rows do not need to be pre-sorted. They are ordered internally by
#'   `epi_id` and `time`.
#'
#'   Each epidemic must form one complete exposure history containing exactly
#'   `max_lag + 1` time points. When more than one time point is present, spacing
#'   must be constant within epidemics and identical across epidemics. Missing
#'   exposure values are not allowed in EpiExposure v1. Because all exposures
#'   are finite columns of the same validated rows, every fitted exposure uses
#'   the same temporal support and history length.
#' @param vars Unique character vector naming exposure variables in `data`.
#' @param max_lag Non-negative integer maximum lag. The legacy form
#'   `c(0, L)` is also accepted, but EpiExposure v1 requires the lag range to
#'   start at 0.
#' @param df_var Positive integer controlling the exposure-response basis when
#'   `fun_var` is `"ns"` or `"bs"`. For `fun_var = "poly"`, `df_var` is used
#'   as the polynomial degree. It is ignored for `fun_var = "lin"`.
#' @param df_lag Positive integer controlling the lag-response basis when
#'   `fun_lag` is `"ns"` or `"bs"`. For `fun_lag = "poly"`, `df_lag` is used
#'   as the polynomial degree. It is ignored for `fun_lag = "lin"`.
#' @param fun_var Character exposure-basis function. Supported unpenalized
#'   EpiExposure v1 choices are `"ns"`, `"bs"`, `"poly"`, and `"lin"`.
#' @param fun_lag Character lag-basis function. Supported unpenalized
#'   EpiExposure v1 choices are `"ns"`, `"bs"`, `"poly"`, and `"lin"`.
#'
#'   Penalized dlnm bases such as `"ps"` are not supported by the current
#'   engine-agnostic EpiExposure fitting contract. A P-spline transformation is
#'   not, by itself, a penalized DLNM: the corresponding penalty must also be
#'   supplied during model fitting (for example through `dlnm::cbPen()` and
#'   `mgcv::gam(..., paraPen = ...)`, or through a dedicated cross-basis smooth).
#'   EpiExposure v1 therefore stops explicitly rather than fitting an
#'   unintentionally unpenalized `"ps"` basis.
#'
#' @return A named list of `dlnm::crossbasis` training templates, one per
#'   exposure. The list carries:
#'
#'   - `attr(x, "spec")`: effective per-exposure basis specification;
#'   - `attr(x, "epiexposure_vars")`: exposure names;
#'   - `attr(x, "epiexposure_max_lag")`: the common fitted maximum lag;
#'   - `attr(x, "epiexposure_history_length")`: the exact required exposure
#'     history length, `max_lag + 1`;
#'   - `attr(x, "epiexposure_history_contract")`:
#'     `"all_fitted_exposures_same_exact_max_lag_plus_one"`;
#'   - `attr(x, "epiexposure_time_step")`: the common temporal spacing, or `NA`
#'     when `max_lag = 0` because one observation has no estimable spacing;
#'   - `attr(x, "epiexposure_template_contract")`:
#'     `"training_grouped_crossbasis"`;
#'   - `attr(x, "epiexposure_profile_order")`: `"chronological"`.
#'
#'   Each individual cross-basis also stores:
#'
#'   - `attr(cb, "cb_colnames")`: native dlnm cross-basis column names;
#'   - `attr(cb, "epiexposure_canonical_cb_colnames")`: the canonical future
#'     design names `cb_<variable>_<index>`;
#'   - `attr(cb, "epiexposure_variable")`: exposure name.
#'
#' @details
#' ## Independent epidemic series
#'
#' Earlier EpiExposure code concatenated epidemic histories and inserted
#' `max_lag` missing values between them before calling `crossbasis()`. The
#' current implementation uses the native `group` argument of
#' `dlnm::crossbasis()` instead.
#'
#' This is the intended dlnm representation for several independent time
#' series. After internal sorting, every epidemic remains consecutive and
#' complete, while exposure values from all training epidemics still contribute
#' to the common exposure-basis definition. Thus, for a spline exposure basis,
#' data-dependent knots are estimated from the pooled **training exposure
#' distribution**, but lagged histories never cross epidemic boundaries.
#'
#' ## Temporal requirements
#'
#' A vector passed to `dlnm::crossbasis()` represents one complete exposure
#' history. `define_exposures()` therefore rejects:
#'
#' - duplicated time values within an epidemic;
#' - irregular spacing within an epidemic when more than one time point exists;
#' - different time steps across epidemics;
#' - epidemics with either fewer **or more** than `max_lag + 1` observations.
#'
#' If the common maximum lag is `L`, every epidemic must contain exactly
#'
#' \deqn{
#'   L + 1
#' }
#'
#' observations. Histories are never truncated to a trailing window, padded,
#' or silently realigned. Because all fitted exposures are columns of the same
#' long-format rows, all variables necessarily have the same temporal length
#' and support within each epidemic.
#'
#' For `max_lag = 0`, the exact history length is one observation. In that
#' special case no time interval exists to estimate, so
#' `epiexposure_time_step` is stored as `NA`.
#'
#' These are the same temporal assumptions enforced downstream by
#' `build_design()`.
#'
#' ## Exposure and lag intercepts
#'
#' The exposure basis is constructed without an intercept. This follows the
#' standard DLNM identifiability requirement and avoids rank deficiency when a
#' model intercept is present.
#'
#' The lag basis includes an intercept. This is the standard cross-basis
#' parameterization for the lag dimension and allows an effect distributed over
#' the full fitted lag window.
#'
#' ## Effective versus requested basis arguments
#'
#' `dlnm::crossbasis()` calls `onebasis()` internally and may normalize or
#' modify basis arguments. Accordingly, EpiExposure stores the **effective**
#' `argvar`, `arglag`, and lag attributes returned by the completed cross-basis,
#' rather than assuming that the original request is the final
#' parameterization.
#'
#' These effective attributes are the authoritative training template reused by
#' `build_design()`, `predict_outcomes()`, cross-validation, scenario functions,
#' and effect summaries. Prediction data never redefine knots or boundary knots.
#'
#' ## Scope of `df_var` and `df_lag`
#'
#' For spline bases, the user-supplied degrees of freedom determine the requested
#' flexibility, but the returned `spec` records the effective basis arguments.
#' For `"poly"`, the corresponding `df_*` argument is interpreted as polynomial
#' degree. For `"lin"`, the corresponding `df_*` value has no effect.
#'
#' @export
define_exposures <- function(
    data,
    vars,
    max_lag,
    df_var = 4,
    df_lag = 4,
    fun_var = "ns",
    fun_lag = "ns"
) {

  # ==========================================================================
  # SMALL VALIDATORS
  # ==========================================================================

  valid_scalar_character <- function(x) {
    is.character(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      nzchar(x)
  }

  is_whole_scalar <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x) &&
      x == as.integer(x)
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
    stop(
      "`data` must be a non-empty data.frame.",
      call. = FALSE
    )
  }

  missing_basic <- setdiff(
    c("epi_id", "time"),
    names(data)
  )

  if (length(missing_basic)) {
    stop(
      "`data` is missing required column(s): ",
      paste(missing_basic, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.character(vars) ||
      !length(vars) ||
      anyNA(vars) ||
      any(!nzchar(vars)) ||
      anyDuplicated(vars)) {
    stop(
      "`vars` must contain unique non-empty exposure names.",
      call. = FALSE
    )
  }

  if (any(vars %in% c("epi_id", "time", "y"))) {
    stop(
      "Exposure variables cannot use the reserved EpiExposure names ",
      "'epi_id', 'time', or 'y'.",
      call. = FALSE
    )
  }

  missing_vars <- setdiff(
    vars,
    names(data)
  )

  if (length(missing_vars)) {
    stop(
      "Exposure variable(s) not found in `data`: ",
      paste(missing_vars, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (anyNA(data$epi_id)) {
    stop(
      "`data$epi_id` cannot contain missing values.",
      call. = FALSE
    )
  }

  if (!is.numeric(data$time) ||
      anyNA(data$time) ||
      any(!is.finite(data$time))) {
    stop(
      "`data$time` must contain only finite numeric values.",
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

    if (length(unique(x)) < 2L) {
      stop(
        "Exposure variable '",
        variable,
        "' is constant. A cross-basis effect cannot be estimated from a ",
        "constant exposure.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # LAG AND BASIS ARGUMENTS
  # ==========================================================================

  if (!is.numeric(max_lag) ||
      !length(max_lag) ||
      length(max_lag) > 2L ||
      anyNA(max_lag) ||
      any(!is.finite(max_lag)) ||
      any(max_lag < 0) ||
      any(max_lag != as.integer(max_lag))) {
    stop(
      "`max_lag` must be one non-negative integer or the legacy lag range ",
      "`c(0, L)`.",
      call. = FALSE
    )
  }

  max_lag <- as.integer(max_lag)

  if (length(max_lag) == 2L &&
      min(max_lag) != 0L) {
    stop(
      "EpiExposure v1 requires the fitted lag range to start at lag 0. ",
      "Use `max_lag = L` or `max_lag = c(0, L)`.",
      call. = FALSE
    )
  }

  maximum_lag <- as.integer(
    max(max_lag)
  )

  if (!is_whole_scalar(df_var) ||
      df_var < 1L) {
    stop(
      "`df_var` must be a positive integer.",
      call. = FALSE
    )
  }
  df_var <- as.integer(df_var)

  if (!is_whole_scalar(df_lag) ||
      df_lag < 1L) {
    stop(
      "`df_lag` must be a positive integer.",
      call. = FALSE
    )
  }
  df_lag <- as.integer(df_lag)

  if (!valid_scalar_character(fun_var)) {
    stop(
      "`fun_var` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!valid_scalar_character(fun_lag)) {
    stop(
      "`fun_lag` must be one non-empty character value.",
      call. = FALSE
    )
  }

  fun_var <- tolower(fun_var)
  fun_lag <- tolower(fun_lag)

  penalized_functions <- c("ps", "cr")

  if (fun_var %in% penalized_functions ||
      fun_lag %in% penalized_functions) {
    requested_penalized <- unique(
      c(
        if (fun_var %in% penalized_functions) fun_var,
        if (fun_lag %in% penalized_functions) fun_lag
      )
    )

    stop(
      "Penalized dlnm basis function(s) ",
      paste(
        paste0("'", requested_penalized, "'"),
        collapse = ", "
      ),
      " are not supported by the engine-agnostic EpiExposure v1 fitting ",
      "contract. `ps`/`cr` require their penalty matrices to be propagated ",
      "during fitting; EpiExposure does not silently fit them as unpenalized ",
      "basis columns.",
      call. = FALSE
    )
  }

  allowed_unpenalized <- c(
    "ns",
    "bs",
    "poly",
    "lin"
  )

  if (!fun_var %in% allowed_unpenalized) {
    stop(
      "Unsupported `fun_var = ",
      sQuote(fun_var),
      "`. Supported EpiExposure v1 exposure bases are: ",
      paste(allowed_unpenalized, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!fun_lag %in% allowed_unpenalized) {
    stop(
      "Unsupported `fun_lag = ",
      sQuote(fun_lag),
      "`. Supported EpiExposure v1 lag bases are: ",
      paste(allowed_unpenalized, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (maximum_lag == 0L &&
      !identical(fun_lag, "lin")) {
    stop(
      "When `max_lag = 0`, use `fun_lag = 'lin'`. A non-linear lag basis ",
      "cannot be identified from a single lag value.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dlnm", quietly = TRUE)) {
    stop(
      "Package 'dlnm' is required by `define_exposures()`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # TEMPORAL REGULARITY
  # ==========================================================================

  data_ordered <- data[
    order(data$epi_id, data$time),
    ,
    drop = FALSE
  ]
  rownames(data_ordered) <- NULL

  group_key <- as.character(
    data_ordered$epi_id
  )

  group_rows <- split(
    seq_len(nrow(data_ordered)),
    group_key,
    drop = TRUE
  )

  if (!length(group_rows)) {
    stop(
      "No epidemic groups were found in `data`.",
      call. = FALSE
    )
  }

  first_rows <- vapply(
    group_rows,
    function(idx) idx[1L],
    integer(1)
  )
  group_rows <- group_rows[
    order(first_rows)
  ]

  common_time_step <- NA_real_
  tolerance <- sqrt(.Machine$double.eps)
  required_history_length <- maximum_lag + 1L
  history_contract <- "all_fitted_exposures_same_exact_max_lag_plus_one"

  for (i in seq_along(group_rows)) {
    idx <- group_rows[[i]]
    current_id <- data_ordered$epi_id[
      idx[1L]
    ]
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
        " observations for max_lag = ",
        maximum_lag,
        "; received ",
        length(tt),
        ". Histories are not truncated, padded, or silently realigned.",
        call. = FALSE
      )
    }

    if (required_history_length > 1L) {
      dt <- diff(tt)

      if (any(dt <= 0) ||
          any(!is.finite(dt))) {
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
        tolerance *
        pmax(
          1,
          abs(dt),
          abs(current_step)
        )
      )) {
        stop(
          "Irregular time spacing was detected within epi_id = ",
          as.character(current_id),
          ". DLNM vector histories must be complete and equally spaced.",
          call. = FALSE
        )
      }

      if (is.na(common_time_step)) {
        common_time_step <- as.numeric(
          current_step
        )
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
          ", while the common time step is ",
          format(common_time_step),
          ".",
          call. = FALSE
        )
      }
    }
  }

  # Build a consecutive factor exactly matching the sorted rows. This is the
  # native dlnm mechanism for several independent time series.
  group_factor <- factor(
    as.character(data_ordered$epi_id),
    levels = unique(
      as.character(data_ordered$epi_id)
    )
  )

  # ==========================================================================
  # BASIS-ARGUMENT CONSTRUCTORS
  # ==========================================================================

  make_argvar <- function(fun, df) {
    if (identical(fun, "ns")) {
      return(
        list(
          fun = "ns",
          df = df,
          intercept = FALSE
        )
      )
    }

    if (identical(fun, "bs")) {
      return(
        list(
          fun = "bs",
          df = df,
          intercept = FALSE
        )
      )
    }

    if (identical(fun, "poly")) {
      return(
        list(
          fun = "poly",
          degree = df,
          intercept = FALSE
        )
      )
    }

    if (identical(fun, "lin")) {
      return(
        list(
          fun = "lin",
          intercept = FALSE
        )
      )
    }

    stop(
      "Internal error while constructing the exposure basis.",
      call. = FALSE
    )
  }

  make_arglag <- function(fun, df) {
    if (identical(fun, "ns")) {
      return(
        list(
          fun = "ns",
          df = df,
          intercept = TRUE
        )
      )
    }

    if (identical(fun, "bs")) {
      return(
        list(
          fun = "bs",
          df = df,
          intercept = TRUE
        )
      )
    }

    if (identical(fun, "poly")) {
      return(
        list(
          fun = "poly",
          degree = df,
          intercept = TRUE
        )
      )
    }

    if (identical(fun, "lin")) {
      return(
        list(
          fun = "lin",
          intercept = TRUE
        )
      )
    }

    stop(
      "Internal error while constructing the lag basis.",
      call. = FALSE
    )
  }

  requested_argvar <- make_argvar(
    fun_var,
    df_var
  )
  requested_arglag <- make_arglag(
    fun_lag,
    df_lag
  )

  # ==========================================================================
  # CREATE TRAINING CROSS-BASIS TEMPLATES
  # ==========================================================================

  cb_templates <- stats::setNames(
    vector("list", length(vars)),
    vars
  )

  epiexposure_spec <- stats::setNames(
    vector("list", length(vars)),
    vars
  )

  for (variable in vars) {
    x <- data_ordered[[variable]]

    if (length(x) != nrow(data_ordered)) {
      stop(
        "Internal exposure-history alignment failed for variable '",
        variable,
        "'. All exposures must occupy the same validated time rows.",
        call. = FALSE
      )
    }

    for (i in seq_along(group_rows)) {
      idx <- group_rows[[i]]

      if (length(x[idx]) != required_history_length) {
        stop(
          "Exposure history for '",
          variable,
          "' and epi_id = ",
          as.character(data_ordered$epi_id[idx[1L]]),
          " must contain exactly ",
          required_history_length,
          " observations for max_lag = ",
          maximum_lag,
          "; received ",
          length(x[idx]),
          ".",
          call. = FALSE
        )
      }
    }

    cb <- tryCatch(
      dlnm::crossbasis(
        x = x,
        lag = c(0L, maximum_lag),
        argvar = requested_argvar,
        arglag = requested_arglag,
        group = group_factor
      ),
      error = function(e) {
        stop(
          "Could not construct the training cross-basis for exposure '",
          variable,
          "': ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (!inherits(cb, "crossbasis")) {
      stop(
        "`dlnm::crossbasis()` did not return a valid `crossbasis` object for ",
        "exposure '",
        variable,
        "'.",
        call. = FALSE
      )
    }

    effective_argvar <- attr(
      cb,
      "argvar",
      exact = TRUE
    )
    effective_arglag <- attr(
      cb,
      "arglag",
      exact = TRUE
    )
    effective_lag <- attr(
      cb,
      "lag",
      exact = TRUE
    )
    effective_range <- attr(
      cb,
      "range",
      exact = TRUE
    )
    effective_df <- attr(
      cb,
      "df",
      exact = TRUE
    )

    if (!is.list(effective_argvar) ||
        !length(effective_argvar)) {
      stop(
        "Constructed cross-basis for exposure '",
        variable,
        "' is missing effective `argvar` metadata.",
        call. = FALSE
      )
    }

    if (!is.list(effective_arglag) ||
        !length(effective_arglag)) {
      stop(
        "Constructed cross-basis for exposure '",
        variable,
        "' is missing effective `arglag` metadata.",
        call. = FALSE
      )
    }

    if (!is.numeric(effective_lag) ||
        length(effective_lag) != 2L ||
        anyNA(effective_lag) ||
        any(!is.finite(effective_lag)) ||
        as.integer(min(effective_lag)) != 0L ||
        as.integer(max(effective_lag)) != maximum_lag) {
      stop(
        "Constructed cross-basis for exposure '",
        variable,
        "' returned an unexpected lag definition.",
        call. = FALSE
      )
    }

    if (nrow(cb) != nrow(data_ordered) ||
        ncol(cb) < 1L) {
      stop(
        "Constructed cross-basis for exposure '",
        variable,
        "' has an unexpected dimension.",
        call. = FALSE
      )
    }

    native_colnames <- colnames(cb)

    if (is.null(native_colnames) ||
        length(native_colnames) != ncol(cb) ||
        anyNA(native_colnames) ||
        any(!nzchar(native_colnames)) ||
        anyDuplicated(native_colnames)) {
      stop(
        "Constructed cross-basis for exposure '",
        variable,
        "' has invalid native column names.",
        call. = FALSE
      )
    }

    canonical_colnames <- paste0(
      "cb_",
      variable,
      "_",
      seq_len(ncol(cb))
    )

    attr(cb, "cb_colnames") <-
      native_colnames

    attr(
      cb,
      "epiexposure_canonical_cb_colnames"
    ) <- canonical_colnames

    attr(
      cb,
      "epiexposure_variable"
    ) <- variable

    attr(
      cb,
      "epiexposure_time_step"
    ) <- common_time_step

    attr(
      cb,
      "epiexposure_max_lag"
    ) <- maximum_lag

    attr(
      cb,
      "epiexposure_history_length"
    ) <- required_history_length

    attr(
      cb,
      "epiexposure_history_contract"
    ) <- history_contract

    cb_templates[[variable]] <- cb

    epiexposure_spec[[variable]] <- list(
      max_lag = maximum_lag,
      history_length = required_history_length,
      history_contract = history_contract,
      lag = as.integer(effective_lag),
      argvar = effective_argvar,
      arglag = effective_arglag,
      range = effective_range,
      df = effective_df,
      native_colnames = native_colnames,
      canonical_colnames = canonical_colnames,
      requested = list(
        fun_var = fun_var,
        df_var = if (
          identical(fun_var, "lin")
        ) {
          NULL
        } else {
          df_var
        },
        fun_lag = fun_lag,
        df_lag = if (
          identical(fun_lag, "lin")
        ) {
          NULL
        } else {
          df_lag
        }
      )
    )
  }

  # ==========================================================================
  # LIST-LEVEL CONTRACT
  # ==========================================================================

  attr(
    cb_templates,
    "spec"
  ) <- epiexposure_spec

  attr(
    cb_templates,
    "epiexposure_vars"
  ) <- vars

  attr(
    cb_templates,
    "epiexposure_max_lag"
  ) <- maximum_lag

  attr(
    cb_templates,
    "epiexposure_history_length"
  ) <- required_history_length

  attr(
    cb_templates,
    "epiexposure_history_contract"
  ) <- history_contract

  attr(
    cb_templates,
    "epiexposure_time_step"
  ) <- common_time_step

  attr(
    cb_templates,
    "epiexposure_template_contract"
  ) <- "training_grouped_crossbasis"

  attr(
    cb_templates,
    "epiexposure_profile_order"
  ) <- "chronological"

  attr(
    cb_templates,
    "epiexposure_group_levels"
  ) <- levels(group_factor)

  cb_templates
}
