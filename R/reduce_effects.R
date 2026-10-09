#' Reduce a fitted DLNM to a one-dimensional association
#'
#' Reduces one fitted EpiExposure distributed lag non-linear model (DLNM)
#' exposure-lag association using `dlnm::crossreduce()`.
#'
#' The function is a **reduction/reparameterization** of the fitted DLNM, not a
#' new model fit. The fitted two-dimensional exposure-by-lag association is
#' re-expressed using the one-dimensional basis for either the exposure
#' dimension or the lag dimension.
#'
#' Critically, `reduce_effects()` uses the exact `crossbasis` object and
#' coefficient mapping stored by the current `fit_epidlnm()` implementation.
#' It never re-estimates knots, boundary knots, basis dimensions, lag
#' parameterization, or coefficient alignment from the supplied `data`.
#'
#' @param fit Fitted model returned by the current `fit_epidlnm()`.
#'
#' @param data Non-empty long-format exposure data containing `group`, `time`,
#'   and **every exposure variable fitted in `fit`**.
#'
#'   `data` is used only to:
#'
#'   - validate the exact EpiExposure temporal-profile contract;
#'   - calculate a method-based reference value when requested;
#'   - validate exposure support for `value`, `at`, and the reference.
#'
#'   It is **not** used to reconstruct or re-estimate the fitted cross-basis.
#'
#'   Every group must contain exactly
#'
#'   \deqn{
#'     profile\_length = max\_lag + 1
#'   }
#'
#'   rows, with finite complete exposure profiles, unique chronological times,
#'   regular spacing within groups, and the same spacing across groups.
#'
#' @param vars Character scalar naming one exposure variable that was fitted in
#'   `fit`.
#'
#' @param group Character scalar naming the grouping column in `data`. Default
#'   is `"epi_id"`.
#'
#' @param time Character scalar naming the chronological time column in `data`.
#'   Default is `"time"`. The column must be finite numeric.
#'
#' @param type Reduction type:
#'
#'   - `"overall"`: overall cumulative exposure-response association across the
#'     complete fitted lag window;
#'   - `"lag"`: exposure-response association at one specified lag;
#'   - `"var"`: lag-response association at one specified exposure value.
#'
#'   These meanings follow `dlnm::crossreduce()`.
#'
#' @param value `NULL` for `type = "overall"`. For `type = "lag"`, one finite
#'   numeric lag coordinate within `0:max_lag`. For `type = "var"`, one finite
#'   numeric exposure value at which the lag-response association is reduced.
#'
#' @param scale Character. One of `"percent"`, `"response"`, or `"link"`.
#'
#'   This argument is retained for backward compatibility with the original
#'   EpiExposure API:
#'
#'   - `"link"` returns the centered DLNM association contrast `eta` unchanged;
#'   - `"response"` is a legacy label for the directly transformed
#'     **association measure**, not an absolute expected response. For a log
#'     link it returns `exp(eta)`, a response ratio; for a logit link it returns
#'     `exp(eta)`, an odds ratio; for an identity link it returns the additive
#'     response-scale difference `eta`;
#'   - `"percent"` returns `100 * (exp(eta) - 1)` for log and logit links.
#'
#'   For a log link, `"percent"` is the percent relative change in the expected
#'   response. For a logit link, it is the percent change in **odds**, not a
#'   percentage-point change in probability.
#'
#'   A centered DLNM contrast alone does not identify an absolute response
#'   probability/mean because the model intercept and the complete joint
#'   reference profile are not part of the reduced contrast. Therefore
#'   `"response"` is unavailable for probit, complementary-log-log, inverse, or
#'   other links for which no direct association-scale transformation is
#'   defined here. Use `summarise_effects()` when absolute response-scale
#'   `baseline`, `predicted`, and `delta` quantities are required.
#'
#' @param uncertainty Logical scalar. If `FALSE`, use the central
#'   fixed/population parameter estimate. If `TRUE`, propagate joint
#'   fixed/population parameter uncertainty draw by draw.
#'
#' @param output Character. `"summary"` or `"samples"`.
#'
#'   `"samples"` requires `uncertainty = TRUE` and returns one row per
#'   parameter draw and reduction coordinate. `"summary"` returns deterministic
#'   central-parameter values when `uncertainty = FALSE`, or median, standard
#'   deviation, and empirical quantiles when `uncertainty = TRUE`.
#'
#' @param n_samples Positive integer number of parameter draws when
#'   `uncertainty = TRUE`. Default is `1000`; at least two draws are required.
#'
#' @param seed `NULL` or one non-negative finite integer controlling parameter
#'   sampling when `uncertainty = TRUE`. The caller's previous random-number
#'   state is restored when the function exits.
#'
#' @param ref Reference exposure specification used as `cen` in
#'   `dlnm::crossreduce()`. Default is
#'   `list(method = "median", value = NULL)`.
#'
#'   Supported method-based forms are:
#'
#'   - `list(method = "median", value = NULL)`;
#'   - `list(method = "percentile", value = p)`, with `0 <= p <= 1`;
#'   - `list(method = "fixed", value = x)`, with finite numeric `x`.
#'
#'   Unlike `summarise_effects()`, this function reduces only one focal
#'   exposure at a time, so a joint multi-exposure reference profile is not
#'   required. The selected scalar reference is passed explicitly to
#'   `crossreduce(cen = ...)`; the function never relies on the implicit
#'   mid-range centering defaults of `dlnm`.
#'
#' @param at Optional finite numeric vector of exposure values at which the
#'   reduced exposure-response association is evaluated for
#'   `type = "overall"` or `type = "lag"`.
#'
#'   If `NULL`, the prediction grid is selected by `dlnm::crossreduce()` from
#'   the stored fitted basis range. `at` is not applicable to `type = "var"`
#'   because that reduction is evaluated over lag coordinates.
#'
#' @param interval_probs Numeric vector of length two defining the empirical
#'   uncertainty interval. Default is `c(0.025, 0.975)`.
#'
#' @param extrapolation Character controlling exposure values outside the
#'   range stored with the fitted cross-basis: `"error"` (default), `"warn"`,
#'   or `"allow"`.
#'
#'   The setting applies to the reference value, `value` when
#'   `type = "var"`, and user-supplied `at`. It never re-estimates the basis.
#'
#' @return A data frame.
#'
#'   The reduction coordinate is stored in `x`:
#'
#'   - for `type = "overall"` and `type = "lag"`, `x` is exposure;
#'   - for `type = "var"`, `x` is lag.
#'
#'   Every output contains:
#'
#'   - `x`: reduction coordinate;
#'   - `eta`: centered reduced association on the fitted link/linear-predictor
#'     scale;
#'   - `effect`: requested representation of that association;
#'   - `type`, `value`, `reference`, `scale`, and `var`.
#'
#'   With `uncertainty = TRUE, output = "summary"`, the result additionally
#'   contains:
#'
#'   - `eta_sd`, `low`, `high`: SD and empirical interval for `eta`;
#'   - `effect_sd`, `low_eff`, `high_eff`: SD and empirical interval for the
#'     draw-by-draw transformed effect.
#'
#'   With `output = "samples"`, `sample` identifies the joint parameter draw.
#'
#'   Attributes store the central reduced one-dimensional coefficients and
#'   basis returned by `dlnm::crossreduce()`, the reference and transformation
#'   contracts, and the inherited EpiExposure temporal metadata.
#'
#' @details
#' ## What `crossreduce()` does
#'
#' A fitted DLNM is parameterized by a two-dimensional cross-basis. The
#' `dlnm::crossreduce()` operation re-expresses that fit using modified
#' coefficients for a one-dimensional basis. It can produce:
#'
#' - an overall cumulative exposure-response summary (`type = "overall"`);
#' - an exposure-response summary at one lag (`type = "lag"`);
#' - a lag-response summary at one exposure value (`type = "var"`).
#'
#' The reduction is algebraic: it does not refit the epidemiological model.
#'
#' ## Stored fitted basis is authoritative
#'
#' The `crossbasis` used here is the object stored in
#' `attr(fit, "epiexposure_basis_objects")`. Canonical EpiExposure cross-basis
#' coefficient names are mapped to the stored native basis order through the
#' metadata created by `define_exposures()`, `build_design()`, and
#' `fit_epidlnm()`.
#'
#' `data` is never pooled to create another `crossbasis`. This prevents
#' prediction/reduction data from redefining data-dependent spline knots,
#' boundary knots, ranges, or column order.
#'
#' ## Exact common-profile contract
#'
#' Current EpiExposure fits satisfy:
#'
#' \deqn{
#'   profile\_length = max\_lag + 1
#' }
#'
#' and every fitted exposure uses the same `max_lag`.
#'
#' `reduce_effects()` validates this metadata and requires every supplied group
#' to contain exactly that many rows. Longer profiles are not truncated and
#' shorter profiles are not padded. Every fitted exposure must be present in
#' `data` on the same rows, so the function cannot silently compare or derive
#' references from profiles with different temporal support.
#'
#' For `max_lag = 0`, `history_length = 1`; no temporal step can be inferred
#' from a one-row profile and the stored time step is expected to be `NA`.
#'
#' ## Centering and reference
#'
#' `dlnm::crossreduce()` has its own default centering rules when `cen` is not
#' supplied. EpiExposure deliberately does not use those implicit defaults.
#' The `ref` argument resolves one explicit focal-exposure reference value and
#' passes it as `cen`.
#'
#' Consequently, `eta = 0` is the neutral centered association and, for log or
#' logit links, `exp(eta) = 1` is the neutral exponentiated association.
#'
#' ## Link scale versus transformed association
#'
#' `eta` is a contrast on the model's additive linear-predictor scale. Applying
#' an inverse link directly to a contrast is generally **not** the same as
#' calculating an expected response.
#'
#' In particular, for a logit model:
#'
#' \deqn{
#'   exp(eta)
#' }
#'
#' is an odds ratio, whereas
#'
#' \deqn{
#'   plogis(eta)
#' }
#'
#' is not the response probability associated with the fitted model unless the
#' omitted baseline linear predictor were exactly zero. Therefore this function
#' never applies `plogis()` directly to the reduced contrast.
#'
#' ## Uncertainty contract
#'
#' When `uncertainty = FALSE`, the reduction uses the central parameter estimate:
#'
#' - fitted fixed coefficients for frequentist engines;
#' - posterior-mean fixed/population coefficients for Bayesian engines.
#'
#' When `uncertainty = TRUE`, EpiExposure obtains one joint
#' fixed/population-parameter draw matrix from the shared engine helpers. The
#' same row of that matrix defines one coherent draw. The focal cross-basis
#' coefficients are then reduced by `crossreduce()` draw by draw.
#'
#' All transformations are performed **before** uncertainty summaries are
#' calculated. Summary output therefore reports the median, SD, and empirical
#' interval of the transformed draw distribution itself.
#'
#' Fitted random effects are excluded. Residual, observation, process,
#' dispersion, and posterior-predictive noise are not added. Uncertainty
#' therefore refers to the fitted association/expected-response parameter
#' structure, not to a future noisy observation.
#'
#' @export
reduce_effects <- function(
    fit,
    data,
    vars,
    group = "epi_id",
    time = "time",
    type = c("overall", "lag", "var"),
    value = NULL,
    scale = c("percent", "response", "link"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    seed = NULL,
    ref = list(method = "median", value = NULL),
    at = NULL,
    interval_probs = c(0.025, 0.975),
    extrapolation = c("error", "warn", "allow")
) {

  # SMALL LOCAL VALIDATORS

  scalar_name <- function(x) {
    is.character(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      nzchar(x)
  }

  scalar_number <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x)
  }

  scalar_integer <- function(x) {
    scalar_number(x) &&
      x == as.integer(x)
  }

  scalar_flag <- function(x) {
    is.logical(x) &&
      length(x) == 1L &&
      !is.na(x)
  }

  # ARGUMENTS

  type <- match.arg(type)
  scale <- match.arg(scale)
  output <- match.arg(output)
  extrapolation <- match.arg(extrapolation)

  if (is.null(fit)) {
    stop(
      "`fit` cannot be NULL.",
      call. = FALSE
    )
  }

  if (!is.data.frame(data) ||
      !nrow(data)) {
    stop(
      "`data` must be a non-empty data.frame.",
      call. = FALSE
    )
  }

  if (!scalar_name(vars)) {
    stop(
      "`vars` must be one non-empty fitted exposure-variable name.",
      call. = FALSE
    )
  }

  if (!scalar_name(group)) {
    stop(
      "`group` must be one non-empty column name.",
      call. = FALSE
    )
  }

  if (!scalar_name(time)) {
    stop(
      "`time` must be one non-empty column name.",
      call. = FALSE
    )
  }

  if (identical(group, time)) {
    stop(
      "`group` and `time` must identify different columns.",
      call. = FALSE
    )
  }

  if (!scalar_flag(uncertainty)) {
    stop(
      "`uncertainty` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!scalar_integer(n_samples) ||
      n_samples < 1L) {
    stop(
      "`n_samples` must be one positive integer.",
      call. = FALSE
    )
  }
  n_samples <- as.integer(n_samples)

  if (uncertainty &&
      n_samples < 2L) {
    stop(
      "`n_samples` must be at least 2 when `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  if (identical(output, "samples") &&
      !uncertainty) {
    stop(
      "`output = 'samples'` requires `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    if (!scalar_integer(seed) ||
        seed < 0 ||
        seed > .Machine$integer.max) {
      stop(
        "`seed` must be NULL or one non-negative integer not exceeding ",
        "`.Machine$integer.max`.",
        call. = FALSE
      )
    }
    seed <- as.integer(seed)

    if (!uncertainty) {
      warning(
        "`seed` is ignored when `uncertainty = FALSE`.",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(interval_probs) ||
      length(interval_probs) != 2L ||
      anyNA(interval_probs) ||
      any(!is.finite(interval_probs)) ||
      any(interval_probs <= 0 |
          interval_probs >= 1) ||
      interval_probs[1L] >=
      interval_probs[2L]) {
    stop(
      "`interval_probs` must contain two finite probabilities satisfying ",
      "0 < interval_probs[1] < interval_probs[2] < 1.",
      call. = FALSE
    )
  }
  interval_probs <- as.numeric(interval_probs)

  # STRICT FITTED-MODEL CONTRACT

  required_internal_helpers <- c(
    ".get_epiexposure_metadata",
    ".extract_central_parameters",
    ".extract_parameter_draws",
    ".epix_cb_cols_for_var",
    ".epix_validate_regular_time"
  )

  function_environment <- environment()

  unavailable_helpers <- required_internal_helpers[
    !vapply(
      required_internal_helpers,
      function(helper) {
        exists(
          helper,
          envir = function_environment,
          mode = "function",
          inherits = TRUE
        )
      },
      logical(1)
    )
  ]

  if (length(unavailable_helpers)) {
    stop(
      "`reduce_effects()` requires the current EpiExposure internal ",
      "engine/prediction helpers. Missing helper(s): ",
      paste(unavailable_helpers, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  metadata <- .get_epiexposure_metadata(
    fit
  )

  if (!is.character(metadata$vars) ||
      !length(metadata$vars) ||
      anyNA(metadata$vars) ||
      any(!nzchar(metadata$vars)) ||
      anyDuplicated(metadata$vars)) {
    stop(
      "Invalid fitted exposure metadata in `fit`.",
      call. = FALSE
    )
  }

  if (!vars %in%
      metadata$vars) {
    stop(
      "Variable '",
      vars,
      "' was not fitted by the model.",
      call. = FALSE
    )
  }

  if (!is.list(metadata$spec) ||
      !all(metadata$vars %in%
           names(metadata$spec))) {
    stop(
      "The fitted model contains incomplete `epiexposure_spec` metadata.",
      call. = FALSE
    )
  }

  if (!is.list(metadata$basis_objects) ||
      !all(metadata$vars %in%
           names(metadata$basis_objects))) {
    stop(
      "`fit` does not contain one stored cross-basis object for every fitted ",
      "exposure. Refit the model with the current `fit_epidlnm()` ",
      "implementation.",
      call. = FALSE
    )
  }

  for (variable in metadata$vars) {
    if (!inherits(
      metadata$basis_objects[[variable]],
      "crossbasis"
    )) {
      stop(
        "Stored basis object for fitted exposure '",
        variable,
        "' must inherit from `crossbasis`.",
        call. = FALSE
      )
    }
  }

  family_name <- metadata$family_name
  link_name <- tolower(
    metadata$link
  )

  # `.get_epiexposure_metadata()` already validates supported v1 families, but
  # retain an explicit public-facing guard for the package-wide ordinal rule.
  if (identical(
    family_name,
    "ordinal"
  )) {
    stop(
      "Ordinal outcomes are not supported in EpiExposure v1.",
      call. = FALSE
    )
  }

  # EXACT COMMON TEMPORAL METADATA

  fit_max_lag <- attr(
    fit,
    "epiexposure_max_lag",
    exact = TRUE
  )

  fit_history_length <- attr(
    fit,
    "epiexposure_history_length",
    exact = TRUE
  )

  fit_history_contract <- attr(
    fit,
    "epiexposure_history_contract",
    exact = TRUE
  )

  fit_time_step <- attr(
    fit,
    "epiexposure_time_step",
    exact = TRUE
  )

  if (!scalar_integer(fit_max_lag) ||
      fit_max_lag < 0L) {
    stop(
      "`fit` is missing valid scalar `epiexposure_max_lag` metadata. Refit ",
      "with the current `fit_epidlnm()` implementation.",
      call. = FALSE
    )
  }
  fit_max_lag <- as.integer(fit_max_lag)

  if (!scalar_integer(fit_history_length) ||
      fit_history_length < 1L) {
    stop(
      "`fit` is missing valid scalar `epiexposure_history_length` metadata.",
      call. = FALSE
    )
  }
  fit_history_length <- as.integer(
    fit_history_length
  )

  if (!identical(
    fit_history_length,
    fit_max_lag + 1L
  )) {
    stop(
      "Inconsistent fitted temporal metadata: `history_length` must equal ",
      "`max_lag + 1`. Expected ",
      fit_max_lag + 1L,
      " but found ",
      fit_history_length,
      ".",
      call. = FALSE
    )
  }

  expected_history_contract <-
    "all_fitted_exposures_same_exact_max_lag_plus_one"

  if (!identical(
    fit_history_contract,
    expected_history_contract
  )) {
    stop(
      "`fit` does not satisfy the current EpiExposure exact common-profile ",
      "contract.",
      call. = FALSE
    )
  }

  for (variable in metadata$vars) {
    current_spec <- metadata$spec[[variable]]

    if (!is.list(current_spec) ||
        !scalar_integer(current_spec$max_lag) ||
        as.integer(current_spec$max_lag) !=
        fit_max_lag) {
      stop(
        "All fitted exposure variables must use the same `max_lag` under the ",
        "EpiExposure exact-profile contract. Invalid metadata detected for '",
        variable,
        "'.",
        call. = FALSE
      )
    }

    if (!is.null(current_spec$history_length)) {
      if (!scalar_integer(
        current_spec$history_length
      ) ||
      as.integer(
        current_spec$history_length
      ) !=
      fit_history_length) {
        stop(
          "Stored history-length metadata for fitted exposure '",
          variable,
          "' is inconsistent with the common fitted profile length of ",
          fit_history_length,
          ".",
          call. = FALSE
        )
      }
    }

    lag_attribute <- attr(
      metadata$basis_objects[[variable]],
      "lag",
      exact = TRUE
    )

    if (!is.numeric(lag_attribute) ||
        length(lag_attribute) != 2L ||
        anyNA(lag_attribute) ||
        any(!is.finite(lag_attribute)) ||
        !isTRUE(all.equal(
          as.numeric(lag_attribute),
          c(0, fit_max_lag),
          tolerance = 0
        ))) {
      stop(
        "Stored cross-basis lag metadata for fitted exposure '",
        variable,
        "' is inconsistent with the common fitted lag range 0:",
        fit_max_lag,
        ".",
        call. = FALSE
      )
    }
  }

  # INPUT DATA: ALL FITTED EXPOSURES + EXACT PROFILE

  required_columns <- unique(
    c(
      group,
      time,
      metadata$vars
    )
  )

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns)) {
    stop(
      "`data` is missing required fitted-model column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (anyNA(data[[group]])) {
    stop(
      "The grouping column '",
      group,
      "' cannot contain missing values.",
      call. = FALSE
    )
  }

  if (is.character(data[[group]]) &&
      any(!nzchar(data[[group]]))) {
    stop(
      "The grouping column '",
      group,
      "' cannot contain empty character identifiers.",
      call. = FALSE
    )
  }

  if (!is.numeric(data[[time]]) ||
      anyNA(data[[time]]) ||
      any(!is.finite(data[[time]]))) {
    stop(
      "`",
      time,
      "` must contain only finite numeric values.",
      call. = FALSE
    )
  }

  for (variable in metadata$vars) {
    if (!is.numeric(data[[variable]]) ||
        anyNA(data[[variable]]) ||
        any(!is.finite(data[[variable]]))) {
      stop(
        "Exposure variable '",
        variable,
        "' must contain only finite numeric values.",
        call. = FALSE
      )
    }
  }

  group_index <- as.character(
    data[[group]]
  )

  group_levels <- unique(
    group_index
  )

  if (!length(group_levels)) {
    stop(
      "`data` contains no exposure-profile groups.",
      call. = FALSE
    )
  }

  group_rows <- split(
    seq_len(nrow(data)),
    factor(
      group_index,
      levels = group_levels
    )
  )

  group_sizes <- vapply(
    group_rows,
    length,
    integer(1)
  )

  invalid_length <- group_sizes !=
    fit_history_length

  if (any(invalid_length)) {
    bad_group <- names(group_sizes)[
      which(invalid_length)[1L]
    ]
    bad_length <- group_sizes[
      which(invalid_length)[1L]
    ]

    stop(
      "Exposure profile for group '",
      bad_group,
      "' must contain exactly ",
      fit_history_length,
      " observations for max_lag = ",
      fit_max_lag,
      "; received ",
      bad_length,
      ". Profiles are never truncated, padded, or realigned silently.",
      call. = FALSE
    )
  }

  # Shared helper validates uniqueness, completeness of temporal spacing within
  # groups, and common spacing across groups.
  .epix_validate_regular_time(
    data = data,
    group_index = group_index,
    time_col = time
  )

  observed_steps <- numeric(0)

  if (fit_history_length >
      1L) {
    observed_steps <- vapply(
      group_rows,
      function(idx) {
        tt <- sort(
          data[[time]][idx]
        )
        diff(tt)[1L]
      },
      numeric(1)
    )
  }

  data_time_step <- if (
    length(observed_steps)
  ) {
    observed_steps[[1L]]
  } else {
    NA_real_
  }

  if (fit_history_length ==
      1L) {
    if (!is.null(fit_time_step) &&
        length(fit_time_step) == 1L &&
        is.finite(fit_time_step)) {
      stop(
        "For `max_lag = 0`, the fitted profile contains one observation and ",
        "`epiexposure_time_step` must be `NA` because no temporal interval can ",
        "be inferred.",
        call. = FALSE
      )
    }
  } else {
    if (!scalar_number(fit_time_step) ||
        fit_time_step <= 0) {
      stop(
        "`fit` contains invalid `epiexposure_time_step` metadata.",
        call. = FALSE
      )
    }

    step_tolerance <- sqrt(
      .Machine$double.eps
    ) *
      max(
        1,
        abs(
          fit_time_step
        )
      )

    if (abs(
      data_time_step -
      fit_time_step
    ) >
    step_tolerance) {
      stop(
        "`data` uses temporal spacing ",
        format(data_time_step),
        " but the fitted model uses spacing ",
        format(fit_time_step),
        ". Reduction data must use the same temporal support as the fitted ",
        "DLNM.",
        call. = FALSE
      )
    }
  }

  # TYPE / VALUE / AT

  if (identical(
    type,
    "overall"
  )) {
    if (!is.null(value)) {
      stop(
        "`value` must be NULL when `type = 'overall'`.",
        call. = FALSE
      )
    }
  } else {
    if (!scalar_number(value)) {
      stop(
        "`value` must be one finite numeric scalar when `type = '",
        type,
        "'`.",
        call. = FALSE
      )
    }
    value <- as.numeric(value)
  }

  if (identical(
    type,
    "lag"
  ) &&
  (
    value < 0 ||
    value >
    fit_max_lag
  )) {
    stop(
      "For `type = 'lag'`, `value` must lie within the fitted lag range 0:",
      fit_max_lag,
      ".",
      call. = FALSE
    )
  }

  if (identical(
    type,
    "var"
  )) {
    if (!is.null(at)) {
      stop(
        "`at` is not used when `type = 'var'`; the reduced curve is evaluated ",
        "over lag.",
        call. = FALSE
      )
    }
  } else if (!is.null(at)) {
    if (!is.numeric(at) ||
        !length(at) ||
        anyNA(at) ||
        any(!is.finite(at))) {
      stop(
        "`at` must be NULL or a non-empty finite numeric vector.",
        call. = FALSE
      )
    }

    if (anyDuplicated(at)) {
      stop(
        "`at` must contain unique exposure values.",
        call. = FALSE
      )
    }

    at <- as.numeric(at)
  }

  # STORED FOCAL CROSS-BASIS AND NAME-BASED COEFFICIENT MAPPING

  basis <- metadata$basis_objects[[vars]]

  cb_names <- .epix_cb_cols_for_var(
    metadata$cb_cols,
    vars
  )

  if (!length(cb_names)) {
    stop(
      "No canonical fitted cross-basis coefficient names were found for ",
      "variable '",
      vars,
      "'.",
      call. = FALSE
    )
  }

  canonical_basis_names <- attr(
    basis,
    "epiexposure_canonical_cb_colnames",
    exact = TRUE
  )

  if (is.null(canonical_basis_names) ||
      !is.character(canonical_basis_names) ||
      !identical(
        canonical_basis_names,
        cb_names
      )) {
    stop(
      "Stored canonical cross-basis column metadata for variable '",
      vars,
      "' is missing or inconsistent. Coefficients are not aligned by ",
      "position. Refit with the current EpiExposure implementation.",
      call. = FALSE
    )
  }

  if (ncol(basis) !=
      length(cb_names)) {
    stop(
      "Stored cross-basis dimension for variable '",
      vars,
      "' does not match its canonical fitted coefficient block.",
      call. = FALSE
    )
  }

  native_basis_names <- colnames(
    basis
  )

  if (is.null(native_basis_names) ||
      length(native_basis_names) !=
      ncol(basis) ||
      anyNA(native_basis_names) ||
      any(!nzchar(native_basis_names)) ||
      anyDuplicated(native_basis_names)) {
    stop(
      "Stored cross-basis for variable '",
      vars,
      "' must contain unique native dlnm column names.",
      call. = FALSE
    )
  }

  argvar <- attr(
    basis,
    "argvar",
    exact = TRUE
  )

  arglag <- attr(
    basis,
    "arglag",
    exact = TRUE
  )

  if (!is.list(argvar) ||
      is.null(argvar$fun) ||
      !is.list(arglag) ||
      is.null(arglag$fun)) {
    stop(
      "Stored cross-basis for variable '",
      vars,
      "' lacks valid effective `argvar`/`arglag` metadata.",
      call. = FALSE
    )
  }

  if (isTRUE(
    argvar$intercept
  )) {
    stop(
      "The stored exposure basis for '",
      vars,
      "' includes an intercept, which violates the EpiExposure v1 ",
      "cross-basis identifiability contract.",
      call. = FALSE
    )
  }

  basis_range <- attr(
    basis,
    "range",
    exact = TRUE
  )

  if (!is.numeric(basis_range) ||
      length(basis_range) != 2L ||
      anyNA(basis_range) ||
      any(!is.finite(basis_range)) ||
      basis_range[1L] >
      basis_range[2L]) {
    stop(
      "Stored exposure range for fitted variable '",
      vars,
      "' is missing or invalid.",
      call. = FALSE
    )
  }

  basis_range <- as.numeric(
    basis_range
  )

  # REFERENCE

  if (!is.list(ref) ||
      !length(ref) ||
      is.null(names(ref)) ||
      anyNA(names(ref)) ||
      any(!nzchar(names(ref))) ||
      anyDuplicated(names(ref))) {
    stop(
      "`ref` must be a non-empty named list such as ",
      "`list(method = 'median', value = NULL)`.",
      call. = FALSE
    )
  }

  unknown_ref_names <- setdiff(
    names(ref),
    c(
      "method",
      "value"
    )
  )

  if (length(unknown_ref_names)) {
    stop(
      "`ref` contains unsupported field(s): ",
      paste(
        unknown_ref_names,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (!"method" %in%
      names(ref) ||
      !scalar_name(ref$method)) {
    stop(
      "`ref$method` must be one non-empty character value.",
      call. = FALSE
    )
  }

  ref_method <- match.arg(
    ref$method,
    choices = c(
      "median",
      "percentile",
      "fixed"
    )
  )

  focal_exposure <- as.numeric(
    data[[vars]]
  )

  if (identical(
    ref_method,
    "median"
  )) {
    if (!is.null(ref$value)) {
      stop(
        "`ref$value` must be NULL when `ref$method = 'median'`.",
        call. = FALSE
      )
    }

    reference_value <- stats::median(
      focal_exposure
    )

  } else if (identical(
    ref_method,
    "percentile"
  )) {
    if (!scalar_number(ref$value) ||
        ref$value < 0 ||
        ref$value > 1) {
      stop(
        "For `ref$method = 'percentile'`, `ref$value` must be one finite ",
        "probability between 0 and 1.",
        call. = FALSE
      )
    }

    reference_value <- as.numeric(
      stats::quantile(
        focal_exposure,
        probs = ref$value,
        names = FALSE,
        type = 7
      )
    )

  } else {
    if (!scalar_number(
      ref$value
    )) {
      stop(
        "For `ref$method = 'fixed'`, `ref$value` must be one finite numeric ",
        "value.",
        call. = FALSE
      )
    }

    reference_value <- as.numeric(
      ref$value
    )
  }

  if (!is.finite(
    reference_value
  )) {
    stop(
      "Could not resolve a finite reference value for variable '",
      vars,
      "'.",
      call. = FALSE
    )
  }

  # EXPOSURE-SUPPORT POLICY

  handle_exposure_support <- function(
    values,
    label
  ) {

    if (is.null(values) ||
        !length(values)) {
      return(
        invisible(TRUE)
      )
    }

    outside <- values <
      basis_range[1L] |
      values >
      basis_range[2L]

    if (!any(outside)) {
      return(
        invisible(TRUE)
      )
    }

    offending <- unique(
      as.numeric(
        values[
          outside
        ]
      )
    )

    message_text <- paste0(
      label,
      " for variable '",
      vars,
      "' extends outside the fitted exposure range [",
      format(basis_range[1L]),
      ", ",
      format(basis_range[2L]),
      "]. Offending value(s): ",
      paste(
        utils::head(
          offending,
          8L
        ),
        collapse = ", "
      ),
      if (length(offending) >
          8L) {
        ", ..."
      } else {
        ""
      },
      "."
    )

    if (identical(
      extrapolation,
      "error"
    )) {
      stop(
        message_text,
        call. = FALSE
      )
    }

    if (identical(
      extrapolation,
      "warn"
    )) {
      warning(
        message_text,
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  handle_exposure_support(
    reference_value,
    "Reference value"
  )

  if (identical(
    type,
    "var"
  )) {
    handle_exposure_support(
      value,
      "`value`"
    )
  }

  if (!is.null(at)) {
    handle_exposure_support(
      at,
      "`at`"
    )
  }

  # SCALE CONTRACT

  if (identical(
    scale,
    "percent"
  ) &&
  !link_name %in%
  c(
    "log",
    "logit"
  )) {
    stop(
      "`scale = 'percent'` is supported only for fitted log and logit links. ",
      "For a log link it represents percent relative change in expected ",
      "response; for a logit link it represents percent change in odds.",
      call. = FALSE
    )
  }

  if (identical(
    scale,
    "response"
  ) &&
  !link_name %in%
  c(
    "identity",
    "log",
    "logit"
  )) {
    stop(
      "`scale = 'response'` cannot be derived from a centered DLNM contrast ",
      "alone for fitted link '",
      link_name,
      "'. Use `scale = 'link'`, or use `summarise_effects()` when absolute ",
      "response-scale baseline/predicted quantities are required.",
      call. = FALSE
    )
  }

  transform_effect <- function(x) {

    if (identical(
      scale,
      "link"
    )) {
      return(
        x
      )
    }

    if (identical(
      scale,
      "response"
    )) {
      if (identical(
        link_name,
        "identity"
      )) {
        return(
          x
        )
      }

      return(
        exp(x)
      )
    }

    # `scale = "percent"` reaches here only for log/logit links.
    100 *
      (
        exp(x) -
          1
      )
  }

  # CENTRAL COEFFICIENTS

  central_parameters <- .extract_central_parameters(
    fit
  )

  missing_central <- setdiff(
    cb_names,
    names(
      central_parameters
    )
  )

  if (length(missing_central)) {
    stop(
      "Central parameter vector is missing cross-basis coefficient(s) for '",
      vars,
      "': ",
      paste(
        missing_central,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  beta_sub <- central_parameters[
    cb_names
  ]

  if (length(beta_sub) !=
      ncol(basis) ||
      anyNA(beta_sub) ||
      any(!is.finite(beta_sub))) {
    stop(
      "Central cross-basis coefficient block for variable '",
      vars,
      "' is invalid.",
      call. = FALSE
    )
  }

  # Map the already name-aligned canonical EpiExposure coefficient block to
  # the stored native dlnm cross-basis names. The ordering is verified above by
  # `epiexposure_canonical_cb_colnames`; it is never inferred from coefficient
  # position alone.
  beta_dlnm <- as.numeric(
    beta_sub
  )
  names(beta_dlnm) <- native_basis_names

  zero_vcov <- matrix(
    0,
    nrow = ncol(basis),
    ncol = ncol(basis),
    dimnames = list(
      native_basis_names,
      native_basis_names
    )
  )

  # CROSSREDUCE WRAPPER

  if (!requireNamespace(
    "dlnm",
    quietly = TRUE
  )) {
    stop(
      "Package 'dlnm' is required by `reduce_effects()`.",
      call. = FALSE
    )
  }

  run_crossreduce <- function(
    coefficients
  ) {

    coefficients <- as.numeric(
      coefficients
    )

    if (length(coefficients) !=
        ncol(basis) ||
        anyNA(coefficients) ||
        any(!is.finite(coefficients))) {
      stop(
        "Cross-basis coefficient vector supplied to the reduction is invalid.",
        call. = FALSE
      )
    }

    names(coefficients) <-
      native_basis_names

    args <- list(
      basis = basis,
      coef = coefficients,
      vcov = zero_vcov,
      model.link = link_name,
      type = type,
      cen = reference_value
    )

    if (type %in%
        c(
          "lag",
          "var"
        )) {
      args$value <- value
    }

    if (!is.null(at) &&
        !identical(
          type,
          "var"
        )) {
      args$at <- at
    }

    out <- tryCatch(
      do.call(
        dlnm::crossreduce,
        args
      ),
      error = function(e) {
        stop(
          "`dlnm::crossreduce()` failed for fitted exposure '",
          vars,
          "': ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (!inherits(
      out,
      "crossreduce"
    )) {
      stop(
        "`dlnm::crossreduce()` did not return a valid `crossreduce` object.",
        call. = FALSE
      )
    }

    out
  }

  point_reduction <- run_crossreduce(
    beta_dlnm
  )

  reduction_coordinate <- function(
    reduction
  ) {

    if (identical(
      type,
      "var"
    )) {
      lag_range <- reduction$lag

      if (!is.numeric(lag_range) ||
          length(lag_range) != 2L ||
          anyNA(lag_range) ||
          any(!is.finite(lag_range)) ||
          !scalar_number(
            reduction$bylag
          ) ||
          reduction$bylag <= 0) {
        stop(
          "The reduced lag-response object contains invalid lag coordinates.",
          call. = FALSE
        )
      }

      x <- seq(
        from = lag_range[1L],
        to = lag_range[2L],
        by = reduction$bylag
      )

    } else {
      x <- reduction$predvar
    }

    if (!is.numeric(x) ||
        !length(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "The reduced DLNM contains invalid prediction coordinates.",
        call. = FALSE
      )
    }

    as.numeric(x)
  }

  x_point <- reduction_coordinate(
    point_reduction
  )

  eta_point <- as.numeric(
    point_reduction$fit
  )

  if (length(eta_point) !=
      length(x_point) ||
      anyNA(eta_point) ||
      any(!is.finite(eta_point))) {
    stop(
      "The central DLNM reduction produced invalid association values.",
      call. = FALSE
    )
  }

  effect_point <- as.numeric(
    transform_effect(
      eta_point
    )
  )

  if (length(effect_point) !=
      length(eta_point) ||
      anyNA(effect_point) ||
      any(!is.finite(effect_point))) {
    stop(
      "The requested effect transformation produced non-finite central values.",
      call. = FALSE
    )
  }

  # DETERMINISTIC RETURN

  if (!uncertainty) {

    result <- data.frame(
      x = x_point,
      eta = eta_point,
      effect = effect_point,
      stringsAsFactors = FALSE
    )

  } else {

    # REPRODUCIBLE JOINT PARAMETER DRAWS

    had_random_seed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )

    if (had_random_seed) {
      previous_random_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    } else {
      previous_random_seed <- NULL
    }

    if (!is.null(seed)) {
      on.exit(
        {
          if (had_random_seed) {
            assign(
              ".Random.seed",
              previous_random_seed,
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

      set.seed(
        seed
      )
    }

    parameter_draws <- .extract_parameter_draws(
      fit = fit,
      n_samples = n_samples
    )

    if (!is.matrix(parameter_draws) ||
        nrow(parameter_draws) < 2L ||
        is.null(colnames(parameter_draws)) ||
        anyNA(parameter_draws) ||
        any(!is.finite(parameter_draws))) {
      stop(
        "The shared EpiExposure parameter-draw helper returned an invalid ",
        "draw matrix.",
        call. = FALSE
      )
    }

    missing_draw_columns <- setdiff(
      cb_names,
      colnames(
        parameter_draws
      )
    )

    if (length(missing_draw_columns)) {
      stop(
        "Parameter-draw matrix is missing cross-basis coefficient(s) for '",
        vars,
        "': ",
        paste(
          missing_draw_columns,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    focal_draws <- parameter_draws[
      ,
      cb_names,
      drop = FALSE
    ]

    n_draws <- nrow(
      focal_draws
    )

    eta_draw_list <- vector(
      "list",
      n_draws
    )

    for (i in seq_len(
      n_draws
    )) {

      draw_reduction <- run_crossreduce(
        focal_draws[
          i,
          ,
          drop = TRUE
        ]
      )

      x_draw <- reduction_coordinate(
        draw_reduction
      )

      if (length(x_draw) !=
          length(x_point) ||
          !isTRUE(all.equal(
            x_draw,
            x_point,
            tolerance = 0
          ))) {
        stop(
          "The reduction coordinate changed across parameter draws. This ",
          "violates the fixed fitted-basis contract.",
          call. = FALSE
        )
      }

      current_eta <- as.numeric(
        draw_reduction$fit
      )

      if (length(current_eta) !=
          length(x_point) ||
          anyNA(current_eta) ||
          any(!is.finite(current_eta))) {
        stop(
          "A draw-specific DLNM reduction produced non-finite association ",
          "values.",
          call. = FALSE
        )
      }

      eta_draw_list[[i]] <-
        current_eta
    }

    # Rows are reduction coordinates; columns are coherent joint parameter
    # draws.
    eta_matrix <- do.call(
      cbind,
      eta_draw_list
    )

    if (!is.matrix(eta_matrix)) {
      eta_matrix <- matrix(
        eta_matrix,
        nrow = length(x_point),
        ncol = n_draws
      )
    }

    if (!identical(
      dim(eta_matrix),
      c(
        length(x_point),
        n_draws
      )
    ) ||
    anyNA(eta_matrix) ||
    any(!is.finite(eta_matrix))) {
      stop(
        "Could not construct a valid draw-by-reduction-coordinate matrix.",
        call. = FALSE
      )
    }

    effect_matrix <- matrix(
      as.numeric(
        transform_effect(
          as.vector(
            eta_matrix
          )
        )
      ),
      nrow = nrow(eta_matrix),
      ncol = ncol(eta_matrix),
      byrow = FALSE
    )

    if (anyNA(effect_matrix) ||
        any(!is.finite(effect_matrix))) {
      stop(
        "Draw-by-draw effect transformation produced non-finite values.",
        call. = FALSE
      )
    }

    # DRAW-LEVEL OUTPUT

    if (identical(
      output,
      "samples"
    )) {

      result <- data.frame(
        x = rep(
          x_point,
          times = n_draws
        ),
        eta = as.vector(
          eta_matrix
        ),
        effect = as.vector(
          effect_matrix
        ),
        sample = rep(
          seq_len(n_draws),
          each = length(x_point)
        ),
        stringsAsFactors = FALSE
      )

    } else {

      # EMPIRICAL SUMMARY OUTPUT

      eta_median <- apply(
        eta_matrix,
        1L,
        stats::median
      )

      eta_sd <- apply(
        eta_matrix,
        1L,
        stats::sd
      )

      eta_low <- apply(
        eta_matrix,
        1L,
        stats::quantile,
        probs = interval_probs[1L],
        names = FALSE,
        type = 7
      )

      eta_high <- apply(
        eta_matrix,
        1L,
        stats::quantile,
        probs = interval_probs[2L],
        names = FALSE,
        type = 7
      )

      effect_median <- apply(
        effect_matrix,
        1L,
        stats::median
      )

      effect_sd <- apply(
        effect_matrix,
        1L,
        stats::sd
      )

      effect_low <- apply(
        effect_matrix,
        1L,
        stats::quantile,
        probs = interval_probs[1L],
        names = FALSE,
        type = 7
      )

      effect_high <- apply(
        effect_matrix,
        1L,
        stats::quantile,
        probs = interval_probs[2L],
        names = FALSE,
        type = 7
      )

      result <- data.frame(
        x = x_point,
        eta = as.numeric(
          eta_median
        ),
        eta_sd = as.numeric(
          eta_sd
        ),
        low = as.numeric(
          eta_low
        ),
        high = as.numeric(
          eta_high
        ),
        effect = as.numeric(
          effect_median
        ),
        effect_sd = as.numeric(
          effect_sd
        ),
        low_eff = as.numeric(
          effect_low
        ),
        high_eff = as.numeric(
          effect_high
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  # ROW-LEVEL METADATA

  result$type <- rep(
    type,
    nrow(result)
  )

  result$value <- rep(
    if (is.null(value)) {
      NA_real_
    } else {
      as.numeric(value)
    },
    nrow(result)
  )

  result$reference <- rep(
    reference_value,
    nrow(result)
  )

  result$scale <- rep(
    scale,
    nrow(result)
  )

  result$var <- rep(
    vars,
    nrow(result)
  )

  # Keep the most interpretable core columns first.
  preferred_order <- c(
    "var",
    "type",
    "value",
    "reference",
    "scale",
    "x",
    "eta",
    "eta_sd",
    "low",
    "high",
    "effect",
    "effect_sd",
    "low_eff",
    "high_eff",
    "sample"
  )

  preferred_order <- preferred_order[
    preferred_order %in%
      names(result)
  ]

  result <- result[
    ,
    c(
      preferred_order,
      setdiff(
        names(result),
        preferred_order
      )
    ),
    drop = FALSE
  ]

  rownames(result) <- NULL

  # OUTPUT ATTRIBUTES

  attr(
    result,
    "epiexposure_reduction_contract"
  ) <- "stored_fitted_crossbasis_crossreduce_no_refit"

  attr(
    result,
    "epiexposure_reduction_type"
  ) <- type

  attr(
    result,
    "epiexposure_reduction_variable"
  ) <- vars

  attr(
    result,
    "epiexposure_reduction_coordinate"
  ) <- if (
    identical(
      type,
      "var"
    )
  ) {
    "lag"
  } else {
    "exposure"
  }

  attr(
    result,
    "epiexposure_reference_value"
  ) <- reference_value

  attr(
    result,
    "epiexposure_reference_method"
  ) <- ref_method

  attr(
    result,
    "epiexposure_effect_scale"
  ) <- scale

  attr(
    result,
    "epiexposure_effect_contract"
  ) <- switch(
    scale,
    link = "centered_dlnm_linear_predictor_contrast",
    response = if (
      identical(
        link_name,
        "identity"
      )
    ) {
      "additive_response_difference_from_centered_contrast"
    } else if (
      identical(
        link_name,
        "log"
      )
    ) {
      "exponentiated_contrast_response_ratio"
    } else {
      "exponentiated_contrast_odds_ratio"
    },
    percent = if (
      identical(
        link_name,
        "log"
      )
    ) {
      "percent_relative_change_in_expected_response"
    } else {
      "percent_change_in_odds"
    }
  )

  attr(
    result,
    "epiexposure_reduced_coefficients"
  ) <- point_reduction$coefficients

  attr(
    result,
    "epiexposure_reduced_basis"
  ) <- point_reduction$basis

  attr(
    result,
    "epiexposure_uncertainty"
  ) <- uncertainty

  attr(
    result,
    "epiexposure_uncertainty_contract"
  ) <- if (uncertainty) {
    "joint_fixed_population_parameter_draws_no_residual_process_noise"
  } else {
    "central_parameter_estimate"
  }

  attr(
    result,
    "epiexposure_summary_center"
  ) <- if (
    uncertainty &&
    identical(
      output,
      "summary"
    )
  ) {
    "median"
  } else {
    NULL
  }

  attr(
    result,
    "epiexposure_interval_probs"
  ) <- if (uncertainty) {
    interval_probs
  } else {
    NULL
  }

  attr(
    result,
    "epiexposure_n_samples"
  ) <- if (uncertainty) {
    if (exists(
      "n_draws",
      inherits = FALSE
    )) {
      n_draws
    } else {
      n_samples
    }
  } else {
    NULL
  }

  attr(
    result,
    "epiexposure_prediction_level"
  ) <- "population"

  attr(
    result,
    "epiexposure_prediction_estimand"
  ) <- "reduced_dlnm_association_contrast"

  attr(
    result,
    "epiexposure_max_lag"
  ) <- fit_max_lag

  attr(
    result,
    "epiexposure_history_length"
  ) <- fit_history_length

  attr(
    result,
    "epiexposure_history_contract"
  ) <- fit_history_contract

  attr(
    result,
    "epiexposure_time_step"
  ) <- data_time_step

  attr(
    result,
    "epiexposure_profile_order"
  ) <- "chronological"

  attr(
    result,
    "epiexposure_extrapolation"
  ) <- extrapolation

  result
}
