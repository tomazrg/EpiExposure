#' Compute Exposure Cumulative Impact
#'
#' Computes a descriptive raw cumulative exposure (`ECI_raw`) and, when a
#' fitted EpiExposure model is supplied, a model-weighted cumulative exposure
#' impact (`ECI_weighted`) for one or more exposure variables.
#'
#' Exposure histories are interpreted chronologically, from the oldest
#' observation to the most recent observation. The most recent value is
#' associated internally with lag 0.
#'
#' The model-weighted ECI is a **contrast relative to an explicit joint exposure
#' reference profile**. It is not the uncentered cross-basis contribution
#' `cb(x) %*% beta`, because that quantity depends on basis parameterization and
#' is not, by itself, an interpretable exposure effect.
#'
#' @param profile Optional exposure profile input. Supply exactly one of
#'   `profile` or `data`.
#'
#'   When one exposure variable is evaluated, `profile` may be a finite numeric
#'   vector. For multiple variables, supply a list with one numeric vector per
#'   requested variable. A named list is recommended. If the list is unnamed,
#'   its order must match `var`.
#'
#'   Profiles must be chronological from oldest to most recent. When `fit` is
#'   supplied, every focal exposure profile must contain **exactly**
#'   `max_lag + 1` observations for that fitted exposure. Histories that are
#'   shorter or longer are rejected; `compute_eci()` never selects a temporal
#'   window silently from a longer profile.
#'
#'   When `fit = NULL`, no fitted `max_lag` exists. In that raw-only mode,
#'   `profile` must be one finite numeric vector and the function simply returns
#'   its descriptive cumulative sum.
#' @param data Optional non-empty long-format data frame containing observed
#'   chronological exposure histories. Supply exactly one of `profile` or
#'   `data`.
#'
#'   `data` must contain the grouping column named by `group`, the time column
#'   named by `time`, and every exposure requested in `var`. Every evaluated
#'   group must contain exactly `max_lag + 1` equally spaced observations for
#'   the requested focal exposure(s). Groups with either fewer or more rows are
#'   rejected. If a method-based reference is used, `data` must additionally
#'   contain every exposure in the fitted model because those variables are
#'   needed to construct the joint reference profile.
#' @param group Character scalar naming the independent-history grouping column
#'   in `data`. Required when `data` is supplied.
#' @param time Character scalar naming the chronological numeric time column in
#'   `data`. Default is `"time"`. Times must be unique and equally spaced within
#'   each evaluated group, and the temporal spacing must be common across
#'   groups.
#' @param group_level Optional grouping value or vector of values to evaluate.
#'   If `NULL`, all group levels are evaluated in first-occurrence order.
#' @param fit Optional fitted model returned by the current `fit_epidlnm()`.
#'   If `NULL`, only `ECI_raw` is calculated and model-weighted ECI,
#'   uncertainty, and response-scale quantities are unavailable.
#' @param var Optional character vector naming fitted exposure variables to
#'   evaluate.
#'
#'   If `fit` contains one exposure and `var = NULL`, that exposure is used.
#'   For a multivariable fit, `var` must be supplied explicitly. Each requested
#'   variable is evaluated separately against the same joint reference profile;
#'   the function does not sum impacts from different exposure variables.
#' @param ref Reference exposure specification used when `fit` is supplied.
#'
#'   Two forms are supported:
#'
#'   \itemize{
#'     \item Method-based:
#'       `list(method = "median", value = NULL)`,
#'       `list(method = "percentile", value = p)`, or
#'       `list(method = "fixed", value = x)`.
#'     \item Exposure-specific: a named list containing exactly one finite
#'       numeric reference value for **every fitted exposure**, for example
#'       `list(tmean = 25, rain = 5, wetness = 10)`.
#'   }
#'
#'   Method-based references require `data` because the reference values are
#'   derived from the supplied exposure data. `"median"` uses each fitted
#'   exposure's own median. `"percentile"` uses the same percentile `p` for
#'   every fitted exposure. `"fixed"` applies the same numeric value to all
#'   fitted exposures and should be used cautiously when variables have
#'   different units.
#'
#'   For direct `profile` input, use the exposure-specific named-list form.
#'
#'   The default is `list(method = "median", value = NULL)`.
#' @param scale Character defining the interpretation and units of the single
#'   returned model-weighted impact column, `ECI_weighted`:
#'
#'   - `"link"`: the centered cumulative DLNM contrast
#'     \eqn{\Delta\eta};
#'   - `"response"`: the absolute response-scale impact
#'     \eqn{P - B}, where `P` is the population expected response for the focal
#'     profile and `B` is the population expected response for the joint
#'     reference profile;
#'   - `"percent"`: `100 * (exp(Delta eta) - 1)`, available only for log and
#'     logit links.
#'
#'   Internally, `compute_eci()` calculates both the centered link-scale
#'   contrast \eqn{\Delta\eta} and the response-scale difference
#'   \eqn{\Delta = P-B}. These internal quantities are intentionally **not
#'   returned as separate `eta` and `delta` columns** because they duplicate
#'   `ECI_weighted` on the corresponding scale: with `scale = "link"`,
#'   `ECI_weighted = Delta eta`; with `scale = "response"`,
#'   `ECI_weighted = P - B`. The internal quantities are still retained during
#'   computation because they are required for scale transformations and
#'   draw-by-draw uncertainty propagation.
#'
#'   For a log link, `"percent"` is the percent relative change in the expected
#'   response. For a logit link, it is the percent change in odds, not the
#'   percentage-point change in probability.
#' @param uncertainty Logical. If `FALSE`, use the harmonized central
#'   fixed/population parameter estimate. If `TRUE`, propagate joint
#'   fixed/population parameter uncertainty draw by draw.
#' @param output Character. `"summary"` returns deterministic values when
#'   `uncertainty = FALSE`, or medians, SDs, and empirical intervals when
#'   `uncertainty = TRUE`. `"samples"` returns one row per parameter draw and
#'   requires `uncertainty = TRUE`.
#' @param n_samples Positive integer number of parameter draws used when
#'   `uncertainty = TRUE`. At least two are required. Default is 1000.
#' @param interval_probs Numeric vector of length two defining the empirical
#'   uncertainty interval. Default `c(0.025, 0.975)` gives a 95 percent
#'   interval.
#' @param seed Optional finite integer used for parameter-draw sampling. The
#'   caller's global random-number state is restored when the function exits.
#' @param extrapolation Character controlling exposure values outside the
#'   fitted cross-basis exposure range: `"error"` (default), `"warn"`, or
#'   `"allow"`. The fitted basis is never re-estimated from ECI input data.
#'
#' @return A data frame.
#'
#'   With a fitted model, principal columns are:
#'
#'   \describe{
#'     \item{`var`}{Focal exposure variable.}
#'     \item{`reference_value`}{Reference exposure value for the focal
#'       variable.}
#'     \item{`max_lag`}{Common fitted maximum lag.}
#'     \item{`n_exposure_values`}{Number of exposure values used; always
#'       `max_lag + 1`.}
#'     \item{`ECI_raw`}{Descriptive sum of the focal exposure values over the
#'       fitted lag window.}
#'     \item{`ECI_raw_centered`}{Descriptive sum of exposure deviations from the
#'       focal reference value over the same lag window.}
#'     \item{`baseline`}{Population-level expected response under the joint
#'       reference exposure profile, with fitted random effects excluded.}
#'     \item{`predicted`}{Population-level expected response when the focal
#'       exposure follows the evaluated profile and every other fitted exposure
#'       remains at its reference value.}
#'     \item{`ECI_weighted`}{The model-weighted exposure impact on the scale
#'       requested by `scale`. It equals the internally calculated
#'       \eqn{\Delta\eta} on the link scale, the internally calculated
#'       \eqn{P-B} on the response scale, or the permitted exponential percent
#'       transformation on the percent scale.}
#'     \item{`scale`}{Requested weighted-ECI scale.}
#'   }
#'
#'   `eta` and `delta` are not returned as separate columns. Their information is
#'   represented by `ECI_weighted` when `scale = "link"` and
#'   `scale = "response"`, respectively.
#'
#'   With `uncertainty = TRUE` and `output = "summary"`, `baseline`,
#'   `predicted`, and `ECI_weighted` are each summarized by their median plus
#'   `_sd`, `_lower`, and `_upper` columns. `ECI_raw` and
#'   `ECI_raw_centered` are exposure-profile descriptors and therefore remain
#'   deterministic.
#'
#'   With `uncertainty = TRUE` and `output = "samples"`, `baseline`,
#'   `predicted`, and `ECI_weighted` are returned draw by draw with a `sample`
#'   column. The internal link-scale and response-scale contrasts are calculated
#'   for every draw but are not duplicated as output columns.
#'
#'   When `fit = NULL`, the result contains only `ECI_raw` and the number of
#'   exposure values used.
#'
#'   The output attribute `epiexposure_eci_history_contract` records
#'   `"exact_max_lag_plus_one"` whenever a fitted model is used.
#'
#' @details
#' ## Raw ECI
#'
#' For the discrete exposure history used by the fitted lag window,
#'
#' \deqn{
#'   ECI_{raw} = \sum_{j=0}^{L} x_j.
#' }
#'
#' This is a discrete cumulative exposure sum, not a continuous-time integral.
#' It retains the units of the exposure multiplied by the number of discrete
#' lag positions. Comparisons are meaningful only for the same exposure
#' variable, units, lag window, and temporal resolution.
#'
#' The centered descriptive counterpart is
#'
#' \deqn{
#'   ECI_{raw,centered} =
#'   \sum_{j=0}^{L}(x_j-x_{ref}).
#' }
#'
#' Neither raw quantity incorporates the fitted exposure-response or lag-response
#' shape.
#'
#' ## Exact exposure-history length
#'
#' When a fitted model is supplied, `compute_eci()` uses the EpiExposure exact
#' temporal-history contract. For a fitted maximum lag `L`, every evaluated
#' focal exposure history must contain exactly
#'
#' \deqn{
#'   L + 1
#' }
#'
#' chronological observations: the first value represents lag `L` and the last
#' value represents lag 0. All fitted exposures share the same `max_lag`.
#' Longer histories are not truncated with `tail()` and shorter histories are
#' not padded. This prevents ambiguity about which exposure window generated the
#' cumulative impact and keeps ECI calculations comparable across profiles,
#' groups, and downstream lag decompositions.
#'
#' In raw-only mode (`fit = NULL`), this restriction cannot be inferred because
#' no fitted `max_lag` is available; the supplied vector is summed as given.
#'
#' ## Centered weighted ECI
#'
#' A DLNM cross-basis is a basis expansion. The uncentered quantity
#'
#' \deqn{
#'   cb(x)^T\beta
#' }
#'
#' is a contribution to the linear predictor under the chosen basis
#' parameterization, but it is not an invariant exposure effect. In particular,
#' applying the inverse link to this contribution alone does not produce a valid
#' expected response because the model intercept and the reference contribution
#' of the other fitted exposures are absent.
#'
#' EpiExposure therefore defines the model-weighted cumulative impact as a
#' contrast between the focal profile and an explicit joint reference profile:
#'
#' \deqn{
#'   \Delta\eta =
#'   \{cb(x)-cb(x_{ref})\}^T\beta_{focal}.
#' }
#'
#' Every non-focal fitted exposure is held at its own reference value throughout
#' its complete lag window.
#'
#' This centering makes `ECI_weighted` consistent with the DLNM contrast logic
#' used by `summarise_effects()`.
#'
#' ## Why only `ECI_weighted` is returned
#'
#' Three mathematical quantities are required internally:
#'
#' \deqn{
#'   \Delta\eta = \eta_{target} - \eta_{reference},
#' }
#'
#' \deqn{
#'   B = g^{-1}(\eta_{reference}), \qquad
#'   P = g^{-1}(\eta_{target}),
#' }
#'
#' and
#'
#' \deqn{
#'   \Delta = P-B.
#' }
#'
#' Earlier output designs could expose `eta`, `delta`, and `ECI_weighted`
#' simultaneously. This is redundant because the selected scale already defines
#' which model-weighted impact is being reported. The current output therefore
#' uses one canonical column:
#'
#' \itemize{
#'   \item `scale = "link"`: `ECI_weighted` is \eqn{\Delta\eta};
#'   \item `scale = "response"`: `ECI_weighted` is \eqn{\Delta = P-B};
#'   \item `scale = "percent"`: `ECI_weighted` is the permitted percent
#'     transformation of \eqn{\Delta\eta}.
#' }
#'
#' `baseline` and `predicted` remain in the output because they are distinct
#' expected-response quantities that provide the response-scale context for the
#' contrast. The intermediate `eta` and `delta` quantities continue to be
#' computed internally, including for every uncertainty draw, but are not
#' duplicated as public output columns.
#'
#' ## Response-scale ECI
#'
#' Let \eqn{X_{ref}} denote the complete fixed/population design under the
#' joint reference exposure profile. The population expected-response baseline
#' is
#'
#' \deqn{
#'   B = g^{-1}(X_{ref}\beta).
#' }
#'
#' Let \eqn{X_{target}} replace only the focal exposure history by the evaluated
#' profile. Then
#'
#' \deqn{
#'   P = g^{-1}(X_{target}\beta)
#' }
#'
#' and
#'
#' \deqn{
#'   \Delta = P-B.
#' }
#'
#' With `scale = "response"`, `ECI_weighted = Delta`.
#'
#' ## Percent scale
#'
#' With a log link,
#'
#' \deqn{
#'   100\{\exp(\Delta\eta)-1\}
#' }
#'
#' is the percent relative change in the population expected response.
#'
#' With a logit link, the same mathematical transformation is the percent change
#' in odds. It is not the percent change or percentage-point change in predicted
#' probability.
#'
#' Other links do not receive this exponential interpretation, so
#' `scale = "percent"` is rejected explicitly.
#'
#' ## Population-level prediction contract
#'
#' Weighted ECI follows the same EpiExposure v1 prediction target as
#' `predict_outcomes()` and `summarise_effects()`:
#'
#' \deqn{
#'   E(Y\mid X,\theta)
#' }
#'
#' using only the fixed/population component. Fitted group-specific random
#' effects are set to zero/excluded. No residual, observation, process,
#' dispersion, or posterior-predictive noise is added.
#'
#' With `uncertainty = FALSE`, frequentist fits use fitted fixed coefficients
#' and Bayesian fits use posterior-mean fixed/population coefficients.
#'
#' ## Uncertainty
#'
#' With `uncertainty = TRUE`, one coherent fixed/population parameter draw is
#' applied jointly to the reference design and every target ECI design. For
#' draw \eqn{s},
#'
#' \deqn{
#'   \Delta\eta_i^{(s)} =
#'   (X_i-X_{ref})\beta^{(s)},
#' }
#'
#' \deqn{
#'   B^{(s)} = g^{-1}(X_{ref}\beta^{(s)}),
#' }
#'
#' \deqn{
#'   P_i^{(s)} = g^{-1}(X_i\beta^{(s)}),
#' }
#'
#' and
#'
#' \deqn{
#'   \Delta_i^{(s)} = P_i^{(s)}-B^{(s)}.
#' }
#'
#' The requested `ECI_weighted` transformation is then applied **within each
#' draw** before medians, SDs, and empirical quantiles are calculated. The same
#' parameter draw determines the baseline and every target row, preserving
#' covariance among variables, groups, baseline, target predictions, and the
#' weighted ECI.
#'
#' Frequentist engines use the joint asymptotic fixed-effect covariance matrix.
#' Bayesian engines use joint posterior or approximate-posterior fixed-effect
#' draws through the centralized EpiExposure prediction helpers. The function
#' does not independently reconstruct engine-specific uncertainty and does not
#' add residual or posterior-predictive outcome noise.
#'
#' ## Method-based reference values
#'
#' Method-based references are resolved once from the complete supplied
#' `data`, not independently within each group. This ensures every evaluated
#' group is compared with the same joint exposure reference condition.
#'
#' @export
compute_eci <- function(
    profile = NULL,
    data = NULL,
    group = NULL,
    time = "time",
    group_level = NULL,
    fit = NULL,
    var = NULL,
    ref = list(
      method = "median",
      value = NULL
    ),
    scale = c(
      "link",
      "response",
      "percent"
    ),
    uncertainty = FALSE,
    output = c(
      "summary",
      "samples"
    ),
    n_samples = 1000,
    interval_probs = c(
      0.025,
      0.975
    ),
    seed = NULL,
    extrapolation = c(
      "error",
      "warn",
      "allow"
    )
) {

  # ==========================================================================
  # SMALL HELPERS
  # ==========================================================================

  valid_name <- function(x) {
    is.character(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      nzchar(x)
  }

  valid_flag <- function(x) {
    is.logical(x) &&
      length(x) == 1L &&
      !is.na(x)
  }

  valid_integer_scalar <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x) &&
      x == as.integer(x)
  }

  validate_interval_probs <- function(x) {
    if (!is.numeric(x) ||
        length(x) != 2L ||
        anyNA(x) ||
        any(!is.finite(x)) ||
        any(x <= 0 | x >= 1) ||
        x[1L] >= x[2L]) {
      stop(
        "`interval_probs` must contain two finite probabilities satisfying ",
        "0 < interval_probs[1] < interval_probs[2] < 1.",
        call. = FALSE
      )
    }

    as.numeric(x)
  }

  summarise_numeric <- function(
    x,
    probs
  ) {
    if (!is.numeric(x) ||
        !length(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "Internal ECI uncertainty values must be finite.",
        call. = FALSE
      )
    }

    c(
      estimate =
        stats::median(x),
      sd =
        stats::sd(x),
      lower =
        stats::quantile(
          x,
          probs = probs[1L],
          names = FALSE,
          type = 7
        ),
      upper =
        stats::quantile(
          x,
          probs = probs[2L],
          names = FALSE,
          type = 7
        )
    )
  }

  attach_eci_attributes <- function(
    out,
    metadata = NULL,
    ref_values = NULL,
    output_kind = "summary"
  ) {
    attr(
      out,
      "epiexposure_eci_raw_contract"
    ) <- "discrete_sum_over_fitted_lag_window"

    attr(
      out,
      "epiexposure_eci_weighted_contract"
    ) <- if (is.null(metadata)) {
      "not_available_without_fit"
    } else {
      "centered_profile_vs_joint_reference"
    }

    attr(
      out,
      "epiexposure_eci_output"
    ) <- output_kind

    attr(
      out,
      "epiexposure_eci_history_contract"
    ) <- if (is.null(metadata)) {
      "raw_only_no_fitted_max_lag"
    } else {
      "exact_max_lag_plus_one"
    }

    if (!is.null(metadata)) {
      attr(
        out,
        "epiexposure_family_name"
      ) <- metadata$family_name

      attr(
        out,
        "epiexposure_link"
      ) <- metadata$link

      attr(
        out,
        "epiexposure_prediction_level"
      ) <- "population"

      attr(
        out,
        "epiexposure_prediction_estimand"
      ) <- "expected_response"

      attr(
        out,
        "epiexposure_point_prediction_contract"
      ) <- "central_expected_response"

      attr(
        out,
        "epiexposure_uncertainty_contract"
      ) <- "draw_by_draw_median_quantiles"

      attr(
        out,
        "epiexposure_eci_scale"
      ) <- scale

      attr(
        out,
        "epiexposure_eci_reference_values"
      ) <- ref_values

      attr(
        out,
        "epiexposure_uncertainty"
      ) <- uncertainty

      if (uncertainty) {
        attr(
          out,
          "epiexposure_n_samples"
        ) <- n_samples

        if (identical(
          output_kind,
          "summary"
        )) {
          attr(
            out,
            "epiexposure_interval_probs"
          ) <- interval_probs

          attr(
            out,
            "epiexposure_summary_center"
          ) <- "median"
        }
      }

      attr(
        out,
        "epiexposure_eci_percent_interpretation"
      ) <- if (!identical(
        scale,
        "percent"
      )) {
        "not_requested"
      } else if (identical(
        metadata$link,
        "log"
      )) {
        "percent_relative_change_in_expected_response"
      } else if (identical(
        metadata$link,
        "logit"
      )) {
        "percent_change_in_odds"
      } else {
        "not_applicable"
      }
    }

    out
  }

  # ==========================================================================
  # ARGUMENT MATCHING AND BASIC VALIDATION
  # ==========================================================================

  scale <- match.arg(scale)
  output <- match.arg(output)
  extrapolation <- match.arg(
    extrapolation
  )

  if (!valid_flag(uncertainty)) {
    stop(
      "`uncertainty` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!uncertainty &&
      identical(
        output,
        "samples"
      )) {
    stop(
      "`output = 'samples'` requires `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  if (!valid_integer_scalar(
    n_samples
  ) ||
  n_samples < 1) {
    stop(
      "`n_samples` must be one positive integer.",
      call. = FALSE
    )
  }
  n_samples <- as.integer(
    n_samples
  )

  if (uncertainty &&
      n_samples < 2L) {
    stop(
      "`n_samples` must be at least 2 when `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  interval_probs <-
    validate_interval_probs(
      interval_probs
    )

  if (!is.null(seed)) {
    if (!valid_integer_scalar(seed)) {
      stop(
        "`seed` must be NULL or one finite integer.",
        call. = FALSE
      )
    }
    seed <- as.integer(seed)
  }

  using_profile <- !is.null(
    profile
  )
  using_data <- !is.null(
    data
  )

  if (identical(
    using_profile,
    using_data
  )) {
    stop(
      "Supply exactly one of `profile` or `data`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # RAW-ONLY MODE WITHOUT A FIT
  # ==========================================================================

  if (is.null(fit)) {
    if (uncertainty) {
      stop(
        "`uncertainty = TRUE` requires a fitted model. With `fit = NULL`, ",
        "`compute_eci()` calculates descriptive raw cumulative exposure only.",
        call. = FALSE
      )
    }

    if (!using_profile) {
      stop(
        "When `fit = NULL`, supply one numeric `profile`.",
        call. = FALSE
      )
    }

    if (!is.numeric(profile) ||
        is.list(profile) ||
        !is.null(dim(profile)) ||
        !length(profile) ||
        anyNA(profile) ||
        any(!is.finite(profile))) {
      stop(
        "When `fit = NULL`, `profile` must be one finite numeric vector.",
        call. = FALSE
      )
    }

    out <- data.frame(
      ECI_raw =
        sum(
          as.numeric(profile)
        ),
      n_exposure_values =
        length(profile),
      stringsAsFactors = FALSE
    )

    out <- attach_eci_attributes(
      out,
      metadata = NULL,
      output_kind = "summary"
    )

    return(out)
  }

  # ==========================================================================
  # STRICT EPIEXPOSURE FIT CONTRACT
  # ==========================================================================

  required_helpers <- c(
    ".get_epiexposure_metadata",
    ".extract_central_parameters",
    ".extract_parameter_draws",
    ".epix_link_object",
    ".build_newdata_basis",
    ".epix_standard_fixed_design",
    ".epix_validate_regular_time"
  )

  # `exists()` without an explicit environment is unsafe inside `vapply()`
  # because its default search frame becomes the iterator's evaluation frame,
  # not necessarily the EpiExposure namespace. Search explicitly from the
  # current function evaluation environment, whose enclosing environment is the
  # package namespace when EpiExposure is loaded normally.
  helper_env <- environment()

  missing_helpers <- required_helpers[
    !vapply(
      required_helpers,
      function(helper) {
        exists(
          helper,
          envir = helper_env,
          mode = "function",
          inherits = TRUE
        )
      },
      logical(1)
    )
  ]

  if (length(missing_helpers)) {
    stop(
      "Current EpiExposure internal prediction helper(s) are unavailable: ",
      paste(
        missing_helpers,
        collapse = ", "
      ),
      ". Load the complete current EpiExposure package.",
      call. = FALSE
    )
  }

  metadata <- .get_epiexposure_metadata(
    fit
  )

  if (!identical(
    metadata$prediction_level,
    "population"
  )) {
    stop(
      "`compute_eci()` requires the population/fixed-component prediction ",
      "contract.",
      call. = FALSE
    )
  }

  if (!identical(
    metadata$prediction_estimand,
    "expected_response"
  )) {
    stop(
      "`compute_eci()` requires the expected-response prediction estimand.",
      call. = FALSE
    )
  }

  if (!identical(
    metadata$point_prediction_contract,
    "central_expected_response"
  ) ||
  !identical(
    metadata$uncertainty_contract,
    "draw_by_draw_median_quantiles"
  )) {
    stop(
      "The fitted model does not satisfy the current EpiExposure point/",
      "uncertainty prediction contract.",
      call. = FALSE
    )
  }

  link_object <- .epix_link_object(
    metadata$link
  )

  if (identical(
    scale,
    "percent"
  ) &&
  !metadata$link %in%
  c(
    "log",
    "logit"
  )) {
    stop(
      "`scale = 'percent'` is defined only for fitted log and logit links. ",
      "Use `scale = 'link'` or `scale = 'response'` for link '",
      metadata$link,
      "'.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # RESOLVE FOCAL VARIABLES
  # ==========================================================================

  if (is.null(var)) {
    if (length(
      metadata$vars
    ) == 1L) {
      variables <-
        metadata$vars
    } else {
      stop(
        "The fitted model contains multiple exposures: ",
        paste(
          metadata$vars,
          collapse = ", "
        ),
        ". Supply `var` explicitly.",
        call. = FALSE
      )
    }

  } else {
    if (!is.character(var) ||
        !length(var) ||
        anyNA(var) ||
        any(!nzchar(var)) ||
        anyDuplicated(var)) {
      stop(
        "`var` must be NULL or contain unique, non-empty exposure names.",
        call. = FALSE
      )
    }

    unknown <- setdiff(
      var,
      metadata$vars
    )

    if (length(unknown)) {
      stop(
        "Exposure variable(s) not found in the fitted model: ",
        paste(
          unknown,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    variables <- var
  }

  focal_history_lengths <- vapply(
    variables,
    function(variable) {
      max_lag_value <- metadata$spec[[variable]]$max_lag

      if (!is.numeric(max_lag_value) ||
          length(max_lag_value) != 1L ||
          is.na(max_lag_value) ||
          !is.finite(max_lag_value) ||
          max_lag_value < 0 ||
          max_lag_value != as.integer(max_lag_value)) {
        stop(
          "Invalid fitted `max_lag` metadata for exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      as.integer(max_lag_value) + 1L
    },
    integer(1)
  )

  # ==========================================================================
  # VALIDATE / NORMALIZE INPUT PROFILES OR DATA
  # ==========================================================================

  direct_profiles <- NULL
  levels_to_use <- NULL
  data_group_index <- NULL

  if (using_profile) {
    if (is.numeric(profile) &&
        is.null(dim(profile))) {
      if (length(
        variables
      ) != 1L) {
        stop(
          "A numeric `profile` can be used only when exactly one focal ",
          "exposure is requested.",
          call. = FALSE
        )
      }

      if (!length(profile) ||
          anyNA(profile) ||
          any(!is.finite(profile))) {
        stop(
          "`profile` must contain only finite numeric values.",
          call. = FALSE
        )
      }

      direct_profiles <-
        stats::setNames(
          list(
            as.numeric(profile)
          ),
          variables
        )

    } else if (is.list(profile)) {
      if (!length(profile)) {
        stop(
          "`profile` cannot be an empty list.",
          call. = FALSE
        )
      }

      if (is.null(
        names(profile)
      )) {
        if (length(profile) !=
            length(variables)) {
          stop(
            "An unnamed `profile` list must contain exactly one vector for ",
            "each requested variable, in the same order as `var`.",
            call. = FALSE
          )
        }

        names(profile) <-
          variables

      } else {
        if (anyNA(
          names(profile)
        ) ||
        any(!nzchar(
          names(profile)
        )) ||
        anyDuplicated(
          names(profile)
        )) {
          stop(
            "A named `profile` list must have unique, non-empty names.",
            call. = FALSE
          )
        }

        missing_profiles <- setdiff(
          variables,
          names(profile)
        )

        extra_profiles <- setdiff(
          names(profile),
          variables
        )

        if (length(missing_profiles) ||
            length(extra_profiles)) {
          details <- c(
            if (length(
              missing_profiles
            )) {
              paste0(
                "missing: ",
                paste(
                  missing_profiles,
                  collapse = ", "
                )
              )
            },
            if (length(
              extra_profiles
            )) {
              paste0(
                "not requested: ",
                paste(
                  extra_profiles,
                  collapse = ", "
                )
              )
            }
          )

          stop(
            "`profile` names must match `var` exactly (",
            paste(
              details,
              collapse = "; "
            ),
            ").",
            call. = FALSE
          )
        }

        profile <- profile[
          variables
        ]
      }

      for (variable in variables) {
        current <-
          profile[[variable]]

        if (!is.numeric(current) ||
            is.list(current) ||
            !is.null(dim(current)) ||
            !length(current) ||
            anyNA(current) ||
            any(!is.finite(current))) {
          stop(
            "Profile for exposure '",
            variable,
            "' must be one finite numeric vector.",
            call. = FALSE
          )
        }

        profile[[variable]] <-
          as.numeric(current)
      }

      direct_profiles <-
        profile

    } else {
      stop(
        "`profile` must be one numeric vector or a list of numeric vectors.",
        call. = FALSE
      )
    }

    for (variable in variables) {
      received_length <- length(
        direct_profiles[[variable]]
      )
      required_length <- focal_history_lengths[[variable]]

      if (received_length != required_length) {
        stop(
          "Exposure profile for '",
          variable,
          "' must contain exactly " ,
          required_length,
          " observations for fitted max_lag = " ,
          required_length - 1L,
          "; received " ,
          received_length,
          ". Longer histories are not truncated and shorter histories are " ,
          "not padded.",
          call. = FALSE
        )
      }
    }

  } else {
    if (!is.data.frame(data) ||
        !nrow(data)) {
      stop(
        "`data` must be a non-empty data.frame.",
        call. = FALSE
      )
    }

    if (!valid_name(group) ||
        !group %in%
        names(data)) {
      stop(
        "When `data` is supplied, `group` must name a grouping column present ",
        "in `data`.",
        call. = FALSE
      )
    }

    if (!valid_name(time) ||
        !time %in%
        names(data)) {
      stop(
        "`time` must name a chronological time column present in `data`.",
        call. = FALSE
      )
    }

    if (anyNA(
      data[[group]]
    )) {
      stop(
        "The grouping column cannot contain missing values.",
        call. = FALSE
      )
    }

    if (!is.numeric(
      data[[time]]
    ) ||
    anyNA(
      data[[time]]
    ) ||
    any(!is.finite(
      data[[time]]
    ))) {
      stop(
        "`time` must contain only finite numeric values.",
        call. = FALSE
      )
    }

    missing_focal <- setdiff(
      variables,
      names(data)
    )

    if (length(missing_focal)) {
      stop(
        "`data` is missing requested exposure variable(s): ",
        paste(
          missing_focal,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    for (variable in variables) {
      if (!is.numeric(
        data[[variable]]
      ) ||
      anyNA(
        data[[variable]]
      ) ||
      any(!is.finite(
        data[[variable]]
      ))) {
        stop(
          "Exposure variable '",
          variable,
          "' must contain only finite numeric values.",
          call. = FALSE
        )
      }
    }

    available_levels <-
      unique(
        data[[group]]
      )

    if (is.null(
      group_level
    )) {
      levels_to_use <-
        available_levels

    } else {
      requested_character <-
        as.character(
          group_level
        )
      available_character <-
        as.character(
          available_levels
        )

      missing_levels <-
        requested_character[
          !requested_character %in%
            available_character
        ]

      if (length(
        missing_levels
      )) {
        stop(
          "Requested group level(s) not found in `data`: ",
          paste(
            unique(
              missing_levels
            ),
            collapse = ", "
          ),
          ".",
          call. = FALSE
        )
      }

      index <- match(
        requested_character,
        available_character
      )

      if (anyDuplicated(index)) {
        stop(
          "`group_level` must not contain duplicated group levels.",
          call. = FALSE
        )
      }

      levels_to_use <-
        available_levels[
          index
        ]
    }

    if (!length(
      levels_to_use
    )) {
      stop(
        "No group levels are available for ECI calculation.",
        call. = FALSE
      )
    }

    keep_groups <-
      as.character(
        data[[group]]
      ) %in%
      as.character(
        levels_to_use
      )

    validation_data <-
      data[
        keep_groups,
        ,
        drop = FALSE
      ]

    .epix_validate_regular_time(
      validation_data,
      group_index =
        as.character(
          validation_data[[
            group
          ]]
        ),
      time_col =
        time
    )

    unique_required_lengths <- unique(
      unname(
        focal_history_lengths
      )
    )

    if (length(unique_required_lengths) != 1L) {
      stop(
        "The requested focal exposures have different fitted `max_lag` values " ,
        "and therefore cannot share one exact-length long-format history in " ,
        "`data`. Requested lengths: " ,
        paste(
          paste0(
            names(focal_history_lengths),
            "=",
            focal_history_lengths
          ),
          collapse = ", "
        ),
        ". Evaluate exposures with different fitted lag windows separately.",
        call. = FALSE
      )
    }

    required_group_length <- unique_required_lengths[[1L]]
    group_counts <- table(
      as.character(
        validation_data[[group]]
      )
    )

    invalid_group_counts <- group_counts[
      group_counts != required_group_length
    ]

    if (length(invalid_group_counts)) {
      examples <- paste0(
        names(invalid_group_counts),
        "=",
        as.integer(invalid_group_counts)
      )

      stop(
        "Every evaluated group must contain exactly " ,
        required_group_length,
        " observations for fitted max_lag = " ,
        required_group_length - 1L,
        ". Group(s) with non-matching history length include: " ,
        paste(
          utils::head(
            examples,
            5L
          ),
          collapse = ", "
        ),
        if (length(examples) > 5L) "; ..." else ".",
        " Longer histories are not truncated and shorter histories are not padded.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # RESOLVE JOINT REFERENCE VALUES
  # ==========================================================================

  ref_is_method <-
    is.list(ref) &&
    identical(
      sort(
        names(ref)
      ),
      sort(
        c(
          "method",
          "value"
        )
      )
    )

  if (ref_is_method) {
    method <-
      ref$method

    if (!valid_name(method) ||
        !method %in%
        c(
          "median",
          "percentile",
          "fixed"
        )) {
      stop(
        "Method-based `ref$method` must be one of 'median', 'percentile', or ",
        "'fixed'.",
        call. = FALSE
      )
    }

    if (using_profile) {
      stop(
        "Method-based `ref` requires `data` so reference exposure values can ",
        "be derived from an observed exposure distribution. With direct ",
        "`profile` input, supply a named reference value for every fitted ",
        "exposure.",
        call. = FALSE
      )
    }

    missing_reference_data <- setdiff(
      metadata$vars,
      names(data)
    )

    if (length(
      missing_reference_data
    )) {
      stop(
        "Method-based `ref` requires every fitted exposure in `data`. Missing: ",
        paste(
          missing_reference_data,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    for (variable in metadata$vars) {
      if (!is.numeric(
        data[[variable]]
      ) ||
      anyNA(
        data[[variable]]
      ) ||
      any(!is.finite(
        data[[variable]]
      ))) {
        stop(
          "Reference exposure variable '",
          variable,
          "' must contain only finite numeric values.",
          call. = FALSE
        )
      }
    }

    if (identical(
      method,
      "median"
    )) {
      if (!is.null(
        ref$value
      )) {
        stop(
          "`ref$value` must be NULL when `ref$method = 'median'`.",
          call. = FALSE
        )
      }

      ref_values <- vapply(
        metadata$vars,
        function(variable) {
          stats::median(
            data[[variable]]
          )
        },
        numeric(1)
      )

    } else if (identical(
      method,
      "percentile"
    )) {
      p <- ref$value

      if (!is.numeric(p) ||
          length(p) != 1L ||
          is.na(p) ||
          !is.finite(p) ||
          p < 0 ||
          p > 1) {
        stop(
          "`ref$value` must be one finite probability between 0 and 1 when ",
          "`ref$method = 'percentile'`.",
          call. = FALSE
        )
      }

      ref_values <- vapply(
        metadata$vars,
        function(variable) {
          stats::quantile(
            data[[variable]],
            probs = p,
            names = FALSE,
            type = 7
          )
        },
        numeric(1)
      )

    } else {
      fixed_value <-
        ref$value

      if (!is.numeric(
        fixed_value
      ) ||
      length(
        fixed_value
      ) != 1L ||
      is.na(
        fixed_value
      ) ||
      !is.finite(
        fixed_value
      )) {
        stop(
          "`ref$value` must be one finite numeric value when ",
          "`ref$method = 'fixed'`.",
          call. = FALSE
        )
      }

      ref_values <-
        stats::setNames(
          rep(
            as.numeric(
              fixed_value
            ),
            length(
              metadata$vars
            )
          ),
          metadata$vars
        )
    }

  } else {
    if (!is.list(ref) ||
        is.null(
          names(ref)
        ) ||
        anyNA(
          names(ref)
        ) ||
        any(!nzchar(
          names(ref)
        )) ||
        anyDuplicated(
          names(ref)
        )) {
      stop(
        "`ref` must be either a method-based list with `method` and `value`, ",
        "or a named list with one finite reference value for every fitted ",
        "exposure.",
        call. = FALSE
      )
    }

    missing_ref <- setdiff(
      metadata$vars,
      names(ref)
    )

    extra_ref <- setdiff(
      names(ref),
      metadata$vars
    )

    if (length(missing_ref) ||
        length(extra_ref)) {
      details <- c(
        if (length(
          missing_ref
        )) {
          paste0(
            "missing: ",
            paste(
              missing_ref,
              collapse = ", "
            )
          )
        },
        if (length(
          extra_ref
        )) {
          paste0(
            "not fitted: ",
            paste(
              extra_ref,
              collapse = ", "
            )
          )
        }
      )

      stop(
        "Exposure-specific `ref` names must match all fitted exposures exactly ",
        "(",
        paste(
          details,
          collapse = "; "
        ),
        ").",
        call. = FALSE
      )
    }

    ref <- ref[
      metadata$vars
    ]

    ref_values <- vapply(
      metadata$vars,
      function(variable) {
        value <-
          ref[[variable]]

        if (!is.numeric(value) ||
            length(value) != 1L ||
            is.na(value) ||
            !is.finite(value)) {
          stop(
            "Reference value for exposure '",
            variable,
            "' must be one finite numeric value.",
            call. = FALSE
          )
        }

        as.numeric(value)
      },
      numeric(1)
    )
  }

  names(ref_values) <-
    metadata$vars

  # ==========================================================================
  # CONSTRUCT THE JOINT REFERENCE DESIGN
  # ==========================================================================

  reference_profiles <-
    stats::setNames(
      vector(
        "list",
        length(
          metadata$vars
        )
      ),
      metadata$vars
    )

  for (variable in metadata$vars) {
    max_lag_use <-
      metadata$spec[[
        variable
      ]]$max_lag

    reference_profiles[[
      variable
    ]] <- rep(
      ref_values[[
        variable
      ]],
      max_lag_use + 1L
    )
  }

  reference_design <-
    .build_newdata_basis(
      fit = fit,
      profiles =
        reference_profiles,
      extrapolation =
        extrapolation
    )

  X_reference <-
    .epix_standard_fixed_design(
      reference_design,
      metadata
    )

  if (nrow(
    X_reference
  ) != 1L) {
    stop(
      "Internal ECI reference design must contain exactly one row.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # BUILD FOCAL TARGET DESIGNS
  # ==========================================================================

  target_records <- list()
  target_designs <- list()
  target_index <- 1L

  add_target <- function(
    variable,
    profile_values,
    group_value = NULL
  ) {
    max_lag_use <-
      metadata$spec[[
        variable
      ]]$max_lag

    required_length <-
      max_lag_use + 1L

    if (!is.numeric(
      profile_values
    ) ||
    !length(
      profile_values
    ) ||
    anyNA(
      profile_values
    ) ||
    any(!is.finite(
      profile_values
    ))) {
      stop(
        "Profile for exposure '",
        variable,
        "' must contain only finite numeric values.",
        call. = FALSE
      )
    }

    if (length(
      profile_values
    ) !=
    required_length) {
      stop(
        "Exposure profile for '",
        variable,
        "' must contain exactly ",
        required_length,
        " observations for fitted max_lag = ",
        max_lag_use,
        "; received ",
        length(
          profile_values
        ),
        ". Longer histories are not truncated and shorter histories are not padded.",
        call. = FALSE
      )
    }

    profile_used <-
      as.numeric(
        profile_values
      )

    target_profiles <-
      reference_profiles

    target_profiles[[
      variable
    ]] <-
      profile_used

    target_design <-
      .build_newdata_basis(
        fit = fit,
        profiles =
          target_profiles,
        extrapolation =
          extrapolation
      )

    X_target <-
      .epix_standard_fixed_design(
        target_design,
        metadata
      )

    if (nrow(
      X_target
    ) != 1L ||
    !identical(
      colnames(
        X_target
      ),
      colnames(
        X_reference
      )
    )) {
      stop(
        "Internal ECI target design is inconsistent with the reference design.",
        call. = FALSE
      )
    }

    record <- data.frame(
      var =
        variable,
      reference_value =
        ref_values[[
          variable
        ]],
      max_lag =
        max_lag_use,
      n_exposure_values =
        required_length,
      ECI_raw =
        sum(
          profile_used
        ),
      ECI_raw_centered =
        sum(
          profile_used -
            ref_values[[
              variable
            ]]
        ),
      stringsAsFactors = FALSE
    )

    if (!is.null(
      group_value
    )) {
      record[[
        group
      ]] <-
        group_value

      record <- record[
        ,
        c(
          group,
          setdiff(
            names(record),
            group
          )
        ),
        drop = FALSE
      ]
    }

    target_records[[
      target_index
    ]] <<-
      record

    target_designs[[
      target_index
    ]] <<-
      X_target

    target_index <<-
      target_index + 1L

    invisible(NULL)
  }

  if (using_profile) {
    for (variable in variables) {
      add_target(
        variable =
          variable,
        profile_values =
          direct_profiles[[
            variable
          ]]
      )
    }

  } else {
    for (level_index in seq_along(
      levels_to_use
    )) {
      current_level <-
        levels_to_use[[
          level_index
        ]]

      group_rows <-
        as.character(
          data[[group]]
        ) ==
        as.character(
          current_level
        )

      group_data <-
        data[
          group_rows,
          ,
          drop = FALSE
        ]

      if (!nrow(
        group_data
      )) {
        stop(
          "No observations found for group level '",
          as.character(
            current_level
          ),
          "'.",
          call. = FALSE
        )
      }

      order_index <-
        order(
          group_data[[
            time
          ]]
        )

      group_data <-
        group_data[
          order_index,
          ,
          drop = FALSE
        ]

      for (variable in variables) {
        add_target(
          variable =
            variable,
          profile_values =
            group_data[[
              variable
            ]],
          group_value =
            current_level
        )
      }
    }
  }

  target_metadata <-
    do.call(
      rbind,
      target_records
    )
  rownames(
    target_metadata
  ) <- NULL

  X_targets <-
    do.call(
      rbind,
      target_designs
    )

  rownames(
    X_targets
  ) <- NULL

  if (nrow(
    X_targets
  ) !=
  nrow(
    target_metadata
  )) {
    stop(
      "Internal ECI target metadata/design row counts do not match.",
      call. = FALSE
    )
  }

  # Verify that each target differs from the reference only in the focal
  # exposure's fitted cross-basis columns.
  for (i in seq_len(
    nrow(
      X_targets
    )
  )) {
    variable <-
      target_metadata$var[[
        i
      ]]

    focal_columns <-
      paste0(
        "cb_",
        variable,
        "_"
      )

    changed <-
      colnames(
        X_targets
      )[
        abs(
          X_targets[
            i,
            ,
            drop = TRUE
          ] -
            X_reference[
              1L,
              ,
              drop = TRUE
            ]
        ) >
          sqrt(
            .Machine$double.eps
          )
      ]

    nonfocal_changed <-
      changed[
        changed !=
          "(Intercept)" &
          !startsWith(
            changed,
            focal_columns
          )
      ]

    if (length(
      nonfocal_changed
    )) {
      stop(
        "Internal ECI construction changed non-focal fitted terms for exposure '",
        variable,
        "': ",
        paste(
          nonfocal_changed,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # DETERMINISTIC CENTRAL-PARAMETER ECI
  # ==========================================================================

  transform_weighted <- function(
    eta,
    delta
  ) {
    if (identical(
      scale,
      "link"
    )) {
      return(
        eta
      )
    }

    if (identical(
      scale,
      "response"
    )) {
      return(
        delta
      )
    }

    100 *
      (
        exp(
          eta
        ) -
          1
      )
  }

  if (!uncertainty) {
    beta <-
      .extract_central_parameters(
        fit
      )

    if (!identical(
      names(beta),
      colnames(
        X_reference
      )
    )) {
      stop(
        "Central parameter names do not align exactly with the ECI design.",
        call. = FALSE
      )
    }

    baseline_eta <-
      as.numeric(
        X_reference %*%
          beta
      )

    target_eta <-
      as.numeric(
        X_targets %*%
          beta
      )

    eta_contrast <-
      target_eta -
      baseline_eta

    baseline_response <-
      as.numeric(
        link_object$linkinv(
          baseline_eta
        )
      )

    predicted_response <-
      as.numeric(
        link_object$linkinv(
          target_eta
        )
      )

    if (any(!is.finite(
      c(
        baseline_response,
        predicted_response
      )
    ))) {
      stop(
        "Inverse-link transformation produced non-finite ECI response values.",
        call. = FALSE
      )
    }

    delta_response <-
      predicted_response -
      baseline_response

    weighted <-
      transform_weighted(
        eta =
          eta_contrast,
        delta =
          delta_response
      )

    if (any(!is.finite(
      weighted
    ))) {
      stop(
        "The requested weighted ECI transformation produced non-finite values.",
        call. = FALSE
      )
    }

    out <-
      target_metadata

    out$baseline <-
      rep(
        baseline_response,
        nrow(
          out
        )
      )
    out$predicted <-
      predicted_response
    out$ECI_weighted <-
      weighted
    out$scale <-
      scale

    out <-
      attach_eci_attributes(
        out,
        metadata =
          metadata,
        ref_values =
          ref_values,
        output_kind =
          "summary"
      )

    return(out)
  }

  # ==========================================================================
  # PARAMETER-DRAW ECI
  # ==========================================================================

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
      if (!is.null(seed)) {
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
      }
    },
    add = TRUE
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  parameter_draws <-
    .extract_parameter_draws(
      fit,
      n_samples =
        n_samples
    )

  if (!is.matrix(
    parameter_draws
  ) ||
  nrow(
    parameter_draws
  ) !=
  n_samples ||
  !identical(
    colnames(
      parameter_draws
    ),
    colnames(
      X_reference
    )
  ) ||
  any(!is.finite(
    parameter_draws
  ))) {
    stop(
      "Parameter draws do not satisfy the current EpiExposure ECI design ",
      "contract.",
      call. = FALSE
    )
  }

  baseline_eta_draws <-
    as.numeric(
      parameter_draws %*%
        as.numeric(
          X_reference[
            1L,
            ,
            drop = TRUE
          ]
        )
    )

  target_eta_draws <-
    parameter_draws %*%
    t(
      X_targets
    )

  if (!is.matrix(
    target_eta_draws
  )) {
    target_eta_draws <-
      matrix(
        target_eta_draws,
        nrow =
          n_samples
      )
  }

  eta_contrast_draws <-
    sweep(
      target_eta_draws,
      MARGIN = 1L,
      STATS =
        baseline_eta_draws,
      FUN = "-"
    )

  baseline_response_draws <-
    as.numeric(
      link_object$linkinv(
        baseline_eta_draws
      )
    )

  predicted_response_draws <-
    matrix(
      link_object$linkinv(
        as.vector(
          target_eta_draws
        )
      ),
      nrow =
        n_samples,
      ncol =
        nrow(
          X_targets
        )
    )

  if (any(!is.finite(
    baseline_response_draws
  )) ||
  any(!is.finite(
    predicted_response_draws
  ))) {
    stop(
      "Inverse-link transformation produced non-finite ECI parameter-draw ",
      "predictions.",
      call. = FALSE
    )
  }

  delta_response_draws <-
    sweep(
      predicted_response_draws,
      MARGIN = 1L,
      STATS =
        baseline_response_draws,
      FUN = "-"
    )

  weighted_draws <- if (
    identical(
      scale,
      "link"
    )
  ) {
    eta_contrast_draws

  } else if (identical(
    scale,
    "response"
  )) {
    delta_response_draws

  } else {
    100 *
      (
        exp(
          eta_contrast_draws
        ) -
          1
      )
  }

  if (any(!is.finite(
    weighted_draws
  ))) {
    stop(
      "The requested weighted ECI transformation produced non-finite parameter ",
      "draw values.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # SAMPLE-LEVEL OUTPUT
  # ==========================================================================

  if (identical(
    output,
    "samples"
  )) {
    sample_blocks <- vector(
      "list",
      nrow(
        target_metadata
      )
    )

    for (i in seq_len(
      nrow(
        target_metadata
      )
    )) {
      current <-
        target_metadata[
          rep(
            i,
            n_samples
          ),
          ,
          drop = FALSE
        ]

      current$sample <-
        seq_len(
          n_samples
        )


      current$baseline <-
        baseline_response_draws

      current$predicted <-
        predicted_response_draws[
          ,
          i
        ]


      current$ECI_weighted <-
        weighted_draws[
          ,
          i
        ]

      current$scale <-
        scale

      # Place sample after identifiers/descriptive columns and before model
      # quantities for predictable downstream use.
      descriptive <- c(
        if (!is.null(group) &&
            group %in%
            names(current)) {
          group
        },
        "var",
        "reference_value",
        "max_lag",
        "n_exposure_values",
        "ECI_raw",
        "ECI_raw_centered",
        "sample",
        "baseline",
        "predicted",
        "ECI_weighted",
        "scale"
      )

      current <-
        current[
          ,
          descriptive,
          drop = FALSE
        ]

      sample_blocks[[
        i
      ]] <-
        current
    }

    out <-
      do.call(
        rbind,
        sample_blocks
      )

    rownames(out) <-
      NULL

    out <-
      attach_eci_attributes(
        out,
        metadata =
          metadata,
        ref_values =
          ref_values,
        output_kind =
          "samples"
      )

    return(out)
  }

  # ==========================================================================
  # UNCERTAINTY SUMMARY OUTPUT
  # ==========================================================================

  summary_rows <- vector(
    "list",
    nrow(
      target_metadata
    )
  )

  for (i in seq_len(
    nrow(
      target_metadata
    )
  )) {
    current <-
      target_metadata[
        i,
        ,
        drop = FALSE
      ]

    quantities <- list(
      baseline =
        baseline_response_draws,
      predicted =
        predicted_response_draws[
          ,
          i
        ],
      ECI_weighted =
        weighted_draws[
          ,
          i
        ]
    )

    for (quantity_name in names(
      quantities
    )) {
      summary_values <-
        summarise_numeric(
          quantities[[
            quantity_name
          ]],
          interval_probs
        )

      current[[
        quantity_name
      ]] <-
        unname(
          summary_values[
            "estimate"
          ]
        )

      current[[
        paste0(
          quantity_name,
          "_sd"
        )
      ]] <-
        unname(
          summary_values[
            "sd"
          ]
        )

      current[[
        paste0(
          quantity_name,
          "_lower"
        )
      ]] <-
        unname(
          summary_values[
            "lower"
          ]
        )

      current[[
        paste0(
          quantity_name,
          "_upper"
        )
      ]] <-
        unname(
          summary_values[
            "upper"
          ]
        )
    }

    current$scale <-
      scale

    summary_rows[[
      i
    ]] <-
      current
  }

  out <-
    do.call(
      rbind,
      summary_rows
    )

  rownames(out) <-
    NULL

  out <-
    attach_eci_attributes(
      out,
      metadata =
        metadata,
      ref_values =
        ref_values,
      output_kind =
        "summary"
    )

  out
}
