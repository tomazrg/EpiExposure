#' Compute lag-specific decomposition of Exposure Cumulative Impact (ECI)
#'
#' Decomposes the centered model-weighted Exposure Cumulative Impact into exact
#' lag-specific contributions on the linear-predictor scale.
#'
#' `compute_ecilag()` is designed to complement `compute_eci()` and
#' `summarise_effects()`:
#'
#' - `summarise_effects(scale = "lag")` evaluates a general exposure-by-lag
#'   effect surface;
#' - `compute_eci()` evaluates the total impact of one complete exposure
#'   trajectory relative to a joint reference profile;
#' - `compute_ecilag()` decomposes that trajectory-specific total impact into
#'   exact lag-specific contributions.
#'
#' Exposure profiles are supplied chronologically from the oldest observation
#' to the most recent observation. The returned decomposition is indexed on the
#' retrospective lag scale, where lag 0 is the most recent exposure and
#' `max_lag` is the oldest exposure.
#'
#' @param profile Optional exposure-profile input. Supply exactly one of
#'   `profile` or `data`.
#'
#'   When exactly one exposure is requested, `profile` may be one finite numeric
#'   vector. For multiple requested exposures, supply a list with exactly one
#'   finite numeric vector per requested variable. A named list is recommended;
#'   an unnamed list must follow the same order as `var`.
#'
#'   Every supplied history must contain exactly `max_lag + 1` observations.
#'   Longer histories are not truncated and shorter histories are not padded.
#'
#'   Under the EpiExposure exact-history contract, all fitted exposure variables
#'   in the model must share the same fitted `max_lag`, and therefore the same
#'   history length.
#' @param data Optional non-empty long-format data frame containing exposure
#'   histories. Supply exactly one of `profile` or `data`.
#'
#'   The data must contain the grouping column named by `group`, the time column
#'   named by `time`, and every exposure requested in `var`. If a method-based
#'   reference is used, `data` must additionally contain every fitted exposure
#'   because those variables are needed to derive the joint reference values.
#'
#'   Every group in the supplied `data` must contain exactly `max_lag + 1`
#'   unique, equally spaced time points. All groups must use the same temporal
#'   spacing.
#' @param group Character scalar naming the independent-history grouping column
#'   in `data`. Required when `data` is supplied.
#' @param time Character scalar naming the chronological numeric time column in
#'   `data`. Default is `"time"`.
#' @param group_level Optional grouping value or vector of values to return. If
#'   `NULL`, all groups are evaluated in first-occurrence order. The exact
#'   temporal-history contract is validated for every group supplied in `data`,
#'   including groups not selected by `group_level`.
#' @param fit Fitted model returned by the current `fit_epidlnm()`.
#' @param var Optional character vector naming fitted exposure variables to
#'   decompose.
#'
#'   If the fitted model contains one exposure and `var = NULL`, that exposure
#'   is used. For a multivariable fit, `var` must be supplied explicitly.
#'
#'   Each requested exposure is decomposed separately. Other fitted exposures
#'   remain fixed at their reference values.
#' @param ref Reference exposure specification.
#'
#'   Two forms are supported:
#'
#'   \itemize{
#'     \item Method-based:
#'       `list(method = "median", value = NULL)`,
#'       `list(method = "percentile", value = p)`, or
#'       `list(method = "fixed", value = x)`.
#'     \item Exposure-specific: a named list containing exactly one finite
#'       numeric reference value for every fitted exposure, for example
#'       `list(tmean = 25, rain = 5, wetness = 10)`.
#'   }
#'
#'   Method-based references require `data`. `"median"` uses each fitted
#'   exposure's own median. `"percentile"` uses the same percentile `p` for
#'   every fitted exposure. `"fixed"` applies the same numeric value to all
#'   fitted exposures and should be used cautiously when exposure variables
#'   have different units.
#'
#'   With direct `profile` input, use the exposure-specific named-list form.
#'
#'   The default is `list(method = "median", value = NULL)`.
#' @param absolute Logical. If `TRUE`, add `abs_contribution` columns
#'   (and their uncertainty summaries).
#' @param uncertainty Logical. If `FALSE`, use the harmonized central
#'   fixed/population parameter estimate. If `TRUE`, propagate joint
#'   fixed/population parameter uncertainty draw by draw.
#' @param output Character. `"summary"` returns deterministic values when
#'   `uncertainty = FALSE`, or medians, SDs, and empirical intervals when
#'   `uncertainty = TRUE`. `"samples"` returns draw-level decomposition and
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
#'   `"allow"`.
#'
#' @return A list.
#'
#'   For one profile-variable combination, the list contains:
#'
#'   \describe{
#'     \item{`ECI_raw`}{Discrete sum of the observed focal exposure values.}
#'     \item{`ECI_raw_centered`}{Discrete sum of exposure deviations from the
#'       focal reference value.}
#'     \item{`ECI_weighted`}{Centered cumulative trajectory effect on the
#'       linear-predictor scale.}
#'     \item{`reference_value`}{Focal exposure reference value.}
#'     \item{`max_lag`}{Common fitted maximum lag.}
#'     \item{`n_exposure_values`}{Exact history length, `max_lag + 1`.}
#'     \item{`by_lag`}{Lag-specific decomposition table.}
#'   }
#'
#'   With uncertainty, `ECI_weighted_sd`, `ECI_weighted_lower`, and
#'   `ECI_weighted_upper` are added. With `output = "samples"`,
#'   `by_lag_samples` and `ECI_weighted_samples` are also returned.
#'
#'   For multiple groups and/or variables, the returned list contains
#'   `eci_summary` and `by_lag`, plus sample-level tables when requested.
#'
#'   `by_lag` contains:
#'
#'   \describe{
#'     \item{`var`}{Focal exposure.}
#'     \item{`lag`}{Retrospective lag; 0 is most recent.}
#'     \item{`exposure`}{Observed/simulated focal exposure at that lag.}
#'     \item{`reference_value`}{Focal reference exposure.}
#'     \item{`exposure_minus_reference`}{Exposure deviation from reference.}
#'     \item{`contribution`}{Exact centered lag-specific DLNM contribution on
#'       the linear-predictor scale.}
#'     \item{`percent_contribution`}{Share of total absolute contribution
#'       magnitude, `100 * abs(contribution) / sum(abs(contribution))`.}
#'   }
#'
#' @details
#' ## Exact history contract
#'
#' EpiExposure uses one complete discrete exposure history per epidemiological
#' unit. For a common fitted maximum lag `L`, every fitted exposure and every
#' evaluated group must therefore use exactly
#'
#' \deqn{
#'   L + 1
#' }
#'
#' observations.
#'
#' Histories with different lengths or fitted exposures with different
#' `max_lag` values are rejected explicitly. `compute_ecilag()` never selects a
#' trailing window, truncates older observations, pads shorter histories, or
#' silently aligns variables with different temporal windows.
#'
#' ## Centered total ECI
#'
#' For focal exposure trajectory
#'
#' \deqn{
#'   X = (x_L, \ldots, x_1, x_0)
#' }
#'
#' supplied chronologically from oldest to most recent, and a constant focal
#' reference trajectory
#'
#' \deqn{
#'   X_{ref} = (r, \ldots, r),
#' }
#'
#' the weighted ECI on the linear-predictor scale is
#'
#' \deqn{
#'   ECI_{weighted}
#'   =
#'   \{CB(X)-CB(X_{ref})\}^T\beta_{focal}.
#' }
#'
#' All non-focal fitted exposures are held at their own reference values.
#'
#' This is the same trajectory-level centered effect used by
#' `compute_eci(scale = "link")`.
#'
#' ## Exact lag-specific contribution
#'
#' For lag `l`, define a trajectory that equals the reference at every lag
#' except lag `l`, where the focal exposure takes its observed value `x_l`.
#' The exact contribution is
#'
#' \deqn{
#'   C_l
#'   =
#'   \eta_l(x_l)-\eta_l(r).
#' }
#'
#' Because a DLNM cross-basis is additive over lag-specific basis
#' contributions on the linear-predictor scale,
#'
#' \deqn{
#'   ECI_{weighted}
#'   =
#'   \sum_{l=0}^{L} C_l.
#' }
#'
#' `compute_ecilag()` verifies this identity numerically at the design-matrix
#' level and again for the central estimate. With uncertainty it additionally
#' verifies the equality draw by draw.
#'
#' For the same exposure value, lag, reference, basis, and coefficients,
#' `C_l` is numerically the same centered lag-specific linear effect represented
#' by `summarise_effects(scale = "lag", effect_measure = "linear")`.
#' The functions differ in purpose: `summarise_effects()` maps a general
#' exposure-by-lag surface, whereas `compute_ecilag()` extracts the values
#' actually occurring along one specific trajectory.
#'
#' ## Percent contribution
#'
#' The magnitude share is
#'
#' \deqn{
#'   P_l
#'   =
#'   100
#'   \frac{|C_l|}
#'   {\sum_j |C_j|}.
#' }
#'
#' The signed effect remains available in `contribution`. If every lag
#' contribution is exactly zero, percent contributions are undefined and
#' returned as `NA`.
#'
#' ## Why decomposition stays on the link scale
#'
#' Additivity is guaranteed on the linear-predictor scale. In general,
#'
#' \deqn{
#'   g^{-1}\left(\eta_{ref}+\sum_l C_l\right)
#'   -
#'   g^{-1}(\eta_{ref})
#' }
#'
#' cannot be decomposed into the sum of separately inverse-linked lag effects
#' when the link is nonlinear.
#'
#' Consequently, `compute_ecilag()` intentionally returns contribution and
#' `ECI_weighted` on the link scale only. Use `compute_eci(scale = "response")`
#' when the scientific target is the total trajectory impact on the expected
#' response scale.
#'
#' ## Population-level model contract
#'
#' The decomposition uses the same EpiExposure v1 fixed/population parameter
#' contract as `compute_eci()` and `predict_outcomes()`. Fitted group-specific
#' random effects are excluded/set to zero. No residual, observation, process,
#' dispersion, or posterior-predictive noise is added.
#'
#' With `uncertainty = FALSE`, frequentist fits use fitted fixed coefficients
#' and Bayesian fits use posterior-mean fixed/population coefficients.
#'
#' ## Uncertainty
#'
#' With `uncertainty = TRUE`, one joint fixed/population parameter draw
#' `beta^(s)` is applied simultaneously to every group, variable, and lag.
#' Thus, for each draw,
#'
#' \deqn{
#'   ECI_{weighted}^{(s)}
#'   =
#'   \sum_l C_l^{(s)}.
#' }
#'
#' Contributions, absolute quantities, and percent contributions are calculated
#' draw by draw before summaries are produced.
#' Summary output uses the median, empirical SD, and empirical quantiles.
#'
#' Because the median is nonlinear,
#'
#' \deqn{
#'   median\left(\sum_l C_l^{(s)}\right)
#' }
#'
#' does not generally equal
#'
#' \deqn{
#'   \sum_l median(C_l^{(s)}).
#' }
#'
#' Therefore, exact additivity is a deterministic and draw-level identity; it is
#' not expected between separately summarized medians.
#'
#' Frequentist and Bayesian parameter draws are obtained exclusively through
#' the centralized EpiExposure prediction helpers. This function does not
#' reconstruct engine-specific coefficient covariance or posterior logic.
#'
#' @export
compute_ecilag <- function(
    profile = NULL,
    data = NULL,
    group = NULL,
    time = "time",
    group_level = NULL,
    fit,
    var = NULL,
    ref = list(
      method = "median",
      value = NULL
    ),
    absolute = TRUE,
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

  summarise_finite <- function(
    x,
    probs
  ) {
    if (!is.numeric(x) ||
        !length(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "Internal ECI-lag draw values expected to be finite were not finite.",
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

  summarise_defined <- function(
    x,
    probs
  ) {
    if (!is.numeric(x) ||
        !length(x)) {
      stop(
        "Internal percent-contribution draw values must be numeric.",
        call. = FALSE
      )
    }

    finite <- is.finite(x)
    n_defined <- sum(finite)

    if (!n_defined) {
      return(
        c(
          estimate = NA_real_,
          sd = NA_real_,
          lower = NA_real_,
          upper = NA_real_,
          n_defined = 0
        )
      )
    }

    values <- x[finite]

    c(
      estimate =
        stats::median(values),
      sd =
        if (length(values) >= 2L) {
          stats::sd(values)
        } else {
          NA_real_
        },
      lower =
        stats::quantile(
          values,
          probs = probs[1L],
          names = FALSE,
          type = 7
        ),
      upper =
        stats::quantile(
          values,
          probs = probs[2L],
          names = FALSE,
          type = 7
        ),
      n_defined =
        n_defined
    )
  }

  nearly_equal <- function(
    observed,
    expected,
    tolerance =
      1000 *
      .Machine$double.eps
  ) {
    abs(
      observed -
        expected
    ) <=
      tolerance *
      pmax(
        1,
        abs(observed),
        abs(expected)
      )
  }

  attach_contract_attributes <- function(
    object,
    metadata,
    ref_values,
    common_max_lag,
    common_history_length,
    output_kind
  ) {
    attr(
      object,
      "epiexposure_ecilag_contract"
    ) <-
      "exact_lag_contributions_sum_to_centered_eci_link"

    attr(
      object,
      "epiexposure_eci_weighted_contract"
    ) <-
      "centered_profile_vs_joint_reference"

    attr(
      object,
      "epiexposure_ecilag_history_contract"
    ) <-
      "all_fitted_exposures_same_exact_max_lag_plus_one"

    attr(
      object,
      "epiexposure_ecilag_scale"
    ) <-
      "link"

    attr(
      object,
      "epiexposure_ecilag_percent_contract"
    ) <-
      "absolute_contribution_share_with_signed_contribution_retained"

    attr(
      object,
      "epiexposure_family_name"
    ) <-
      metadata$family_name

    attr(
      object,
      "epiexposure_link"
    ) <-
      metadata$link

    attr(
      object,
      "epiexposure_prediction_level"
    ) <-
      "population"

    attr(
      object,
      "epiexposure_prediction_estimand"
    ) <-
      "expected_response"

    attr(
      object,
      "epiexposure_point_prediction_contract"
    ) <-
      "central_expected_response"

    attr(
      object,
      "epiexposure_uncertainty_contract"
    ) <-
      "draw_by_draw_median_quantiles"

    attr(
      object,
      "epiexposure_ecilag_reference_values"
    ) <-
      ref_values

    attr(
      object,
      "epiexposure_ecilag_max_lag"
    ) <-
      common_max_lag

    attr(
      object,
      "epiexposure_ecilag_history_length"
    ) <-
      common_history_length

    attr(
      object,
      "epiexposure_ecilag_output"
    ) <-
      output_kind

    attr(
      object,
      "epiexposure_uncertainty"
    ) <-
      uncertainty

    if (uncertainty) {
      attr(
        object,
        "epiexposure_n_samples"
      ) <-
        n_samples

      if (identical(
        output_kind,
        "summary"
      )) {
        attr(
          object,
          "epiexposure_interval_probs"
        ) <-
          interval_probs

        attr(
          object,
          "epiexposure_summary_center"
        ) <-
          "median"
      }
    }

    object
  }

  # ==========================================================================
  # ARGUMENT MATCHING AND BASIC VALIDATION
  # ==========================================================================

  output <- match.arg(output)
  extrapolation <- match.arg(
    extrapolation
  )

  if (missing(fit) ||
      is.null(fit)) {
    stop(
      "`fit` must be provided.",
      call. = FALSE
    )
  }

  if (!valid_flag(absolute)) {
    stop(
      "`absolute` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

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
  # STRICT CURRENT EPIEXPOSURE HELPER CONTRACT
  # ==========================================================================

  required_helpers <- c(
    ".get_epiexposure_metadata",
    ".extract_central_parameters",
    ".extract_parameter_draws",
    ".build_newdata_basis",
    ".epix_standard_fixed_design",
    ".epix_validate_regular_time",
    ".epix_basis_definition",
    ".epix_handle_extrapolation",
    ".epix_cb_cols_for_var"
  )

  missing_helpers <- required_helpers[
    !vapply(
      required_helpers,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
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
      "`compute_ecilag()` requires the population/fixed-component prediction ",
      "contract.",
      call. = FALSE
    )
  }

  if (!identical(
    metadata$prediction_estimand,
    "expected_response"
  )) {
    stop(
      "`compute_ecilag()` requires the expected-response prediction estimand.",
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

  if (is.null(
    metadata$basis_objects
  )) {
    stop(
      "`compute_ecilag()` requires the exact stored cross-basis objects from ",
      "`fit_epidlnm()`. The fitted object is missing ",
      "`epiexposure_basis_objects` metadata.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # PACKAGE-WIDE EXACT COMMON LAG/HISTORY CONTRACT
  # ==========================================================================

  fitted_max_lags <- vapply(
    metadata$vars,
    function(variable) {
      current <-
        metadata$spec[[
          variable
        ]]$max_lag

      if (!is.numeric(current) ||
          length(current) != 1L ||
          is.na(current) ||
          !is.finite(current) ||
          current < 0 ||
          current != as.integer(current)) {
        stop(
          "Invalid fitted `max_lag` metadata for exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      as.integer(current)
    },
    integer(1)
  )

  if (length(
    unique(
      fitted_max_lags
    )
  ) != 1L) {
    stop(
      "All fitted exposure variables must use the same `max_lag` under the ",
      "EpiExposure exact-history contract. Fitted values are: ",
      paste(
        paste0(
          names(fitted_max_lags),
          "=",
          fitted_max_lags
        ),
        collapse = ", "
      ),
      ". Refit the model with a common lag window before using ",
      "`compute_ecilag()`.",
      call. = FALSE
    )
  }

  common_max_lag <-
    unname(
      fitted_max_lags[
        1L
      ]
    )

  common_history_length <-
    common_max_lag +
    1L

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
        "`var` must be NULL or contain unique, non-empty fitted exposure ",
        "names.",
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

  # ==========================================================================
  # NORMALIZE PROFILE INPUT OR VALIDATE LONG DATA
  # ==========================================================================

  direct_profiles <- NULL
  levels_to_use <- NULL

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

        profile <-
          profile[
            variables
          ]
      }

      for (variable in variables) {
        current <-
          profile[[
            variable
          ]]

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

        profile[[
          variable
        ]] <-
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

    supplied_lengths <- vapply(
      direct_profiles,
      length,
      integer(1)
    )

    if (length(
      unique(
        supplied_lengths
      )
    ) != 1L) {
      stop(
        "All supplied exposure profiles must have the same number of ",
        "observations. Received: ",
        paste(
          paste0(
            names(supplied_lengths),
            "=",
            supplied_lengths
          ),
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    invalid_lengths <-
      supplied_lengths !=
      common_history_length

    if (any(invalid_lengths)) {
      stop(
        "Every exposure profile must contain exactly ",
        common_history_length,
        " observations for fitted max_lag = ",
        common_max_lag,
        ". Received: ",
        paste(
          paste0(
            names(supplied_lengths),
            "=",
            supplied_lengths
          ),
          collapse = ", "
        ),
        ". Longer histories are not truncated and shorter histories are not ",
        "padded.",
        call. = FALSE
      )
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
      data[[
        group
      ]]
    )) {
      stop(
        "The grouping column cannot contain missing values.",
        call. = FALSE
      )
    }

    if (!is.numeric(
      data[[
        time
      ]]
    ) ||
    anyNA(
      data[[
        time
      ]]
    ) ||
    any(!is.finite(
      data[[
        time
      ]]
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
        data[[
          variable
        ]]
      ) ||
      anyNA(
        data[[
          variable
        ]]
      ) ||
      any(!is.finite(
        data[[
          variable
        ]]
      ))) {
        stop(
          "Exposure variable '",
          variable,
          "' must contain only finite numeric values.",
          call. = FALSE
        )
      }
    }

    .epix_validate_regular_time(
      data,
      group_index =
        as.character(
          data[[
            group
          ]]
        ),
      time_col =
        time
    )

    all_group_counts <- table(
      as.character(
        data[[
          group
        ]]
      )
    )

    invalid_group_counts <-
      all_group_counts !=
      common_history_length

    if (any(invalid_group_counts)) {
      bad <-
        all_group_counts[
          invalid_group_counts
        ]

      examples <- paste0(
        names(bad),
        "=",
        as.integer(bad)
      )

      stop(
        "Every group supplied in `data` must contain exactly ",
        common_history_length,
        " observations for fitted max_lag = ",
        common_max_lag,
        ". Non-matching group(s) include: ",
        paste(
          utils::head(
            examples,
            5L
          ),
          collapse = ", "
        ),
        if (length(examples) > 5L) {
          "; ..."
        } else {
          "."
        },
        " Longer histories are not truncated and shorter histories are not ",
        "padded.",
        call. = FALSE
      )
    }

    available_levels <-
      unique(
        data[[
          group
        ]]
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
        "No group levels are available for ECI-lag calculation.",
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
        data[[
          variable
        ]]
      ) ||
      anyNA(
        data[[
          variable
        ]]
      ) ||
      any(!is.finite(
        data[[
          variable
        ]]
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
            data[[
              variable
            ]]
          )
        },
        numeric(1)
      )

    } else if (identical(
      method,
      "percentile"
    )) {
      p <-
        ref$value

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
            data[[
              variable
            ]],
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

      if (length(
        metadata$vars
      ) > 1L) {
        warning(
          "`ref$method = 'fixed'` applies the same numeric reference value to ",
          "every fitted exposure. For multivariable models with different ",
          "exposure units, an exposure-specific named `ref` list is usually ",
          "more interpretable.",
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

    ref <-
      ref[
        metadata$vars
      ]

    ref_values <- vapply(
      metadata$vars,
      function(variable) {
        value <-
          ref[[
            variable
          ]]

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

  names(
    ref_values
  ) <-
    metadata$vars

  # ==========================================================================
  # CONSTRUCT RECORDS TO EVALUATE
  # ==========================================================================

  records <- list()
  record_index <- 1L

  if (using_profile) {
    for (variable in variables) {
      records[[
        record_index
      ]] <- list(
        variable =
          variable,
        chronology =
          as.numeric(
            direct_profiles[[
              variable
            ]]
          ),
        group_value =
          NULL
      )

      record_index <-
        record_index +
        1L
    }

  } else {
    for (current_level in levels_to_use) {
      group_rows <-
        as.character(
          data[[
            group
          ]]
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

      group_data <-
        group_data[
          order(
            group_data[[
              time
            ]]
          ),
          ,
          drop = FALSE
        ]

      if (nrow(
        group_data
      ) !=
      common_history_length) {
        stop(
          "Internal exact-history validation failed for group '",
          as.character(
            current_level
          ),
          "'.",
          call. = FALSE
        )
      }

      for (variable in variables) {
        records[[
          record_index
        ]] <- list(
          variable =
            variable,
          chronology =
            as.numeric(
              group_data[[
                variable
              ]]
            ),
          group_value =
            current_level
        )

        record_index <-
          record_index +
          1L
      }
    }
  }

  if (!length(records)) {
    stop(
      "No ECI-lag profiles are available to evaluate.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # EXTRAPOLATION CHECKS
  # ==========================================================================

  for (variable in variables) {
    variable_records <-
      records[
        vapply(
          records,
          function(current) {
            identical(
              current$variable,
              variable
            )
          },
          logical(1)
        )
      ]

    all_actual <-
      unlist(
        lapply(
          variable_records,
          `[[`,
          "chronology"
        ),
        use.names = FALSE
      )

    definition <-
      .epix_basis_definition(
        metadata,
        variable
      )

    .epix_handle_extrapolation(
      values =
        unique(
          c(
            ref_values[[
              variable
            ]],
            all_actual
          )
        ),
      exposure_range =
        definition$exposure_range,
      variable =
        variable,
      extrapolation =
        extrapolation
    )
  }

  # Check reference values for fitted non-focal exposures too.
  for (variable in metadata$vars) {
    definition <-
      .epix_basis_definition(
        metadata,
        variable
      )

    .epix_handle_extrapolation(
      values =
        ref_values[[
          variable
        ]],
      exposure_range =
        definition$exposure_range,
      variable =
        variable,
      extrapolation =
        extrapolation
    )
  }

  # ==========================================================================
  # REFERENCE PROFILE AND FIXED DESIGN
  # ==========================================================================

  reference_profiles <-
    stats::setNames(
      lapply(
        metadata$vars,
        function(variable) {
          rep(
            ref_values[[
              variable
            ]],
            common_history_length
          )
        }
      ),
      metadata$vars
    )

  reference_newdata <-
    .build_newdata_basis(
      fit = fit,
      profiles =
        reference_profiles,
      extrapolation =
        "allow"
    )

  X_reference <-
    .epix_standard_fixed_design(
      reference_newdata,
      metadata
    )

  if (nrow(
    X_reference
  ) != 1L) {
    stop(
      "Internal ECI-lag reference design must contain exactly one row.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # CROSSPRED-BASED EXACT LAG-EFFECT DESIGN CACHE
  # ==========================================================================

  if (!requireNamespace(
    "dlnm",
    quietly = TRUE
  )) {
    stop(
      "Package 'dlnm' is required by `compute_ecilag()`.",
      call. = FALSE
    )
  }

  build_effect_design_cache <- function(
    variable,
    at_values
  ) {
    basis <-
      metadata$basis_objects[[
        variable
      ]]

    if (!inherits(
      basis,
      "crossbasis"
    )) {
      stop(
        "Stored basis for exposure '",
        variable,
        "' must inherit from `crossbasis`.",
        call. = FALSE
      )
    }

    variable_columns <-
      .epix_cb_cols_for_var(
        metadata$cb_cols,
        variable
      )

    p <-
      ncol(
        basis
      )

    if (length(
      variable_columns
    ) !=
    p) {
      stop(
        "Stored canonical coefficient columns and cross-basis dimensions do ",
        "not match for exposure '",
        variable,
        "'.",
        call. = FALSE
      )
    }

    at_values <-
      unique(
        as.numeric(
          at_values
        )
      )

    if (!length(
      at_values
    ) ||
    anyNA(
      at_values
    ) ||
    any(!is.finite(
      at_values
    ))) {
      stop(
        "Internal ECI-lag evaluation values for exposure '",
        variable,
        "' are invalid.",
        call. = FALSE
      )
    }

    zero_vcov <-
      matrix(
        0,
        nrow =
          p,
        ncol =
          p
      )

    design_array <-
      array(
        NA_real_,
        dim = c(
          length(
            at_values
          ),
          common_history_length,
          p
        )
      )

    lag_index_reference <-
      NULL

    for (j in seq_len(
      p
    )) {
      unit_coef <-
        numeric(
          p
        )
      unit_coef[[
        j
      ]] <-
        1

      prediction <-
        dlnm::crosspred(
          basis =
            basis,
          coef =
            unit_coef,
          vcov =
            zero_vcov,
          model.link =
            metadata$link,
          at =
            at_values,
          cen =
            ref_values[[
              variable
            ]],
          bylag =
            1,
          cumul =
            FALSE
        )

      eta <-
        as.matrix(
          prediction$matfit
        )

      if (nrow(
        eta
      ) !=
      length(
        at_values
      ) ||
      ncol(
        eta
      ) !=
      common_history_length ||
      anyNA(
        eta
      ) ||
      any(!is.finite(
        eta
      ))) {
        stop(
          "`dlnm::crosspred()` returned an invalid lag-specific effect design ",
          "for exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      lag_index <-
        suppressWarnings(
          as.integer(
            gsub(
              "^lag",
              "",
              colnames(
                eta
              )
            )
          )
        )

      if (is.null(
        colnames(
          eta
        )
      ) ||
      anyNA(
        lag_index
      )) {
        lag_index <-
          seq.int(
            0L,
            ncol(
              eta
            ) -
              1L
          )
      }

      lag_order <-
        order(
          lag_index
        )

      eta <-
        eta[
          ,
          lag_order,
          drop = FALSE
        ]

      lag_index <-
        lag_index[
          lag_order
        ]

      if (!identical(
        as.integer(
          lag_index
        ),
        seq.int(
          0L,
          common_max_lag
        )
      )) {
        stop(
          "The lag-specific effect design for exposure '",
          variable,
          "' does not span exactly lags 0 through ",
          common_max_lag,
          ".",
          call. = FALSE
        )
      }

      if (is.null(
        lag_index_reference
      )) {
        lag_index_reference <-
          lag_index

      } else if (!identical(
        lag_index,
        lag_index_reference
      )) {
        stop(
          "Lag ordering changed while reconstructing the effect design for ",
          "exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      design_array[
        ,
        ,
        j
      ] <-
        eta
    }

    dimnames(
      design_array
    ) <- list(
      NULL,
      paste0(
        "lag",
        lag_index_reference
      ),
      variable_columns
    )

    list(
      at_values =
        at_values,
      lag =
        as.integer(
          lag_index_reference
        ),
      design =
        design_array,
      cb_cols =
        variable_columns
    )
  }

  effect_cache <-
    stats::setNames(
      vector(
        "list",
        length(
          variables
        )
      ),
      variables
    )

  for (variable in variables) {
    variable_records <-
      records[
        vapply(
          records,
          function(current) {
            identical(
              current$variable,
              variable
            )
          },
          logical(1)
        )
      ]

    actual_values <-
      unique(
        unlist(
          lapply(
            variable_records,
            `[[`,
            "chronology"
          ),
          use.names = FALSE
        )
      )

    effect_cache[[
      variable
    ]] <-
      build_effect_design_cache(
        variable =
          variable,
        at_values =
          actual_values
      )
  }

  # ==========================================================================
  # CENTRAL PARAMETERS AND ONE JOINT PARAMETER-DRAW MATRIX
  # ==========================================================================

  central_parameters <-
    .extract_central_parameters(
      fit
    )

  if (!identical(
    names(
      central_parameters
    ),
    colnames(
      X_reference
    )
  )) {
    stop(
      "Central parameter names do not align exactly with the ECI-lag fixed ",
      "design.",
      call. = FALSE
    )
  }

  parameter_draws <-
    NULL

  if (uncertainty) {
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
      set.seed(
        seed
      )
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
    anyNA(
      parameter_draws
    ) ||
    any(!is.finite(
      parameter_draws
    ))) {
      stop(
        "Parameter draws do not satisfy the current EpiExposure ECI-lag ",
        "fixed-design contract.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # PROCESS ONE PROFILE-VARIABLE RECORD
  # ==========================================================================

  process_record <- function(
    current
  ) {
    variable <-
      current$variable

    chronology <-
      as.numeric(
        current$chronology
      )

    if (length(
      chronology
    ) !=
    common_history_length) {
      stop(
        "Internal ECI-lag history for exposure '",
        variable,
        "' does not contain exactly ",
        common_history_length,
        " values.",
        call. = FALSE
      )
    }

    reference_value <-
      ref_values[[
        variable
      ]]

    lag_for_position <-
      seq.int(
        common_max_lag,
        0L
      )

    cache <-
      effect_cache[[
        variable
      ]]

    p_var <-
      length(
        cache$cb_cols
      )

    contribution_design <-
      matrix(
        NA_real_,
        nrow =
          common_history_length,
        ncol =
          p_var,
        dimnames = list(
          NULL,
          cache$cb_cols
        )
      )

    for (position in seq_len(
      common_history_length
    )) {
      lag_value <-
        lag_for_position[[
          position
        ]]

      lag_column <-
        match(
          lag_value,
          cache$lag
        )

      actual_value <-
        chronology[[
          position
        ]]

      actual_row <-
        match(
          actual_value,
          cache$at_values
        )

      if (is.na(
        actual_row
      ) ||
      is.na(
        lag_column
      )) {
        stop(
          "Internal ECI-lag value/lag matching failed for exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      contribution_design[
        position,
      ] <-
        cache$design[
          actual_row,
          lag_column,
          ,
          drop = TRUE
        ]
    }

    if (anyNA(
      contribution_design
    ) ||
    any(!is.finite(
      contribution_design
    ))) {
      stop(
        "The reconstructed ECI-lag contribution design contains non-finite ",
        "values for exposure '",
        variable,
        "'.",
        call. = FALSE
      )
    }

    # ------------------------------------------------------------------------
    # Independently reconstruct the complete trajectory contrast using the
    # centralized prediction-basis helper. This must equal the sum of the exact
    # lag-specific centered effect designs.
    # ------------------------------------------------------------------------

    full_profiles <-
      reference_profiles

    full_profiles[[
      variable
    ]] <-
      chronology

    full_newdata <-
      .build_newdata_basis(
        fit = fit,
        profiles =
          full_profiles,
        extrapolation =
          "allow"
      )

    X_full <-
      .epix_standard_fixed_design(
        full_newdata,
        metadata
      )

    if (nrow(
      X_full
    ) != 1L ||
    !identical(
      colnames(
        X_full
      ),
      colnames(
        X_reference
      )
    )) {
      stop(
        "Internal full-profile ECI-lag design is inconsistent with the ",
        "reference design.",
        call. = FALSE
      )
    }

    full_contrast_design <-
      as.numeric(
        X_full[
          1L,
          ,
          drop = TRUE
        ] -
          X_reference[
            1L,
            ,
            drop = TRUE
          ]
      )

    names(
      full_contrast_design
    ) <-
      colnames(
        X_reference
      )

    nonfocal_names <-
      setdiff(
        names(
          full_contrast_design
        ),
        cache$cb_cols
      )

    nonfocal_change <-
      abs(
        full_contrast_design[
          nonfocal_names
        ]
      )

    if (length(
      nonfocal_change
    ) &&
    any(
      nonfocal_change >
      1000 *
      .Machine$double.eps
    )) {
      bad_names <-
        nonfocal_names[
          nonfocal_change >
            1000 *
            .Machine$double.eps
        ]

      stop(
        "Internal ECI-lag target construction changed non-focal fixed terms: ",
        paste(
          bad_names,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    summed_lag_design <-
      colSums(
        contribution_design
      )

    focal_full_design <-
      full_contrast_design[
        cache$cb_cols
      ]

    design_difference <-
      summed_lag_design -
      focal_full_design

    design_tolerance <-
      1e4 *
      .Machine$double.eps *
      pmax(
        1,
        abs(
          summed_lag_design
        ),
        abs(
          focal_full_design
        )
      )

    if (any(
      abs(
        design_difference
      ) >
      design_tolerance
    )) {
      stop(
        "Exact lag-specific effect designs do not sum to the complete centered ",
        "trajectory design for exposure '",
        variable,
        "'. This indicates an internal cross-basis decomposition mismatch.",
        call. = FALSE
      )
    }

    # ------------------------------------------------------------------------
    # Deterministic central-parameter values
    # ------------------------------------------------------------------------

    beta_focal <-
      central_parameters[
        cache$cb_cols
      ]

    contribution_point <-
      as.numeric(
        contribution_design %*%
          beta_focal
      )

    eci_weighted_point <-
      as.numeric(
        sum(
          full_contrast_design *
            central_parameters
        )
      )

    if (!nearly_equal(
      sum(
        contribution_point
      ),
      eci_weighted_point
    )) {
      stop(
        "Central lag-specific contributions do not sum to the centered ",
        "`ECI_weighted` for exposure '",
        variable,
        "'.",
        call. = FALSE
      )
    }

    denominator_point <-
      sum(
        abs(
          contribution_point
        )
      )

    percent_point <- if (
      is.finite(
        denominator_point
      ) &&
      denominator_point >
      0
    ) {
      100 *
        abs(
          contribution_point
        ) /
        denominator_point
    } else {
      rep(
        NA_real_,
        common_history_length
      )
    }

    eci_raw <-
      sum(
        chronology
      )

    eci_raw_centered <-
      sum(
        chronology -
          reference_value
      )

    # ------------------------------------------------------------------------
    # Deterministic output
    # ------------------------------------------------------------------------

    if (!uncertainty) {
      by_lag <-
        data.frame(
          var =
            variable,
          lag =
            lag_for_position,
          exposure =
            chronology,
          reference_value =
            rep(
              reference_value,
              common_history_length
            ),
          exposure_minus_reference =
            chronology -
            reference_value,
          contribution =
            contribution_point,
          percent_contribution =
            percent_point,
          stringsAsFactors = FALSE
        )

      if (absolute) {
        by_lag$abs_contribution <-
          abs(
            by_lag$contribution
          )
      }

      by_lag <-
        by_lag[
          order(
            by_lag$lag
          ),
          ,
          drop = FALSE
        ]

      rownames(
        by_lag
      ) <-
        NULL

      summary_row <-
        data.frame(
          var =
            variable,
          reference_value =
            reference_value,
          max_lag =
            common_max_lag,
          n_exposure_values =
            common_history_length,
          ECI_raw =
            eci_raw,
          ECI_raw_centered =
            eci_raw_centered,
          ECI_weighted =
            eci_weighted_point,
          stringsAsFactors = FALSE
        )

      return(
        list(
          summary =
            summary_row,
          by_lag =
            by_lag,
          by_lag_samples =
            NULL,
          ECI_weighted_samples =
            NULL
        )
      )
    }

    # ------------------------------------------------------------------------
    # Draw-level quantities. The same global parameter_draws matrix is reused
    # for every record, preserving joint parameter covariance across the entire
    # function result.
    # ------------------------------------------------------------------------

    focal_draws <-
      parameter_draws[
        ,
        cache$cb_cols,
        drop = FALSE
      ]

    contribution_draws <-
      focal_draws %*%
      t(
        contribution_design
      )

    eci_weighted_draws <-
      as.numeric(
        parameter_draws %*%
          full_contrast_design
      )

    if (anyNA(
      contribution_draws
    ) ||
    any(!is.finite(
      contribution_draws
    )) ||
    anyNA(
      eci_weighted_draws
    ) ||
    any(!is.finite(
      eci_weighted_draws
    ))) {
      stop(
        "Non-finite draw-level ECI-lag quantities were produced for exposure '",
        variable,
        "'.",
        call. = FALSE
      )
    }

    row_sum_contribution_draws <-
      rowSums(
        contribution_draws
      )

    draw_tolerance <-
      1e4 *
      .Machine$double.eps *
      pmax(
        1,
        abs(
          row_sum_contribution_draws
        ),
        abs(
          eci_weighted_draws
        )
      )

    if (any(
      abs(
        row_sum_contribution_draws -
        eci_weighted_draws
      ) >
      draw_tolerance
    )) {
      stop(
        "Draw-level lag-specific contributions do not sum to draw-level ",
        "`ECI_weighted` for exposure '",
        variable,
        "'.",
        call. = FALSE
      )
    }

    absolute_denominator_draws <-
      rowSums(
        abs(
          contribution_draws
        )
      )

    percent_draws <-
      matrix(
        NA_real_,
        nrow =
          n_samples,
        ncol =
          common_history_length
      )

    defined_percent_draws <-
      is.finite(
        absolute_denominator_draws
      ) &
      absolute_denominator_draws >
      0

    if (any(
      defined_percent_draws
    )) {
      percent_draws[
        defined_percent_draws,
      ] <-
        100 *
        abs(
          contribution_draws[
            defined_percent_draws,
            ,
            drop = FALSE
          ]
        ) /
        absolute_denominator_draws[
          defined_percent_draws
        ]
    }

    if (sum(
      defined_percent_draws
    ) <
    n_samples) {
      warning(
        "Percent lag contribution is undefined for ",
        n_samples -
          sum(
            defined_percent_draws
          ),
        " parameter draw(s) for exposure '",
        variable,
        "' because every lag contribution was zero in those draw(s). ",
        "Those draws are retained as NA for `percent_contribution` and are not ",
        "silently converted to zero.",
        call. = FALSE
      )
    }

    # ------------------------------------------------------------------------
    # Summarize draw-level lag quantities
    # ------------------------------------------------------------------------

    contribution_summary <-
      lapply(
        seq_len(
          common_history_length
        ),
        function(j) {
          summarise_finite(
            contribution_draws[
              ,
              j
            ],
            interval_probs
          )
        }
      )

    percent_summary <-
      lapply(
        seq_len(
          common_history_length
        ),
        function(j) {
          summarise_defined(
            percent_draws[
              ,
              j
            ],
            interval_probs
          )
        }
      )

    contribution_summary <-
      do.call(
        rbind,
        contribution_summary
      )

    percent_summary <-
      do.call(
        rbind,
        percent_summary
      )

    by_lag <-
      data.frame(
        var =
          variable,
        lag =
          lag_for_position,
        exposure =
          chronology,
        reference_value =
          rep(
            reference_value,
            common_history_length
          ),
        exposure_minus_reference =
          chronology -
          reference_value,
        contribution =
          contribution_summary[
            ,
            "estimate"
          ],
        contribution_sd =
          contribution_summary[
            ,
            "sd"
          ],
        contribution_lower =
          contribution_summary[
            ,
            "lower"
          ],
        contribution_upper =
          contribution_summary[
            ,
            "upper"
          ],
        percent_contribution =
          percent_summary[
            ,
            "estimate"
          ],
        percent_contribution_sd =
          percent_summary[
            ,
            "sd"
          ],
        percent_contribution_lower =
          percent_summary[
            ,
            "lower"
          ],
        percent_contribution_upper =
          percent_summary[
            ,
            "upper"
          ],
        percent_contribution_n_defined =
          as.integer(
            percent_summary[
              ,
              "n_defined"
            ]
          ),
        stringsAsFactors = FALSE
      )

    if (absolute) {
      abs_contribution_summary <-
        do.call(
          rbind,
          lapply(
            seq_len(
              common_history_length
            ),
            function(j) {
              summarise_finite(
                abs(
                  contribution_draws[
                    ,
                    j
                  ]
                ),
                interval_probs
              )
            }
          )
        )

      by_lag$abs_contribution <-
        abs_contribution_summary[
          ,
          "estimate"
        ]

      by_lag$abs_contribution_sd <-
        abs_contribution_summary[
          ,
          "sd"
        ]

      by_lag$abs_contribution_lower <-
        abs_contribution_summary[
          ,
          "lower"
        ]

      by_lag$abs_contribution_upper <-
        abs_contribution_summary[
          ,
          "upper"
        ]
    }

    by_lag <-
      by_lag[
        order(
          by_lag$lag
        ),
        ,
        drop = FALSE
      ]

    rownames(
      by_lag
    ) <-
      NULL

    eci_summary <-
      summarise_finite(
        eci_weighted_draws,
        interval_probs
      )

    summary_row <-
      data.frame(
        var =
          variable,
        reference_value =
          reference_value,
        max_lag =
          common_max_lag,
        n_exposure_values =
          common_history_length,
        ECI_raw =
          eci_raw,
        ECI_raw_centered =
          eci_raw_centered,
        ECI_weighted =
          unname(
            eci_summary[
              "estimate"
            ]
          ),
        ECI_weighted_sd =
          unname(
            eci_summary[
              "sd"
            ]
          ),
        ECI_weighted_lower =
          unname(
            eci_summary[
              "lower"
            ]
          ),
        ECI_weighted_upper =
          unname(
            eci_summary[
              "upper"
            ]
          ),
        stringsAsFactors = FALSE
      )

    # ------------------------------------------------------------------------
    # Optional sample-level tables
    # ------------------------------------------------------------------------

    by_lag_samples <-
      NULL
    weighted_samples <-
      NULL

    if (identical(
      output,
      "samples"
    )) {
      sample_blocks <-
        vector(
          "list",
          n_samples
        )

      for (sample_index in seq_len(
        n_samples
      )) {
        sample_data <-
          data.frame(
            sample =
              sample_index,
            var =
              variable,
            lag =
              lag_for_position,
            exposure =
              chronology,
            reference_value =
              rep(
                reference_value,
                common_history_length
              ),
            exposure_minus_reference =
              chronology -
              reference_value,
            contribution =
              contribution_draws[
                sample_index,
                ,
                drop = TRUE
              ],
            percent_contribution =
              percent_draws[
                sample_index,
                ,
                drop = TRUE
              ],
            stringsAsFactors = FALSE
          )

        if (absolute) {
          sample_data$abs_contribution <-
            abs(
              sample_data$contribution
            )
        }

        sample_data <-
          sample_data[
            order(
              sample_data$lag
            ),
            ,
            drop = FALSE
          ]

        sample_blocks[[
          sample_index
        ]] <-
          sample_data
      }

      by_lag_samples <-
        do.call(
          rbind,
          sample_blocks
        )

      rownames(
        by_lag_samples
      ) <-
        NULL

      weighted_samples <-
        data.frame(
          sample =
            seq_len(
              n_samples
            ),
          ECI_weighted =
            eci_weighted_draws,
          stringsAsFactors = FALSE
        )
    }

    list(
      summary =
        summary_row,
      by_lag =
        by_lag,
      by_lag_samples =
        by_lag_samples,
      ECI_weighted_samples =
        weighted_samples
    )
  }

  # ==========================================================================
  # PROCESS EVERY RECORD
  # ==========================================================================

  results <-
    lapply(
      records,
      process_record
    )

  # ==========================================================================
  # ADD GROUP IDENTIFIERS
  # ==========================================================================

  if (using_data) {
    for (i in seq_along(
      results
    )) {
      current_group <-
        records[[
          i
        ]]$group_value

      results[[
        i
      ]]$summary[[
        group
      ]] <-
        current_group

      results[[
        i
      ]]$summary <-
        results[[
          i
        ]]$summary[
          ,
          c(
            group,
            setdiff(
              names(
                results[[
                  i
                ]]$summary
              ),
              group
            )
          ),
          drop = FALSE
        ]

      results[[
        i
      ]]$by_lag[[
        group
      ]] <-
        current_group

      results[[
        i
      ]]$by_lag <-
        results[[
          i
        ]]$by_lag[
          ,
          c(
            group,
            setdiff(
              names(
                results[[
                  i
                ]]$by_lag
              ),
              group
            )
          ),
          drop = FALSE
        ]

      if (uncertainty &&
          identical(
            output,
            "samples"
          )) {
        results[[
          i
        ]]$by_lag_samples[[
          group
        ]] <-
          current_group

        results[[
          i
        ]]$by_lag_samples <-
          results[[
            i
          ]]$by_lag_samples[
            ,
            c(
              group,
              setdiff(
                names(
                  results[[
                    i
                  ]]$by_lag_samples
                ),
                group
              )
            ),
            drop = FALSE
          ]

        results[[
          i
        ]]$ECI_weighted_samples[[
          group
        ]] <-
          current_group

        results[[
          i
        ]]$ECI_weighted_samples <-
          results[[
            i
          ]]$ECI_weighted_samples[
            ,
            c(
              group,
              setdiff(
                names(
                  results[[
                    i
                  ]]$ECI_weighted_samples
                ),
                group
              )
            ),
            drop = FALSE
          ]
      }
    }
  }

  # ==========================================================================
  # SINGLE COMBINATION: PRESERVE THE HISTORICAL LIST SURFACE
  # ==========================================================================

  if (length(
    results
  ) == 1L &&
  !using_data) {
    one <-
      results[[
        1L
      ]]

    summary_row <-
      one$summary

    out <-
      list(
        ECI_raw =
          summary_row$ECI_raw[[
            1L
          ]],
        ECI_raw_centered =
          summary_row$ECI_raw_centered[[
            1L
          ]],
        ECI_weighted =
          summary_row$ECI_weighted[[
            1L
          ]],
        reference_value =
          summary_row$reference_value[[
            1L
          ]],
        max_lag =
          summary_row$max_lag[[
            1L
          ]],
        n_exposure_values =
          summary_row$n_exposure_values[[
            1L
          ]],
        by_lag =
          one$by_lag
      )

    if (uncertainty) {
      out$ECI_weighted_sd <-
        summary_row$ECI_weighted_sd[[
          1L
        ]]

      out$ECI_weighted_lower <-
        summary_row$ECI_weighted_lower[[
          1L
        ]]

      out$ECI_weighted_upper <-
        summary_row$ECI_weighted_upper[[
          1L
        ]]

      if (identical(
        output,
        "samples"
      )) {
        out$by_lag_samples <-
          one$by_lag_samples

        out$ECI_weighted_samples <-
          one$ECI_weighted_samples
      }
    }

    out <-
      attach_contract_attributes(
        out,
        metadata =
          metadata,
        ref_values =
          ref_values,
        common_max_lag =
          common_max_lag,
        common_history_length =
          common_history_length,
        output_kind =
          if (uncertainty) {
            output
          } else {
            "summary"
          }
      )

    return(out)
  }

  # ==========================================================================
  # MULTIPLE GROUPS / VARIABLES
  # ==========================================================================

  eci_summary <-
    do.call(
      rbind,
      lapply(
        results,
        `[[`,
        "summary"
      )
    )

  rownames(
    eci_summary
  ) <-
    NULL

  by_lag <-
    do.call(
      rbind,
      lapply(
        results,
        `[[`,
        "by_lag"
      )
    )

  rownames(
    by_lag
  ) <-
    NULL

  out <-
    list(
      eci_summary =
        eci_summary,
      by_lag =
        by_lag
    )

  if (uncertainty &&
      identical(
        output,
        "samples"
      )) {
    by_lag_samples <-
      do.call(
        rbind,
        lapply(
          results,
          `[[`,
          "by_lag_samples"
        )
      )

    rownames(
      by_lag_samples
    ) <-
      NULL

    weighted_samples <-
      vector(
        "list",
        length(
          results
        )
      )

    for (i in seq_along(
      results
    )) {
      current <-
        results[[
          i
        ]]$ECI_weighted_samples

      current$var <-
        records[[
          i
        ]]$variable

      if (using_data) {
        current[[
          group
        ]] <-
          records[[
            i
          ]]$group_value

        current <-
          current[
            ,
            c(
              group,
              "var",
              "sample",
              "ECI_weighted"
            ),
            drop = FALSE
          ]

      } else {
        current <-
          current[
            ,
            c(
              "var",
              "sample",
              "ECI_weighted"
            ),
            drop = FALSE
          ]
      }

      weighted_samples[[
        i
      ]] <-
        current
    }

    weighted_samples <-
      do.call(
        rbind,
        weighted_samples
      )

    rownames(
      weighted_samples
    ) <-
      NULL

    out$by_lag_samples <-
      by_lag_samples

    out$ECI_weighted_samples <-
      weighted_samples
  }

  out <-
    attach_contract_attributes(
      out,
      metadata =
        metadata,
      ref_values =
        ref_values,
      common_max_lag =
        common_max_lag,
      common_history_length =
        common_history_length,
      output_kind =
        if (uncertainty) {
          output
        } else {
          "summary"
        }
    )

  out
}
