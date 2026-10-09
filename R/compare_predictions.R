#' Compare predicted outcomes between exposure scenarios
#'
#' Compares population-level expected predictions from two or more explicit
#' exposure scenarios using `predict_outcomes()` as the computational backend.
#'
#' All scenarios are submitted to `predict_outcomes()` in **one joint batch**.
#' This is essential when `uncertainty = TRUE`: the same fixed/population
#' parameter draw is applied to every scenario before differences, ratios, or
#' other contrasts are calculated.
#'
#' Exposure profiles follow the EpiExposure chronological convention: each
#' profile is supplied from the oldest/earliest exposure observation to the most
#' recent observation, with the final value corresponding to lag 0.
#'
#' @param fit Fitted model returned by the current `fit_epidlnm()`.
#' @param profiles Named list containing at least two scenarios. Top-level names
#'   are scenario labels.
#'
#'   Each scenario must contain the complete fitted exposure profile structure
#'   expected by `predict_outcomes()`. For each fitted exposure, the scenario may
#'   supply:
#'
#'   \itemize{
#'     \item one finite numeric scalar, interpreted as a constant profile and
#'       expanded internally to that exposure's fitted `max_lag + 1` positions;
#'     \item one explicit numeric chronological profile;
#'     \item a matrix/data frame whose rows are multiple chronological profiles;
#'     \item a list of multiple numeric chronological profiles;
#'     \item an object returned directly by the audited `simulate_exposures()`.
#'   }
#'
#'   For an audited `simulate_exposures()` object, all `n` simulated exposure
#'   profiles in its canonical `$profiles` collection are used directly.
#'
#'   Thus, for example:
#'
#'   ```
#'   scenario_A <- list(
#'     tmean = 25,
#'     rain = simulate_exposures(..., n = 100),
#'     wetness = simulate_exposures(..., n = 100)
#'   )
#'   ```
#'
#'   creates 100 multivariable profiles: the constant temperature profile is
#'   recycled, while rain profile 1 is paired with wetness profile 1, rain
#'   profile 2 with wetness profile 2, and so forth.
#'
#'   Within a scenario, fitted exposures are paired by profile index exactly as
#'   in `predict_outcomes()`. Across scenarios, the number of matched profile
#'   sets must be either one or a common maximum. A scenario containing one
#'   profile set is recycled across that common profile index when another
#'   scenario contains multiple matched profile sets.
#'
#'   Thus multiple profiles are compared by matched index; a Cartesian product
#'   of scenario profiles is never created.
#'
#'   Explicit non-scalar profile must contain exactly the fitted
#'   `max_lag + 1` positions for their exposure. A `simulate_exposures()` object
#'   must have been generated with the same `max_lag` stored for that fitted
#'   exposure. This strict check prevents extra or missing positions in the profile
#'   from being silently ignored.
#' @param profiles1 Legacy first scenario profile. Used only with `profiles2`
#'   when `profiles = NULL`.
#' @param profiles2 Legacy second scenario profile. Used only with `profiles1`
#'   when `profiles = NULL`.
#' @param type Prediction scale:
#'
#'   - `"response"` (default): population-level expected outcome on its natural
#'     response scale;
#'   - `"link"`: population-level linear predictor.
#'
#'   `diff = pred2 - pred1` is available on either scale. Relative metrics
#'   (`ratio` and `percent_change`) are calculated only for response-scale
#'   predictions.
#' @param uncertainty Logical. If `FALSE`, comparisons use deterministic
#'   predictions from the central fixed/population parameter estimate. If
#'   `TRUE`, comparisons are calculated draw by draw from joint parameter
#'   uncertainty.
#' @param output Character. `"summary"` or `"samples"`.
#'
#'   With `uncertainty = FALSE`, only `"summary"` is valid.
#'
#'   With `uncertainty = TRUE`, `"samples"` returns draw-level pairwise
#'   comparisons and `"summary"` summarizes the **draw-level comparisons** by
#'   median, SD, and empirical interval.
#' @param n_samples Positive integer number of fixed/population parameter draws
#'   used when `uncertainty = TRUE`. At least two are required.
#' @param probs Numeric vector of length two defining the lower and upper
#'   empirical uncertainty probabilities. The default `c(0.025, 0.975)` gives
#'   a 95 percent interval.
#' @param seed Optional finite integer seed passed once to the joint
#'   `predict_outcomes()` call. It controls **model-parameter draws only** and
#'   does not resimulate exposure profiles already produced by
#'   `simulate_exposures()`. The prediction backend restores the caller's
#'   random-number state after seeded sampling.
#' @param extrapolation Behavior when any scenario contains exposure values
#'   outside the fitted exposure range: `"warn"` (default), `"error"`, or
#'   `"allow"`. The fitted DLNM basis is never re-estimated from scenario data.
#' @param eps Positive finite tolerance used only to decide whether a
#'   response-scale reference prediction is sufficiently above zero for a
#'   relative comparison.
#'
#'   No epsilon is added to a prediction. When `pred1 <= eps`, `ratio` and
#'   `percent_change` are returned as `NA` rather than modifying the
#'   denominator.
#'
#' @return A data frame containing pairwise scenario comparisons.
#'
#'   Deterministic output contains:
#'
#'   - `scenario1`, `scenario2`;
#'   - `profile` when more than one matched profile set is compared;
#'   - `pred1`, `pred2`;
#'   - `diff = pred2 - pred1`;
#'   - `abs_diff = abs(diff)`;
#'   - `ratio = pred2 / pred1` when defined on the response scale;
#'   - `percent_change = 100 * (ratio - 1)` when defined;
#'   - `relative_change_defined`;
#'   - `percentage_point_change = 100 * diff` for Beta and Binomial
#'     response-scale predictions, otherwise `NA`.
#'
#'   With `uncertainty = TRUE` and `output = "samples"`, the same quantities are
#'   returned for every matched parameter draw, together with `sample`.
#'
#'   With `uncertainty = TRUE` and `output = "summary"`, `pred1`, `pred2`,
#'   `diff`, `abs_diff`, and any applicable relative/probability contrasts are
#'   summarized by their median plus `_sd`, `_lower`, and `_upper` columns.
#'
#'   For relative metrics, `relative_defined_fraction` gives the fraction of
#'   draws in which `pred1 > eps`. A posterior relative summary is returned only
#'   when the relative metric is defined for **all** matched draws; otherwise
#'   its median and interval are `NA` rather than being calculated from a
#'   truncated subset of draws.
#'
#'   Output attributes record the population expected-response prediction
#'   contract, family, link, prediction scale, scenario order, profile-input
#'   sources, uncertainty contract, and comparison direction.
#'
#' @details
#' ## EpiExposure v1 prediction target
#'
#' `compare_predictions()` compares the same prediction estimand used throughout
#' the current package:
#'
#' \deqn{E(Y \mid X, \theta)}
#'
#' using the fixed/population component with fitted random effects excluded.
#' Models may have been fitted with random effects, but this function does not
#' provide conditional/group-specific prediction.
#'
#' For deterministic prediction, each scenario is evaluated with the central
#' parameter estimate:
#'
#' \deqn{\hat\mu_j = g^{-1}(X_j\hat\beta)}
#'
#' for `type = "response"`, or \eqn{X_j\hat\beta} for `type = "link"`.
#'
#' No residual, observation, process, dispersion, or posterior-predictive noise
#' is added.
#'
#' ## Pairwise direction
#'
#' For every ordered pair generated from the input scenario order,
#'
#' \deqn{diff = pred_2 - pred_1.}
#'
#' A positive difference therefore means that scenario 2 has the larger
#' predicted expected outcome on the selected prediction scale.
#'
#' ## Relative response-scale contrasts
#'
#' When `type = "response"` and `pred1 > eps`,
#'
#' \deqn{ratio = pred_2/pred_1}
#'
#' and
#'
#' \deqn{percent\_change = 100(ratio - 1).}
#'
#' Earlier EpiExposure code used `pred2 / (pred1 + eps)` and divided differences
#' by `abs(pred1) + eps`. Those formulas changed valid predictions and gave
#' unusual behavior for non-positive denominators. The current implementation
#' never modifies predicted values.
#'
#' A numerical response ratio is scientifically interpretable as a relative
#' outcome change only when the outcome scale has an appropriate meaningful
#' zero. This is usually natural for expected counts, positive means, and
#' probabilities, but may not be meaningful for every Gaussian response.
#'
#' With `type = "link"`, `ratio` and `percent_change` are intentionally `NA`.
#' Dividing two linear predictors generally does not define a meaningful model
#' contrast. The valid link-scale contrast is `diff`.
#'
#' ## Probability and proportion outcomes
#'
#' For Beta and Binomial response-scale predictions,
#'
#' \deqn{percentage\_point\_change = 100(pred_2 - pred_1).}
#'
#' This is distinct from `percent_change`. For example, a probability increase
#' from 0.20 to 0.30 is +10 percentage points but +50 percent relative to 0.20.
#'
#' ## Joint uncertainty across scenarios
#'
#' A central requirement of scenario comparison is that scenario contrasts use
#' the same parameter draw. For draw \eqn{s},
#'
#' \deqn{
#'   diff^{(s)} =
#'   \mu_2^{(s)} - \mu_1^{(s)}.
#' }
#'
#' `compare_predictions()` therefore constructs one combined prediction batch
#' containing every scenario and calls `predict_outcomes()` only once. The
#' prediction backend applies one joint \eqn{\beta^{(s)}} to every row in that
#' batch.
#'
#' The function then matches scenario 1 and scenario 2 within the same
#' `profile` and `sample`, calculates all contrasts draw by draw, and only then
#' summarizes them. This preserves covariance among scenario predictions.
#'
#' Calling `predict_outcomes()` independently for each scenario and pairing
#' rows with the same integer sample label would not establish that those rows
#' came from the same parameter draw and is therefore not used.
#'
#' ## Direct `simulate_exposures()` integration
#'
#' An audited `simulate_exposures()` object contains the canonical `$profiles`
#' collection and metadata describing chronological order and `max_lag`.
#' `compare_predictions()` recognizes this object directly and uses all `n`
#' simulated exposure profiles.
#'
#' For example:
#'
#' ```
#' scenario_A <- list(
#'   tmean = simulate_exposures(
#'     max_lag = 85,
#'     n = 100,
#'     mode = "profile",
#'     background = list(dist = "normal", mean = 25, sd = 4)
#'   ),
#'   rain = 5,
#'   wetness = simulate_exposures(
#'     max_lag = 85,
#'     n = 100,
#'     mode = "profile",
#'     background = list(dist = "normal", mean = 10, sd = 2)
#'   )
#' )
#' ```
#'
#' The 100 tmean and wetness profiles are paired by index; the constant rain
#' profile is recycled across them.
#'
#' `n` in `simulate_exposures()` controls exposure-profile variability only.
#' Model-parameter uncertainty requested through `uncertainty = TRUE` in this
#' function is a separate layer. Thus 100 simulated exposure profiles and 1000
#' coefficient draws represent different sources of variation.
#'
#' The direct integration also supports the immediately preceding
#' `simulate_exposures()` return structure for backward compatibility, provided
#' its metadata explicitly identify chronological profile order. New package
#' code should use the standardized audited object.
#'
#' ## Multiple profile sets within scenarios
#'
#' If all scenarios contain one profile set, the output has one row per
#' scenario pair.
#'
#' If at least one scenario contains multiple matched profile sets, comparisons
#' are made by profile index. A scenario with one profile set is recycled to the
#' common number of profile sets; any other non-singleton scenario must contain
#' that same common number.
#'
#' @export
compare_predictions <- function(
    fit,
    profiles = NULL,
    profiles1 = NULL,
    profiles2 = NULL,
    type = c("response", "link"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    probs = c(0.025, 0.975),
    seed = NULL,
    extrapolation = c("warn", "error", "allow"),
    eps = 1e-12
) {

  # ARGUMENT MATCHING AND SMALL VALIDATORS

  type <- match.arg(type)
  output <- match.arg(output)
  extrapolation <- match.arg(extrapolation)

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

  validate_probs <- function(x) {
    if (!is.numeric(x) ||
        length(x) != 2L ||
        anyNA(x) ||
        any(!is.finite(x)) ||
        any(x <= 0 | x >= 1) ||
        x[1L] >= x[2L]) {
      stop(
        "`probs` must contain two finite probabilities satisfying ",
        "0 < probs[1] < probs[2] < 1.",
        call. = FALSE
      )
    }

    as.numeric(x)
  }

  summarise_numeric <- function(x, probs_use) {
    if (!is.numeric(x) ||
        !length(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      return(
        c(
          estimate = NA_real_,
          sd = NA_real_,
          lower = NA_real_,
          upper = NA_real_
        )
      )
    }

    c(
      estimate = stats::median(x),
      sd = if (length(x) > 1L) {
        stats::sd(x)
      } else {
        NA_real_
      },
      lower = stats::quantile(
        x,
        probs = probs_use[1L],
        names = FALSE,
        type = 7
      ),
      upper = stats::quantile(
        x,
        probs = probs_use[2L],
        names = FALSE,
        type = 7
      )
    )
  }

  # STRICT FIT CONTRACT

  if (is.null(fit)) {
    stop(
      "`fit` cannot be NULL.",
      call. = FALSE
    )
  }

  if (!exists(
    ".get_epiexposure_metadata",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "Internal EpiExposure prediction helper `.get_epiexposure_metadata()` ",
      "is unavailable. Load the complete current EpiExposure package before ",
      "calling `compare_predictions()`.",
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
      "`compare_predictions()` requires population/fixed-component prediction.",
      call. = FALSE
    )
  }

  if (!identical(
    metadata$prediction_estimand,
    "expected_response"
  )) {
    stop(
      "`compare_predictions()` requires the expected-response estimand.",
      call. = FALSE
    )
  }

  if (!identical(
    metadata$point_prediction_contract,
    "central_expected_response"
  )) {
    stop(
      "The fitted model does not satisfy the EpiExposure central point-",
      "prediction contract.",
      call. = FALSE
    )
  }

  if (!identical(
    metadata$uncertainty_contract,
    "draw_by_draw_median_quantiles"
  )) {
    stop(
      "The fitted model does not satisfy the EpiExposure draw-by-draw ",
      "uncertainty contract.",
      call. = FALSE
    )
  }

  family_name <- metadata$family_name
  link_name <- metadata$link
  fitted_vars <- metadata$vars
  fitted_spec <- metadata$spec

  supported_families <- c(
    "beta",
    "binomial",
    "poisson",
    "gamma",
    "gaussian",
    "negative_binomial"
  )

  if (!family_name %in%
      supported_families) {
    stop(
      "Unsupported EpiExposure v1 family: '",
      family_name,
      "'.",
      call. = FALSE
    )
  }

  # GENERAL ARGUMENT VALIDATION

  if (!valid_flag(uncertainty)) {
    stop(
      "`uncertainty` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!uncertainty &&
      identical(output, "samples")) {
    stop(
      "`output = 'samples'` requires `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  if (!valid_integer_scalar(n_samples) ||
      n_samples <= 0) {
    stop(
      "`n_samples` must be one positive finite integer.",
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

  probs <- validate_probs(probs)

  if (!is.null(seed)) {
    if (!valid_integer_scalar(seed)) {
      stop(
        "`seed` must be NULL or one finite integer.",
        call. = FALSE
      )
    }
    seed <- as.integer(seed)
  }

  if (!is.numeric(eps) ||
      length(eps) != 1L ||
      is.na(eps) ||
      !is.finite(eps) ||
      eps <= 0) {
    stop(
      "`eps` must be one positive finite numeric tolerance.",
      call. = FALSE
    )
  }

  # LEGACY TWO-SCENARIO INPUT

  if (is.null(profiles)) {
    if (is.null(profiles1) ||
        is.null(profiles2)) {
      stop(
        "Provide either `profiles` or both legacy arguments `profiles1` and ",
        "`profiles2`.",
        call. = FALSE
      )
    }

    profiles <- list(
      scenario1 = profiles1,
      scenario2 = profiles2
    )

  } else if (
    !is.null(profiles1) ||
    !is.null(profiles2)
  ) {
    stop(
      "When `profiles` is supplied, legacy arguments `profiles1` and ",
      "`profiles2` must both be NULL.",
      call. = FALSE
    )
  }

  # TOP-LEVEL SCENARIO VALIDATION

  if (!is.list(profiles) ||
      length(profiles) < 2L) {
    stop(
      "`profiles` must be a named list containing at least two scenarios.",
      call. = FALSE
    )
  }

  scenario_names <- names(profiles)

  if (is.null(scenario_names) ||
      anyNA(scenario_names) ||
      any(!nzchar(scenario_names)) ||
      anyDuplicated(scenario_names)) {
    stop(
      "Every scenario in `profiles` must have a unique, non-empty name.",
      call. = FALSE
    )
  }

  # NORMALIZE EACH SCENARIO TO MATCHED PROFILE SETS

  looks_like_legacy_simulation <- function(x) {
    if (!is.list(x) ||
        is.null(names(x)) ||
        !"meta" %in% names(x)) {
      return(FALSE)
    }

    has_single <-
      "profile" %in% names(x) &&
      is.numeric(x$profile)

    has_multiple <-
      "profiles" %in% names(x) &&
      is.list(x$profiles) &&
      length(x$profiles) >= 1L

    has_single || has_multiple
  }

  is_simulated_exposure_object <- function(x) {
    inherits(
      x,
      "epiexposure_simulated_exposures"
    ) ||
      looks_like_legacy_simulation(x)
  }

  extract_simulated_profiles <- function(
    x,
    variable,
    scenario_name,
    expected_max_lag
  ) {

    new_contract <- inherits(
      x,
      "epiexposure_simulated_exposures"
    )

    if (new_contract) {
      if (!is.list(x$profiles) ||
          !length(x$profiles)) {
        stop(
          "The `simulate_exposures()` object used for exposure '",
          variable,
          "' in scenario '",
          scenario_name,
          "' is missing its canonical `$profiles` collection.",
          call. = FALSE
        )
      }

      meta <- x$meta

      if (!is.list(meta) ||
          !identical(
            meta$profile_order,
            "chronological"
          ) ||
          !identical(
            meta$profile_contract,
            "canonical_profiles_list_chronological"
          ) ||
          !identical(
            meta$simulation_contract,
            "n_exposure_profiles_no_model_parameter_uncertainty"
          )) {
        stop(
          "The `simulate_exposures()` object for exposure '",
          variable,
          "' in scenario '",
          scenario_name,
          "' does not satisfy the audited chronological exposure-simulation ",
          "contract.",
          call. = FALSE
        )
      }

      simulation_max_lag <-
        meta$max_lag

      profile_collection <-
        x$profiles

      source_contract <-
        "simulate_exposures_audited"

    } else {

      # Compatibility with the immediately preceding EpiExposure development
      # structure. This route is intentionally strict and relies on explicit
      # chronological metadata rather than guessing from vector names.
      meta <- x$meta

      if ("profile" %in% names(x) &&
          is.numeric(x$profile)) {

        if (!is.list(meta) ||
            !identical(
              meta$profile_order,
              "chronological"
            ) ||
            is.null(meta$max_lag)) {
          stop(
            "Legacy simulated exposure object for '",
            variable,
            "' in scenario '",
            scenario_name,
            "' lacks explicit chronological/max-lag metadata.",
            call. = FALSE
          )
        }

        profile_collection <-
          list(
            legacy_profile =
              as.numeric(
                x$profile
              )
          )

        simulation_max_lag <-
          meta$max_lag

      } else {

        if (!is.list(x$profiles) ||
            !length(x$profiles) ||
            !is.list(meta) ||
            !length(meta)) {
          stop(
            "Legacy simulated exposure object for '",
            variable,
            "' in scenario '",
            scenario_name,
            "' is malformed.",
            call. = FALSE
          )
        }

        meta_max_lag <- vapply(
          meta,
          function(m) {
            if (!is.list(m) ||
                !identical(
                  m$profile_order,
                  "chronological"
                ) ||
                is.null(m$max_lag)) {
              return(NA_real_)
            }
            as.numeric(m$max_lag)[1L]
          },
          numeric(1)
        )

        if (anyNA(meta_max_lag) ||
            length(
              unique(
                meta_max_lag
              )
            ) != 1L) {
          stop(
            "Legacy simulated profiles for exposure '",
            variable,
            "' in scenario '",
            scenario_name,
            "' do not share one explicit chronological/max-lag contract.",
            call. = FALSE
          )
        }

        profile_collection <-
          x$profiles

        simulation_max_lag <-
          meta_max_lag[1L]

      }

      source_contract <-
        "simulate_exposures_legacy"
    }

    if (!is.numeric(
      simulation_max_lag
    ) ||
    length(
      simulation_max_lag
    ) != 1L ||
    is.na(
      simulation_max_lag
    ) ||
    !is.finite(
      simulation_max_lag
    ) ||
    simulation_max_lag !=
    as.integer(
      simulation_max_lag
    )) {
      stop(
        "The simulated exposure object for '",
        variable,
        "' in scenario '",
        scenario_name,
        "' has invalid `max_lag` metadata.",
        call. = FALSE
      )
    }

    simulation_max_lag <-
      as.integer(
        simulation_max_lag
      )

    if (!identical(
      simulation_max_lag,
      as.integer(
        expected_max_lag
      )
    )) {
      stop(
        "The simulated exposure object for '",
        variable,
        "' in scenario '",
        scenario_name,
        "' was generated with max_lag = ",
        simulation_max_lag,
        ", but the fitted exposure requires max_lag = ",
        expected_max_lag,
        ". Regenerate the exposure profile using the fitted lag window.",
        call. = FALSE
      )
    }

    list(
      profiles =
        profile_collection,
      source =
        source_contract
    )
  }

  as_profile_list <- function(
    x,
    variable,
    scenario_name,
    expected_max_lag
  ) {

    required_length <-
      as.integer(
        expected_max_lag
      ) + 1L

    source <- NULL
    if (is_simulated_exposure_object(x)) {

      extracted <-
        extract_simulated_profiles(
          x = x,
          variable = variable,
          scenario_name = scenario_name,
          expected_max_lag =
            expected_max_lag
        )

      profiles_out <-
        extracted$profiles

      source <-
        extracted$source

    } else if (
      is.numeric(x) &&
      is.null(dim(x))
    ) {

      if (length(x) == 1L) {
        if (is.na(x) ||
            !is.finite(x)) {
          stop(
            "Constant exposure value for '",
            variable,
            "' in scenario '",
            scenario_name,
            "' must be finite.",
            call. = FALSE
          )
        }

        profiles_out <- list(
          rep(
            as.numeric(x),
            required_length
          )
        )

        source <-
          "constant_scalar"

      } else {

        profiles_out <- list(
          as.numeric(x)
        )

        source <-
          "manual_vector"
      }

    } else if (
      is.matrix(x) ||
      is.data.frame(x)
    ) {

      x_matrix <- as.matrix(x)

      if (!is.numeric(x_matrix) ||
          !nrow(x_matrix)) {
        stop(
          "Matrix/data.frame profiles for variable '",
          variable,
          "' in scenario '",
          scenario_name,
          "' must contain at least one numeric row.",
          call. = FALSE
        )
      }

      profiles_out <- lapply(
        seq_len(nrow(x_matrix)),
        function(i) {
          as.numeric(
            x_matrix[i, ]
          )
        }
      )

      source <-
        "manual_matrix"

    } else if (is.list(x)) {

      if (!length(x)) {
        stop(
          "No profiles were supplied for variable '",
          variable,
          "' in scenario '",
          scenario_name,
          "'.",
          call. = FALSE
        )
      }

      profiles_out <- lapply(
        seq_along(x),
        function(i) {
          current <- x[[i]]

          if (is.matrix(current) ||
              is.data.frame(current)) {

            current_matrix <-
              as.matrix(current)

            if (nrow(current_matrix) !=
                1L ||
                !is.numeric(
                  current_matrix
                )) {
              stop(
                "Nested profile ",
                i,
                " for variable '",
                variable,
                "' in scenario '",
                scenario_name,
                "' must contain exactly one numeric row.",
                call. = FALSE
              )
            }

            current <-
              as.numeric(
                current_matrix[1L, ]
              )
          }

          if (!is.numeric(current) ||
              is.list(current) ||
              !is.null(dim(current))) {
            stop(
              "Profile ",
              i,
              " for variable '",
              variable,
              "' in scenario '",
              scenario_name,
              "' must be a numeric vector.",
              call. = FALSE
            )
          }

          as.numeric(current)
        }
      )

      source <-
        "manual_list"

    } else {

      stop(
        "Exposure '",
        variable,
        "' in scenario '",
        scenario_name,
        "' must be supplied as a finite scalar, numeric profile, profile ",
        "matrix/list, or `simulate_exposures()` object.",
        call. = FALSE
      )
    }

    if (!is.list(
      profiles_out
    ) ||
    !length(
      profiles_out
    )) {
      stop(
        "No usable profiles were obtained for exposure '",
        variable,
        "' in scenario '",
        scenario_name,
        "'.",
        call. = FALSE
      )
    }

    for (profile_index in seq_along(
      profiles_out
    )) {

      current <-
        profiles_out[[
          profile_index
        ]]

      if (!is.numeric(current) ||
          !length(current) ||
          anyNA(current) ||
          any(!is.finite(current))) {
        stop(
          "Scenario '",
          scenario_name,
          "', exposure '",
          variable,
          "', profile ",
          profile_index,
          " must contain only finite numeric values.",
          call. = FALSE
        )
      }

      if (length(current) !=
          required_length) {
        stop(
          "Scenario '",
          scenario_name,
          "', exposure '",
          variable,
          "', profile ",
          profile_index,
          " has length ",
          length(current),
          ", but the fitted max_lag = ",
          expected_max_lag,
          " requires exactly ",
          required_length,
          " chronological values. Use a scalar for a constant profile or ",
          "supply/regenerate the complete fitted lag profile.",
          call. = FALSE
        )
      }

      profiles_out[[
        profile_index
      ]] <- as.numeric(
        current
      )
    }

    list(
      profiles =
        profiles_out,
      source =
        source,
      n_profiles =
        length(
          profiles_out
        )
    )
  }

  normalize_scenario <- function(
    scenario,
    scenario_name
  ) {

    # Single-exposure convenience: the entire scenario can itself be a scalar,
    # explicit profile collection, or simulate_exposures() object.
    if (length(fitted_vars) == 1L) {
      fitted_var <- fitted_vars[[1L]]

      if (is_simulated_exposure_object(
        scenario
      ) ||
      is.numeric(scenario) ||
      is.matrix(scenario) ||
      is.data.frame(scenario)) {

        scenario <-
          stats::setNames(
            list(scenario),
            fitted_var
          )

      } else if (is.list(scenario)) {

        scenario_names_local <-
          names(scenario)

        is_explicit_variable_map <-
          length(scenario) == 1L &&
          !is.null(
            scenario_names_local
          ) &&
          identical(
            scenario_names_local,
            fitted_var
          )

        if (!is_explicit_variable_map) {
          scenario <-
            stats::setNames(
              list(scenario),
              fitted_var
            )
        }
      }
    }

    if (!is.list(scenario) ||
        !length(scenario) ||
        is.null(names(scenario)) ||
        anyNA(names(scenario)) ||
        any(!nzchar(names(scenario))) ||
        anyDuplicated(names(scenario))) {
      stop(
        "Scenario '",
        scenario_name,
        "' must provide a named profile structure for the fitted exposure ",
        "variables.",
        call. = FALSE
      )
    }

    missing_vars <- setdiff(
      fitted_vars,
      names(scenario)
    )

    extra_vars <- setdiff(
      names(scenario),
      fitted_vars
    )

    if (length(missing_vars) ||
        length(extra_vars)) {

      message_parts <-
        character(0)

      if (length(missing_vars)) {
        message_parts <- c(
          message_parts,
          paste0(
            "missing: ",
            paste(
              missing_vars,
              collapse = ", "
            )
          )
        )
      }

      if (length(extra_vars)) {
        message_parts <- c(
          message_parts,
          paste0(
            "not fitted: ",
            paste(
              extra_vars,
              collapse = ", "
            )
          )
        )
      }

      stop(
        "Scenario '",
        scenario_name,
        "' exposure names must match the fitted model exactly (",
        paste(
          message_parts,
          collapse = "; "
        ),
        ").",
        call. = FALSE
      )
    }

    scenario <- scenario[
      fitted_vars
    ]

    normalized_variables <-
      stats::setNames(
        vector(
          "list",
          length(
            fitted_vars
          )
        ),
        fitted_vars
      )

    source_rows <-
      vector(
        "list",
        length(
          fitted_vars
        )
      )

    for (var_index in seq_along(
      fitted_vars
    )) {

      variable <-
        fitted_vars[[
          var_index
        ]]

      spec <- fitted_spec[[
        variable
      ]]

      if (is.null(spec) ||
          is.null(
            spec$max_lag
          )) {
        stop(
          "Missing fitted max-lag specification for exposure '",
          variable,
          "'.",
          call. = FALSE
        )
      }

      max_lag_use <-
        as.integer(
          max(
            spec$max_lag
          )
        )

      normalized <-
        as_profile_list(
          x =
            scenario[[
              variable
            ]],
          variable =
            variable,
          scenario_name =
            scenario_name,
          expected_max_lag =
            max_lag_use
        )

      normalized_variables[[
        variable
      ]] <-
        normalized

      source_rows[[
        var_index
      ]] <- data.frame(
        scenario =
          scenario_name,
        variable =
          variable,
        source =
          normalized$source,
        n_profiles_input =
          normalized$n_profiles,
        max_lag =
          max_lag_use,
        stringsAsFactors = FALSE
      )
    }

    profile_counts <- vapply(
      normalized_variables,
      function(x) {
        x$n_profiles
      },
      integer(1)
    )

    n_profiles <- max(
      profile_counts
    )

    incompatible <-
      profile_counts != 1L &
      profile_counts != n_profiles

    if (any(incompatible)) {
      stop(
        "Within scenario '",
        scenario_name,
        "', every fitted exposure must provide either one profile or the ",
        "common maximum number of profiles. Received: ",
        paste(
          paste0(
            names(
              profile_counts
            ),
            "=",
            profile_counts
          ),
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    profile_sets <-
      stats::setNames(
        vector(
          "list",
          length(
            fitted_vars
          )
        ),
        fitted_vars
      )

    for (variable in fitted_vars) {

      current_profiles <-
        normalized_variables[[
          variable
        ]]$profiles

      if (length(
        current_profiles
      ) == 1L &&
      n_profiles > 1L) {

        current_profiles <-
          rep(
            current_profiles,
            n_profiles
          )
      }

      profile_sets[[
        variable
      ]] <- current_profiles
    }

    list(
      sets =
        profile_sets,
      n_profiles =
        n_profiles,
      sources =
        do.call(
          rbind,
          source_rows
        )
    )
  }

  normalized_scenarios <- lapply(
    scenario_names,
    function(scenario_name) {
      normalize_scenario(
        profiles[[scenario_name]],
        scenario_name
      )
    }
  )
  names(normalized_scenarios) <-
    scenario_names

  scenario_profile_counts <- vapply(
    normalized_scenarios,
    function(x) x$n_profiles,
    integer(1)
  )

  profile_source_metadata <- do.call(
    rbind,
    lapply(
      normalized_scenarios,
      function(x) {
        x$sources
      }
    )
  )
  rownames(profile_source_metadata) <- NULL

  common_profile_count <- max(
    scenario_profile_counts
  )

  incompatible_scenarios <-
    scenario_profile_counts != 1L &
    scenario_profile_counts !=
    common_profile_count

  if (any(incompatible_scenarios)) {
    stop(
      "Top-level scenarios must contain either one matched profile set or the ",
      "common maximum number of matched profile sets. Received: ",
      paste(
        paste0(
          names(
            scenario_profile_counts
          ),
          "=",
          scenario_profile_counts
        ),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  # Recycle one-profile scenarios across the common profile index.
  if (common_profile_count > 1L) {
    for (scenario_name in scenario_names) {
      if (normalized_scenarios[[
        scenario_name
      ]]$n_profiles == 1L) {

        for (variable in fitted_vars) {
          normalized_scenarios[[
            scenario_name
          ]]$sets[[variable]] <-
            rep(
              normalized_scenarios[[
                scenario_name
              ]]$sets[[variable]],
              common_profile_count
            )
        }

        normalized_scenarios[[
          scenario_name
        ]]$n_profiles <-
          common_profile_count
      }
    }
  }

  # BUILD ONE JOINT PREDICTION BATCH FOR ALL SCENARIOS

  batch_profiles <- stats::setNames(
    vector(
      "list",
      length(fitted_vars)
    ),
    fitted_vars
  )

  for (variable in fitted_vars) {
    batch_profiles[[variable]] <-
      vector(
        "list",
        length(scenario_names) *
          common_profile_count
      )
  }

  batch_key_rows <- vector(
    "list",
    length(scenario_names) *
      common_profile_count
  )

  batch_index <- 1L

  for (scenario_name in scenario_names) {
    for (profile_index in seq_len(
      common_profile_count
    )) {

      for (variable in fitted_vars) {
        batch_profiles[[variable]][[
          batch_index
        ]] <-
          normalized_scenarios[[
            scenario_name
          ]]$sets[[variable]][[
            profile_index
          ]]
      }

      batch_key_rows[[batch_index]] <-
        data.frame(
          batch_profile = batch_index,
          scenario = scenario_name,
          profile = profile_index,
          stringsAsFactors = FALSE
        )

      batch_index <- batch_index + 1L
    }
  }

  batch_keys <- do.call(
    rbind,
    batch_key_rows
  )
  rownames(batch_keys) <- NULL

  n_batch_predictions <-
    nrow(batch_keys)

  # ONE AND ONLY ONE PREDICT_OUTCOMES() CALL

  backend_output <- if (
    uncertainty
  ) {
    "samples"
  } else {
    "summary"
  }

  backend <- predict_outcomes(
    fit = fit,
    profiles = batch_profiles,
    type = type,
    uncertainty = uncertainty,
    output = backend_output,
    n_samples = n_samples,
    probs = probs,
    seed = seed,
    extrapolation = extrapolation
  )

  if (!is.data.frame(backend) ||
      !nrow(backend)) {
    stop(
      "`predict_outcomes()` did not return a non-empty prediction data frame.",
      call. = FALSE
    )
  }

  required_backend_columns <- c(
    "profile",
    "prediction"
  )

  if (uncertainty) {
    required_backend_columns <-
      c(
        required_backend_columns,
        "sample"
      )
  }

  missing_backend <- setdiff(
    required_backend_columns,
    names(backend)
  )

  if (length(missing_backend)) {
    stop(
      "`predict_outcomes()` is missing required joint-batch output column(s): ",
      paste(
        missing_backend,
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  if (!is.numeric(backend$prediction) ||
      anyNA(backend$prediction) ||
      any(!is.finite(
        backend$prediction
      ))) {
    stop(
      "`predict_outcomes()` returned non-finite predictions.",
      call. = FALSE
    )
  }

  if (!is.numeric(backend$profile) ||
      anyNA(backend$profile) ||
      any(!is.finite(backend$profile)) ||
      any(
        backend$profile !=
        as.integer(
          backend$profile
        )
      )) {
    stop(
      "`predict_outcomes()` returned invalid batch profile indices.",
      call. = FALSE
    )
  }

  backend$profile <-
    as.integer(
      backend$profile
    )

  if (any(
    backend$profile < 1L |
    backend$profile >
    n_batch_predictions
  )) {
    stop(
      "`predict_outcomes()` returned batch profile indices outside the expected ",
      "range.",
      call. = FALSE
    )
  }

  # Strictly verify the backend contract rather than assuming that another
  # version of predict_outcomes() has equivalent semantics.
  if (!identical(
    attr(
      backend,
      "epiexposure_prediction_level",
      exact = TRUE
    ),
    "population"
  ) ||
  !identical(
    attr(
      backend,
      "epiexposure_prediction_estimand",
      exact = TRUE
    ),
    "expected_response"
  ) ||
  !identical(
    attr(
      backend,
      "epiexposure_prediction_type",
      exact = TRUE
    ),
    type
  ) ||
  !identical(
    attr(
      backend,
      "epiexposure_uncertainty",
      exact = TRUE
    ),
    uncertainty
  )) {
    stop(
      "`predict_outcomes()` returned metadata inconsistent with the current ",
      "EpiExposure comparison contract.",
      call. = FALSE
    )
  }

  # MAP JOINT BATCH ROWS BACK TO SCENARIO + WITHIN-SCENARIO PROFILE

  key_match <- match(
    backend$profile,
    batch_keys$batch_profile
  )

  if (anyNA(key_match)) {
    stop(
      "Could not map one or more backend predictions to their scenario profile.",
      call. = FALSE
    )
  }

  predictions_long <- data.frame(
    scenario =
      batch_keys$scenario[
        key_match
      ],
    profile =
      batch_keys$profile[
        key_match
      ],
    stringsAsFactors = FALSE
  )

  if (uncertainty) {
    if (!is.numeric(backend$sample) ||
        anyNA(backend$sample) ||
        any(!is.finite(
          backend$sample
        )) ||
        any(
          backend$sample !=
          as.integer(
            backend$sample
          )
        ) ||
        any(backend$sample < 1L)) {
      stop(
        "`predict_outcomes()` returned invalid parameter-draw sample IDs.",
        call. = FALSE
      )
    }

    predictions_long$sample <-
      as.integer(
        backend$sample
      )
  }

  predictions_long$prediction <-
    as.numeric(
      backend$prediction
    )

  # STRICT SUPPORT VALIDATION

  if (!uncertainty) {
    expected_rows <-
      length(scenario_names) *
      common_profile_count

    if (nrow(predictions_long) !=
        expected_rows) {
      stop(
        "The deterministic joint prediction batch returned ",
        nrow(predictions_long),
        " rows; expected ",
        expected_rows,
        ".",
        call. = FALSE
      )
    }

    if (anyDuplicated(
      predictions_long[
        ,
        c(
          "scenario",
          "profile"
        ),
        drop = FALSE
      ]
    )) {
      stop(
        "Deterministic joint predictions are not unique by scenario and ",
        "profile.",
        call. = FALSE
      )
    }

  } else {

    sample_ids <- sort(
      unique(
        predictions_long$sample
      )
    )

    if (!identical(
      sample_ids,
      seq_len(n_samples)
    )) {
      stop(
        "Joint prediction sample IDs are incomplete or do not equal ",
        "`1:n_samples`.",
        call. = FALSE
      )
    }

    if (anyDuplicated(
      predictions_long[
        ,
        c(
          "scenario",
          "profile",
          "sample"
        ),
        drop = FALSE
      ]
    )) {
      stop(
        "Joint prediction draws are not unique by scenario, profile, and sample.",
        call. = FALSE
      )
    }

    expected_rows <-
      length(scenario_names) *
      common_profile_count *
      n_samples

    if (nrow(predictions_long) !=
        expected_rows) {
      stop(
        "The joint uncertainty batch returned ",
        nrow(predictions_long),
        " rows; expected ",
        expected_rows,
        ".",
        call. = FALSE
      )
    }

    support_key <- interaction(
      predictions_long$scenario,
      predictions_long$profile,
      drop = TRUE,
      lex.order = TRUE
    )

    sample_sets <- split(
      predictions_long$sample,
      support_key
    )

    incomplete_support <- vapply(
      sample_sets,
      function(x) {
        !identical(
          sort(unique(x)),
          seq_len(n_samples)
        )
      },
      logical(1)
    )

    if (any(incomplete_support)) {
      stop(
        "At least one scenario/profile combination is missing one or more joint ",
        "parameter draws.",
        call. = FALSE
      )
    }
  }

  # PAIRWISE COMPARISON HELPER

  add_comparison_metrics <- function(
    comparison_data
  ) {

    comparison_data$diff <-
      comparison_data$pred2 -
      comparison_data$pred1

    comparison_data$abs_diff <-
      abs(
        comparison_data$diff
      )

    # Relative outcome changes are intentionally response-scale only. A
    # positive denominator is required because the returned percent change is a
    # conventional relative change, not merely an arbitrary algebraic quotient.
    comparison_data$relative_change_defined <-
      identical(
        type,
        "response"
      ) &
      comparison_data$pred1 >
      eps

    comparison_data$ratio <-
      NA_real_

    comparison_data$percent_change <-
      NA_real_

    relative_rows <-
      comparison_data$relative_change_defined

    if (any(relative_rows)) {
      comparison_data$ratio[
        relative_rows
      ] <-
        comparison_data$pred2[
          relative_rows
        ] /
        comparison_data$pred1[
          relative_rows
        ]

      if (any(
        !is.finite(
          comparison_data$ratio[
            relative_rows
          ]
        )
      )) {
        stop(
          "A finite response-scale scenario ratio could not be calculated.",
          call. = FALSE
        )
      }

      comparison_data$percent_change[
        relative_rows
      ] <-
        100 * (
          comparison_data$ratio[
            relative_rows
          ] - 1
        )
    }

    probability_family <-
      family_name %in%
      c(
        "beta",
        "binomial"
      )

    comparison_data$percentage_point_change <-
      if (
        probability_family &&
        identical(
          type,
          "response"
        )
      ) {
        100 *
          comparison_data$diff
      } else {
        rep(
          NA_real_,
          nrow(
            comparison_data
          )
        )
      }

    comparison_data
  }

  # PAIR ALL SCENARIOS USING IDENTICAL PROFILE / SAMPLE SUPPORT

  combinations <- utils::combn(
    scenario_names,
    2,
    simplify = FALSE
  )

  comparison_blocks <- lapply(
    combinations,
    function(comparison) {

      scenario1 <- comparison[1L]
      scenario2 <- comparison[2L]

      first <- predictions_long[
        predictions_long$scenario ==
          scenario1,
        ,
        drop = FALSE
      ]

      second <- predictions_long[
        predictions_long$scenario ==
          scenario2,
        ,
        drop = FALSE
      ]

      matching_columns <- "profile"

      if (uncertainty) {
        matching_columns <- c(
          matching_columns,
          "sample"
        )
      }

      if (anyDuplicated(
        first[
          ,
          matching_columns,
          drop = FALSE
        ]
      ) ||
      anyDuplicated(
        second[
          ,
          matching_columns,
          drop = FALSE
        ]
      )) {
        stop(
          "Prediction support is not unique within scenarios '",
          scenario1,
          "' and '",
          scenario2,
          "'.",
          call. = FALSE
        )
      }

      first_key <- do.call(
        paste,
        c(
          first[
            ,
            matching_columns,
            drop = FALSE
          ],
          sep = "\r"
        )
      )

      second_key <- do.call(
        paste,
        c(
          second[
            ,
            matching_columns,
            drop = FALSE
          ],
          sep = "\r"
        )
      )

      if (!setequal(
        first_key,
        second_key
      )) {
        stop(
          "Scenarios '",
          scenario1,
          "' and '",
          scenario2,
          "' do not share identical matched prediction support.",
          call. = FALSE
        )
      }

      first_keep <- first[
        ,
        c(
          matching_columns,
          "prediction"
        ),
        drop = FALSE
      ]

      second_keep <- second[
        ,
        c(
          matching_columns,
          "prediction"
        ),
        drop = FALSE
      ]

      names(first_keep)[
        names(first_keep) ==
          "prediction"
      ] <- "pred1"

      names(second_keep)[
        names(second_keep) ==
          "prediction"
      ] <- "pred2"

      comparison_data <- merge(
        first_keep,
        second_keep,
        by = matching_columns,
        all = TRUE,
        sort = FALSE
      )

      if (nrow(comparison_data) !=
          nrow(first_keep) ||
          nrow(comparison_data) !=
          nrow(second_keep) ||
          anyNA(
            comparison_data$pred1
          ) ||
          anyNA(
            comparison_data$pred2
          )) {
        stop(
          "Predictions could not be matched one-to-one between scenarios '",
          scenario1,
          "' and '",
          scenario2,
          "'.",
          call. = FALSE
        )
      }

      comparison_data$scenario1 <-
        scenario1
      comparison_data$scenario2 <-
        scenario2

      comparison_data <-
        add_comparison_metrics(
          comparison_data
        )

      output_columns <- c(
        "scenario1",
        "scenario2",
        "profile",
        if (uncertainty) "sample",
        "pred1",
        "pred2",
        "diff",
        "abs_diff",
        "ratio",
        "percent_change",
        "relative_change_defined",
        "percentage_point_change"
      )

      comparison_data[
        ,
        output_columns,
        drop = FALSE
      ]
    }
  )

  comparison_samples <- do.call(
    rbind,
    comparison_blocks
  )
  rownames(comparison_samples) <-
    NULL

  comparison_samples <- comparison_samples[
    order(
      match(
        comparison_samples$scenario1,
        scenario_names
      ),
      match(
        comparison_samples$scenario2,
        scenario_names
      ),
      comparison_samples$profile,
      if (uncertainty) {
        comparison_samples$sample
      } else {
        rep(
          1L,
          nrow(
            comparison_samples
          )
        )
      }
    ),
    ,
    drop = FALSE
  ]
  rownames(comparison_samples) <-
    NULL

  multiple_profile_sets <-
    common_profile_count > 1L

  # OUTPUT ATTRIBUTE HELPER

  attach_comparison_attributes <- function(
    out,
    output_kind
  ) {

    attr(
      out,
      "epiexposure_family_name"
    ) <- family_name

    attr(
      out,
      "epiexposure_link"
    ) <- link_name

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
      "epiexposure_prediction_type"
    ) <- type

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
      "epiexposure_uncertainty"
    ) <- uncertainty

    attr(
      out,
      "epiexposure_comparison_output"
    ) <- output_kind

    attr(
      out,
      "epiexposure_scenario_order"
    ) <- scenario_names

    attr(
      out,
      "epiexposure_scenario_profile_counts_input"
    ) <- scenario_profile_counts

    attr(
      out,
      "epiexposure_profile_sources"
    ) <- profile_source_metadata

    attr(
      out,
      "epiexposure_n_matched_profiles"
    ) <- common_profile_count

    attr(
      out,
      "epiexposure_profile_input_contract"
    ) <- "scalar_manual_or_simulated_exposures_matched_by_index"

    attr(
      out,
      "epiexposure_comparison_direction"
    ) <- "scenario2_minus_scenario1"

    attr(
      out,
      "epiexposure_relative_change_contract"
    ) <- if (
      identical(
        type,
        "response"
      )
    ) {
      "pred2_over_pred1_when_pred1_gt_eps"
    } else {
      "not_defined_on_link_scale"
    }

    attr(
      out,
      "epiexposure_probability_difference_contract"
    ) <- if (
      family_name %in%
      c(
        "beta",
        "binomial"
      ) &&
      identical(
        type,
        "response"
      )
    ) {
      "100_times_pred2_minus_pred1_percentage_points"
    } else {
      "not_applicable"
    }

    attr(
      out,
      "epiexposure_joint_scenario_draws"
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
          "epiexposure_probs"
        ) <- probs

        attr(
          out,
          "epiexposure_summary_center"
        ) <- "median"
      }
    }

    out
  }

  # DETERMINISTIC RETURN

  if (!uncertainty) {
    deterministic_out <-
      comparison_samples

    if (!multiple_profile_sets) {
      deterministic_out$profile <-
        NULL
    }

    deterministic_out <-
      attach_comparison_attributes(
        deterministic_out,
        output_kind = "summary"
      )

    return(
      deterministic_out
    )
  }

  # SAMPLE-LEVEL RETURN

  if (identical(
    output,
    "samples"
  )) {
    samples_out <-
      comparison_samples

    if (!multiple_profile_sets) {
      samples_out$profile <-
        NULL
    }

    samples_out <-
      attach_comparison_attributes(
        samples_out,
        output_kind = "samples"
      )

    return(
      samples_out
    )
  }

  # UNCERTAINTY SUMMARY AFTER DRAW-BY-DRAW COMPARISON

  grouping_columns <- c(
    "scenario1",
    "scenario2",
    "profile"
  )

  group_key <- interaction(
    comparison_samples$scenario1,
    comparison_samples$scenario2,
    comparison_samples$profile,
    drop = TRUE,
    lex.order = TRUE
  )

  row_groups <- split(
    seq_len(
      nrow(
        comparison_samples
      )
    ),
    group_key
  )

  summary_rows <- vector(
    "list",
    length(row_groups)
  )

  summary_index <- 1L

  base_quantities <- c(
    "pred1",
    "pred2",
    "diff",
    "abs_diff"
  )

  probability_difference_applicable <-
    family_name %in%
    c(
      "beta",
      "binomial"
    ) &&
    identical(
      type,
      "response"
    )

  for (idx in row_groups) {
    current <- comparison_samples[
      idx,
      ,
      drop = FALSE
    ]

    if (!identical(
      sort(
        unique(
          current$sample
        )
      ),
      seq_len(
        n_samples
      )
    )) {
      stop(
        "A pairwise comparison group is missing one or more joint parameter ",
        "draws.",
        call. = FALSE
      )
    }

    base <- current[
      1L,
      grouping_columns,
      drop = FALSE
    ]

    for (quantity in base_quantities) {
      values <- summarise_numeric(
        current[[quantity]],
        probs_use = probs
      )

      base[[quantity]] <-
        unname(
          values["estimate"]
        )
      base[[paste0(
        quantity,
        "_sd"
      )]] <-
        unname(
          values["sd"]
        )
      base[[paste0(
        quantity,
        "_lower"
      )]] <-
        unname(
          values["lower"]
        )
      base[[paste0(
        quantity,
        "_upper"
      )]] <-
        unname(
          values["upper"]
        )
    }

    if (probability_difference_applicable) {
      values <- summarise_numeric(
        current$percentage_point_change,
        probs_use = probs
      )

      base$percentage_point_change <-
        unname(
          values["estimate"]
        )
      base$percentage_point_change_sd <-
        unname(
          values["sd"]
        )
      base$percentage_point_change_lower <-
        unname(
          values["lower"]
        )
      base$percentage_point_change_upper <-
        unname(
          values["upper"]
        )

    } else {
      base$percentage_point_change <-
        NA_real_
      base$percentage_point_change_sd <-
        NA_real_
      base$percentage_point_change_lower <-
        NA_real_
      base$percentage_point_change_upper <-
        NA_real_
    }

    relative_fraction <- mean(
      current$relative_change_defined
    )

    base$relative_defined_fraction <-
      relative_fraction

    base$relative_change_defined <-
      isTRUE(
        all(
          current$relative_change_defined
        )
      )

    if (base$relative_change_defined) {
      ratio_values <- summarise_numeric(
        current$ratio,
        probs_use = probs
      )

      percent_values <- summarise_numeric(
        current$percent_change,
        probs_use = probs
      )

      base$ratio <-
        unname(
          ratio_values["estimate"]
        )
      base$ratio_sd <-
        unname(
          ratio_values["sd"]
        )
      base$ratio_lower <-
        unname(
          ratio_values["lower"]
        )
      base$ratio_upper <-
        unname(
          ratio_values["upper"]
        )

      base$percent_change <-
        unname(
          percent_values["estimate"]
        )
      base$percent_change_sd <-
        unname(
          percent_values["sd"]
        )
      base$percent_change_lower <-
        unname(
          percent_values["lower"]
        )
      base$percent_change_upper <-
        unname(
          percent_values["upper"]
        )

    } else {
      base$ratio <- NA_real_
      base$ratio_sd <- NA_real_
      base$ratio_lower <- NA_real_
      base$ratio_upper <- NA_real_

      base$percent_change <- NA_real_
      base$percent_change_sd <- NA_real_
      base$percent_change_lower <- NA_real_
      base$percent_change_upper <- NA_real_
    }

    summary_rows[[summary_index]] <-
      base
    summary_index <- summary_index + 1L
  }

  summary_out <- do.call(
    rbind,
    summary_rows
  )
  rownames(summary_out) <- NULL

  summary_out <- summary_out[
    order(
      match(
        summary_out$scenario1,
        scenario_names
      ),
      match(
        summary_out$scenario2,
        scenario_names
      ),
      summary_out$profile
    ),
    ,
    drop = FALSE
  ]
  rownames(summary_out) <- NULL

  if (!multiple_profile_sets) {
    summary_out$profile <- NULL
  }

  attach_comparison_attributes(
    summary_out,
    output_kind = "summary"
  )
}
