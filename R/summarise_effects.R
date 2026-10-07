#' Summarise DLNM exposure-lag effects on link and response scales
#'
#' Summarises lag-specific or period-specific distributed lag nonlinear model
#' (DLNM) associations for models fitted with `fit_epidlnm()`. The function
#' preserves the centering logic of `dlnm::crosspred()`: each exposure-lag
#' association is expressed relative to a reference exposure value (`cen`) for
#' the focal exposure. EpiExposure then extends that DLNM contrast by anchoring
#' it to a joint reference exposure profile, allowing the same association to
#' be reported as `baseline`, `predicted`, and `delta` on the response scale.
#'
#' Exposure histories are supplied chronologically: the oldest exposure is
#' first and the most recent exposure is last. Internally, EpiExposure converts
#' these histories to the retrospective lag scale used by DLNMs, where lag 0
#' represents the most recent exposure observation and increasing lag values
#' represent progressively older exposure observations. All lag-specific and
#' period-specific summaries reported by `summarise_effects()` are therefore
#' indexed on this retrospective lag scale.
#'
#' @param fit Fitted model returned by the current `fit_epidlnm()`.
#' @param data Long-format exposure data containing `epi_id`, `time`, and every
#'   exposure variable used in the fitted model. `data` is used to define
#'   default exposure grids and reference values and to validate temporal
#'   coverage. Under the EpiExposure exact-history contract, every `epi_id`
#'   must contain exactly the common fitted `max_lag + 1` time points, and all
#'   fitted exposures therefore share the same temporal-history length. Longer
#'   histories are not truncated and shorter histories are not padded. `data`
#'   is **not** used to re-estimate spline knots, boundary knots, lag bases, or
#'   any other fitted cross-basis component.
#' @param vars Exposure-variable name or names to summarise, or `NULL` to
#'   summarise all fitted exposures.
#' @param scale Character. `"lag"` returns lag-specific associations.
#'   `"period"` aggregates lag-specific linear-predictor contrasts over
#'   intervals.
#' @param lag_periods Optional data frame containing `period`, `lag_start`, and
#'   `lag_end`. Required when `scale = "period"` and `incremental = FALSE`.
#'   Both `lag_start` and `lag_end` are interpreted on the retrospective lag
#'   scale, where lag 0 is the most recent exposure and increasing lag values
#'   correspond to older exposure observations. Period bounds are inclusive.
#'   Periods may be non-overlapping or overlapping; each period is interpreted
#'   independently.
#' @param probs Unique probabilities used only to construct the default exposure
#'   evaluation grid when `at` is not supplied for a variable. Defaults to
#'   `seq(0.05, 0.95, by = 0.01)`.
#' @param at Optional exposure values at which associations are evaluated.
#'   `NULL` uses the quantiles specified by `probs`. A numeric vector can be used
#'   when exactly one exposure is requested. For multiple exposures, supply a
#'   named list, for example
#'   `list(tmean = seq(20, 35, by = 0.25), rain = seq(0, 40, by = 1))`.
#'   A named list may be partial; omitted requested variables use `probs`.
#' @param ref Reference exposure specification.
#'
#'   Two forms are supported:
#'
#'   \itemize{
#'     \item Method-based: `list(method = "median", value = NULL)`,
#'       `list(method = "percentile", value = p)`, or
#'       `list(method = "fixed", value = x)`.
#'     \item Exposure-specific: a named list with exactly one finite reference
#'       value for every fitted exposure, for example
#'       `list(tmean = 25, rain = 5, wetness = 10)`.
#'   }
#'
#'   Method-based references are applied consistently to **every** fitted
#'   exposure. Thus `method = "median"` uses each exposure's own median;
#'   `method = "percentile"` uses the same percentile for each exposure; and
#'   `method = "fixed"` applies the same numeric value to every exposure.
#'   Because exposures often have different units, the exposure-specific form
#'   is recommended for multivariable models when fixed reference values are
#'   desired.
#'
#'   The focal exposure's reference value is passed to the DLNM centering
#'   operation (`cen`). The complete set of reference values simultaneously
#'   defines the joint reference profile used for `baseline`.
#' @param effect_measure Character. `"linear"` returns the DLNM contrast on the
#'   linear-predictor scale. `"exponentiated"` returns `exp(eta)` and `"percent"`
#'   returns `100 * (exp(eta) - 1)`.
#'
#'   Exponentiation is allowed only for fitted `"log"` and `"logit"` links,
#'   matching the convention used by `dlnm::crosspred()`. With a log link,
#'   `exp(eta)` is a response-scale ratio and `"percent"` is the corresponding
#'   percent relative change in the expected response. With a logit link,
#'   `exp(eta)` is an odds ratio and `"percent"` is the percent change in odds;
#'   neither is a percentage-point change in probability. For links such as
#'   `"identity"`, `"probit"`, `"cloglog"`, or `"inverse"`, use
#'   `effect_measure = "linear"`.
#' @param incremental Logical. Used only with `scale = "period"`. If `TRUE`,
#'   return cumulative effects from lag 0 through each successive lag. If
#'   `FALSE`, aggregate over `lag_periods`.
#' @param uncertainty Logical. If `FALSE`, use the harmonized central
#'   fixed/population parameter estimate. If `TRUE`, propagate joint
#'   fixed/population parameter uncertainty draw by draw.
#' @param output Character. `"summary"` returns deterministic values when
#'   `uncertainty = FALSE`, or median, SD, and empirical intervals when
#'   `uncertainty = TRUE`. `"samples"` returns one row per parameter draw and
#'   exposure-lag/period combination and therefore requires
#'   `uncertainty = TRUE`.
#' @param n_samples Positive integer number of parameter draws used when
#'   `uncertainty = TRUE`. At least two are required.
#' @param interval_probs Numeric vector of length two defining the empirical
#'   uncertainty interval. The default `c(0.025, 0.975)` gives a 95 percent
#'   interval.
#' @param diagnostics Logical. If `TRUE`, additionally return lag-ranking and
#'   lag-contribution diagnostics. Diagnostics are available only for
#'   `scale = "lag"` and are always based on the additive linear-predictor
#'   contrast `eta`, regardless of `effect_measure`. This avoids treating the
#'   neutral value 1 of an exponentiated effect as if it were an additive
#'   contribution.
#' @param seed Optional finite integer used for parameter sampling. The caller's
#'   global random-number state is restored when the function exits.
#' @param extrapolation Character controlling values outside the exposure range
#'   stored with the fitted cross-basis: `"error"` (default), `"warn"`, or
#'   `"allow"`. This does not re-estimate the basis; it only controls whether
#'   extrapolation is rejected, warned about, or allowed.
#'
#' @return A data frame, or when `diagnostics = TRUE`, an object of class
#'   `"epiexposure_effects"` containing `effects` and `diagnostics`.
#'
#'   The principal columns have the following meanings:
#'
#'   \describe{
#'     \item{`eta`}{DLNM association contrast on the fitted link/linear-predictor
#'       scale, centered at the focal exposure reference. For
#'       `scale = "lag"`, this is the contribution of the specified exposure
#'       value at that lag relative to the same lag at the reference exposure.
#'       For period summaries, lag-specific `eta` contrasts are summed on the
#'       additive linear-predictor scale.}
#'     \item{`effect`}{`eta` itself for `effect_measure = "linear"`, or its
#'       permitted log/logit transformation for `"exponentiated"` or
#'       `"percent"`.}
#'     \item{`baseline`}{Population-level expected outcome under the **joint
#'       reference exposure profile**: every fitted exposure is held at its own
#'       reference value at every fitted lag and fitted random effects are set
#'       to zero. This is an EpiExposure response-scale reference prediction; it
#'       is not the DLNM contrast itself and should not be confused with the
#'       neutral DLNM values `eta = 0` or `exp(eta) = 1`.}
#'     \item{`predicted`}{Population-level expected outcome after applying the
#'       focal DLNM contrast to that joint baseline. For a lag-specific row,
#'       only the focal exposure at that lag is changed from its reference to
#'       `value`; all other lags and all other fitted exposures remain at their
#'       reference values. For a period row, the focal exposure is changed to
#'       `value` throughout the specified lag interval while all remaining
#'       exposure-lag positions remain at reference.}
#'     \item{`delta`}{Absolute response-scale difference
#'       `predicted - baseline`. For Gaussian/Gamma-type mean responses this is
#'       a difference in expected means; for Poisson/negative-binomial outcomes
#'       it is a difference in expected counts; and for Beta/Binomial outcomes
#'       it is a difference in expected proportions/probabilities. Multiplying
#'       `delta` by 100 for Beta or Binomial models gives percentage-point
#'       differences.}
#'   }
#'
#' @details
#' ## DLNM contrast and the EpiExposure response-scale extension
#'
#' `eta` is constructed using the centering mechanism of `dlnm::crosspred()`.
#' The fitted cross-basis is never redefined from `data`: EpiExposure uses the
#' stored `argvar`, `arglag`, maximum lag, and, when available, the original
#' stored `crossbasis` object. Consequently, `eta` retains the standard DLNM
#' interpretation as an association relative to `cen`.
#'
#' EpiExposure additionally defines a joint response-scale reference prediction
#'
#' \deqn{B = g^{-1}(X_{ref}\beta),}
#'
#' where every fitted exposure history is held at its selected reference value.
#' A lag- or period-specific DLNM contrast \eqn{\Delta\eta_{x,l}} is then
#' translated to the response scale as
#'
#' \deqn{P_{x,l} = g^{-1}\{g(B) + \Delta\eta_{x,l}\},}
#'
#' and
#'
#' \deqn{D_{x,l} = P_{x,l} - B.}
#'
#' Thus `baseline`, `predicted`, and `delta` do not replace the DLNM contrast;
#' they provide an additional outcome-scale interpretation anchored to the same
#' reference condition.
#'
#' ## Uncertainty
#'
#' With `uncertainty = TRUE`, one **joint** fixed/population coefficient draw is
#' used consistently across the complete requested effect surface. For draw
#' \eqn{s},
#'
#' \deqn{B^{(s)} = g^{-1}(X_{ref}\beta^{(s)}),}
#'
#' \deqn{P^{(s)}_{x,l} =
#'   g^{-1}\{X_{ref}\beta^{(s)} + \Delta\eta^{(s)}_{x,l}\},}
#'
#' and
#'
#' \deqn{D^{(s)}_{x,l} = P^{(s)}_{x,l} - B^{(s)}.}
#'
#' The same draw therefore determines `baseline`, `eta`, `effect`, `predicted`,
#' and `delta`, preserving their covariance. `baseline` changes across parameter
#' draws because the fitted parameters are uncertain, but for a given draw it
#' remains the same across all exposure-lag/period rows that use the same joint
#' reference profile. In summary output, the same baseline uncertainty summary
#' is therefore repeated across those rows.
#'
#' Parameter uncertainty includes the joint uncertainty of the fitted
#' fixed/population coefficients. Fitted group-specific random effects are
#' excluded because EpiExposure v1 reports population-level effects. Residual,
#' observation, process, dispersion, and posterior-predictive noise are not
#' added. The uncertainty interval therefore describes uncertainty in the
#' expected response and in the exposure-lag association, not the dispersion of
#' a future individual observation.
#'
#' For uncertainty summaries, all response-scale transformations are performed
#' draw by draw before medians, SDs, and empirical quantiles are calculated.
#'
#' ## Exact common lag/history contract
#'
#' `summarise_effects()` requires the fitted model to use one common
#' `max_lag` across all fitted exposure variables. If that common maximum lag
#' is `L`, every epidemic history supplied in `data` must contain exactly
#'
#' \deqn{
#'   L + 1
#' }
#'
#' equally spaced observations. Histories with fewer or more observations are
#' rejected explicitly. Because all exposure variables are columns of the same
#' validated long-format history and missing/non-finite exposure values are not
#' permitted, all fitted exposures necessarily use the same time support.
#'
#' ## Period effects
#'
#' Period definitions always use the retrospective lag scale employed by
#' EpiExposure and DLNMs. Thus, lag 0 corresponds to the most recent exposure
#' observation and increasing lag values correspond to progressively older
#' exposure observations. For example, if `max_lag = 85`, a period defined as
#' `0--14` represents the 15 most recent exposure observations, whereas a
#' period defined as `71--85` represents the oldest portion of the fitted
#' exposure history.
#'
#' DLNM contributions are additive on the linear-predictor scale. Period
#' summaries therefore sum lag-specific `eta` contrasts first and only then
#' transform the resulting contrast, if requested. With `incremental = TRUE`,
#' periods are `0-0`, `0-1`, ..., `0-L`. Otherwise each row of `lag_periods`
#' defines an independent lag interval.
#'
#' @export
summarise_effects <- function(
    fit,
    data,
    vars = NULL,
    scale = c("lag", "period"),
    lag_periods = NULL,
    probs = seq(0.05, 0.95, by = 0.01),
    at = NULL,
    ref = list(method = "median", value = NULL),
    effect_measure = c("linear", "exponentiated", "percent"),
    incremental = FALSE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    interval_probs = c(0.025, 0.975),
    diagnostics = FALSE,
    seed = NULL,
    extrapolation = c("error", "warn", "allow")
) {

  # ARGUMENTS AND STRICT FIT CONTRACT

  scale <- match.arg(scale)
  effect_measure <- match.arg(effect_measure)
  output <- match.arg(output)
  extrapolation <- match.arg(extrapolation)

  if (is.null(fit)) {
    stop("`fit` cannot be NULL.", call. = FALSE)
  }

  metadata <- .get_epiexposure_metadata(fit)

  if (!is.data.frame(data) || !nrow(data)) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }

  required_data <- unique(c("epi_id", "time", metadata$vars))
  missing_data <- setdiff(required_data, names(data))
  if (length(missing_data)) {
    stop(
      "`data` is missing required column(s): ",
      paste(missing_data, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (anyNA(data$epi_id)) {
    stop("`data$epi_id` cannot contain missing values.", call. = FALSE)
  }

  if (!is.numeric(data$time) || anyNA(data$time) ||
      any(!is.finite(data$time))) {
    stop("`data$time` must contain only finite numeric values.", call. = FALSE)
  }

  for (variable in metadata$vars) {
    if (!is.numeric(data[[variable]]) || anyNA(data[[variable]]) ||
        any(!is.finite(data[[variable]]))) {
      stop(
        "Exposure variable '", variable,
        "' must contain only finite numeric values.",
        call. = FALSE
      )
    }
  }

  # The same temporal assumptions used by the prediction basis layer apply
  # here. This validates spacing without re-estimating any basis component.
  .epix_validate_regular_time(
    data = data,
    group_index = as.character(data$epi_id),
    time_col = "time"
  )

  if (!identical(
    metadata$history_contract,
    "all_fitted_exposures_same_exact_max_lag_plus_one"
  )) {
    stop(
      "`summarise_effects()` requires the EpiExposure exact-history contract.",
      call. = FALSE
    )
  }

  required_history_length <- metadata$history_length

  if (!is.numeric(required_history_length) ||
      length(required_history_length) != 1L ||
      is.na(required_history_length) ||
      !is.finite(required_history_length) ||
      required_history_length < 1L ||
      required_history_length != as.integer(required_history_length) ||
      !identical(
        as.integer(required_history_length),
        as.integer(metadata$max_lag + 1L)
      )) {
    stop(
      "Stored EpiExposure history-length metadata are inconsistent with ",
      "`max_lag + 1`.",
      call. = FALSE
    )
  }

  required_history_length <- as.integer(required_history_length)

  group_counts <- table(as.character(data$epi_id))
  invalid_groups <- names(group_counts)[
    group_counts != required_history_length
  ]

  if (length(invalid_groups)) {
    examples <- paste0(
      invalid_groups,
      "=",
      as.integer(group_counts[invalid_groups])
    )

    stop(
      "Every epidemic in `data` must contain exactly ",
      required_history_length,
      " observations for fitted max_lag = ",
      metadata$max_lag,
      ". Non-matching epi_id(s) include: ",
      paste(utils::head(examples, 5L), collapse = ", "),
      if (length(examples) > 5L) "; ..." else ".",
      " Histories are not truncated, padded, or silently realigned.",
      call. = FALSE
    )
  }

  if (!is.logical(incremental) || length(incremental) != 1L ||
      is.na(incremental)) {
    stop("`incremental` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(uncertainty) || length(uncertainty) != 1L ||
      is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(diagnostics) || length(diagnostics) != 1L ||
      is.na(diagnostics)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }

  if (diagnostics && !identical(scale, "lag")) {
    stop(
      "`diagnostics = TRUE` is available only when `scale = 'lag'`.",
      call. = FALSE
    )
  }

  if (!uncertainty && identical(output, "samples")) {
    stop(
      "`output = 'samples'` requires `uncertainty = TRUE`.",
      call. = FALSE
    )
  }

  if (uncertainty) {
    n_samples <- .epix_validate_n_samples(n_samples)
    if (n_samples < 2L) {
      stop(
        "`n_samples` must be at least 2 when `uncertainty = TRUE`.",
        call. = FALSE
      )
    }
  }

  interval_probs <- .epix_validate_probs(interval_probs)

  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs < 0 | probs > 1)) {
    stop(
      "`probs` must contain finite probabilities between 0 and 1.",
      call. = FALSE
    )
  }
  if (anyDuplicated(probs)) {
    stop("`probs` must contain unique probabilities.", call. = FALSE)
  }
  probs <- sort(as.numeric(probs))

  if (!is.null(at) && !is.numeric(at) && !is.list(at)) {
    stop(
      "`at` must be NULL, a numeric vector, or a named list.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed != as.integer(seed)) {
      stop("`seed` must be NULL or one finite integer.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  if (effect_measure %in% c("exponentiated", "percent") &&
      !metadata$link %in% c("log", "logit")) {
    stop(
      "`effect_measure = '", effect_measure, "'` is not defined by ",
      "`summarise_effects()` for the fitted link '", metadata$link, "'. ",
      "Following the dlnm convention, exponentiated effects are supported ",
      "only for log and logit links. Use `effect_measure = 'linear'` for ",
      "identity, probit, cloglog, inverse, or other links.",
      call. = FALSE
    )
  }

  # VARIABLE SELECTION

  if (is.null(vars)) {
    variables <- metadata$vars
  } else {
    if (!is.character(vars) || !length(vars) || anyNA(vars) ||
        any(!nzchar(vars)) || anyDuplicated(vars)) {
      stop(
        "`vars` must be NULL or a character vector of unique non-empty ",
        "fitted exposure names.",
        call. = FALSE
      )
    }
    unknown <- setdiff(vars, metadata$vars)
    if (length(unknown)) {
      stop(
        "Variable(s) not found in the fitted model: ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    variables <- vars
  }

  # PERIOD VALIDATION

  if (identical(scale, "period") && !incremental && is.null(lag_periods)) {
    stop(
      "`lag_periods` is required when `scale = 'period'` and ",
      "`incremental = FALSE`.",
      call. = FALSE
    )
  }

  if (!is.null(lag_periods)) {
    required_period_columns <- c("period", "lag_start", "lag_end")

    if (!is.data.frame(lag_periods) ||
        !all(required_period_columns %in% names(lag_periods))) {
      stop(
        "`lag_periods` must be a data.frame containing `period`, ",
        "`lag_start`, and `lag_end`.",
        call. = FALSE
      )
    }

    if (anyNA(lag_periods$period) ||
        any(!nzchar(as.character(lag_periods$period))) ||
        anyDuplicated(as.character(lag_periods$period))) {
      stop(
        "`lag_periods$period` must contain unique non-empty labels.",
        call. = FALSE
      )
    }

    lag_periods$period <- as.character(lag_periods$period)

    for (nm in c("lag_start", "lag_end")) {
      x <- lag_periods[[nm]]
      if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
          any(x < 0) || any(x != as.integer(x))) {
        stop(
          "`lag_periods$", nm,
          "` must contain non-negative finite integers.",
          call. = FALSE
        )
      }
      lag_periods[[nm]] <- as.integer(x)
    }

    if (any(lag_periods$lag_start > lag_periods$lag_end)) {
      stop(
        "Every lag period must satisfy `lag_start <= lag_end`.",
        call. = FALSE
      )
    }
  }

  # RNG: LOCAL AND REPRODUCIBLE

  if (!is.null(seed)) {
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
          assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
        } else if (exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    set.seed(seed)
  }

  # SMALL LOCAL HELPERS

  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::sd(x)
  }

  safe_quantile <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(NA_real_, NA_real_))
    as.numeric(stats::quantile(
      x,
      probs = interval_probs,
      na.rm = TRUE,
      names = FALSE
    ))
  }

  transform_effect <- function(eta) {
    if (identical(effect_measure, "linear")) {
      return(eta)
    }
    if (identical(effect_measure, "exponentiated")) {
      return(exp(eta))
    }
    100 * (exp(eta) - 1)
  }

  handle_extrapolation <- function(values, variable, label) {
    definition <- .epix_basis_definition(metadata, variable)
    fitted_range <- definition$exposure_range

    if (is.null(fitted_range)) {
      fitted_range <- range(data[[variable]], na.rm = TRUE)
    }

    outside <- values < fitted_range[1] | values > fitted_range[2]
    if (!any(outside)) return(invisible(TRUE))

    msg <- paste0(
      "`", label, "` for variable '", variable,
      "' contains value(s) outside the fitted exposure range [",
      format(fitted_range[1]), ", ", format(fitted_range[2]), "]."
    )

    if (identical(extrapolation, "error")) {
      stop(msg, call. = FALSE)
    }
    if (identical(extrapolation, "warn")) {
      warning(msg, call. = FALSE)
    }

    invisible(TRUE)
  }

  # Reconstruct a crossbasis object only when the exact stored object is not
  # available. The fitted argvar/arglag are reused; they are never estimated
  # from `data`.
  get_effect_basis <- function(variable) {
    if (!is.null(metadata$basis_objects) &&
        inherits(metadata$basis_objects[[variable]], "crossbasis")) {
      return(metadata$basis_objects[[variable]])
    }

    if (!requireNamespace("dlnm", quietly = TRUE)) {
      stop(
        "Package 'dlnm' is required by `summarise_effects()`.",
        call. = FALSE
      )
    }

    definition <- .epix_basis_definition(metadata, variable)
    ids <- unique(as.character(data$epi_id))
    pieces <- vector("list", length(ids))

    for (i in seq_along(ids)) {
      idx <- which(as.character(data$epi_id) == ids[i])
      idx <- idx[order(data$time[idx])]
      values <- data[[variable]][idx]

      # NA separation prevents an exposure history from one epidemic from
      # contributing to the beginning of the next epidemic when a vector
      # crossbasis is reconstructed only to recover the fitted basis object.
      pieces[[i]] <- c(
        values,
        rep(NA_real_, definition$max_lag)
      )
    }

    pooled <- unlist(pieces, use.names = FALSE)

    dlnm::crossbasis(
      pooled,
      lag = definition$max_lag,
      argvar = definition$argvar,
      arglag = definition$arglag
    )
  }

  get_at_values <- function(variable) {
    custom <- NULL

    if (is.numeric(at)) {
      if (length(variables) != 1L) {
        stop(
          "A numeric `at` can be used only when exactly one exposure is ",
          "requested. For multiple variables, use a named list.",
          call. = FALSE
        )
      }
      custom <- at
    } else if (is.list(at) && variable %in% names(at)) {
      custom <- at[[variable]]
    }

    if (is.null(custom)) {
      values <- sort(unique(as.numeric(stats::quantile(
        data[[variable]],
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      ))))
    } else {
      if (!is.numeric(custom) || !length(custom) || anyNA(custom) ||
          any(!is.finite(custom))) {
        stop(
          "`at[['", variable,
          "']]` must contain finite numeric values.",
          call. = FALSE
        )
      }
      if (anyDuplicated(custom)) {
        stop(
          "`at[['", variable, "']]` must contain unique values.",
          call. = FALSE
        )
      }
      values <- sort(as.numeric(custom))
    }

    if (!length(values) || any(!is.finite(values))) {
      stop(
        "No valid exposure evaluation values are available for variable '",
        variable, "'.",
        call. = FALSE
      )
    }

    handle_extrapolation(values, variable, "at")
    values
  }

  # VALIDATE `at`

  if (is.list(at)) {
    if (!length(at) || is.null(names(at)) || anyNA(names(at)) ||
        any(!nzchar(names(at))) || anyDuplicated(names(at))) {
      stop(
        "When supplied as a list, `at` must be a non-empty named list with ",
        "unique, non-empty variable names.",
        call. = FALSE
      )
    }

    unknown_at <- setdiff(names(at), variables)
    if (length(unknown_at)) {
      stop(
        "`at` contains variable(s) not requested in `vars`: ",
        paste(unknown_at, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # JOINT REFERENCE PROFILE

  if (!is.list(ref) || !length(ref)) {
    stop("`ref` must be a non-empty list.", call. = FALSE)
  }

  ref_is_method <- "method" %in% names(ref)

  if (ref_is_method) {
    if (!is.character(ref$method) || length(ref$method) != 1L ||
        is.na(ref$method) || !nzchar(ref$method)) {
      stop(
        "Method-based `ref` must contain one valid character `method`.",
        call. = FALSE
      )
    }

    ref_method <- match.arg(
      ref$method,
      choices = c("median", "percentile", "fixed")
    )

    if (identical(ref_method, "percentile")) {
      if (!is.numeric(ref$value) || length(ref$value) != 1L ||
          is.na(ref$value) || !is.finite(ref$value) ||
          ref$value < 0 || ref$value > 1) {
        stop(
          "For `ref$method = 'percentile'`, `ref$value` must be one ",
          "probability between 0 and 1.",
          call. = FALSE
        )
      }
    }

    if (identical(ref_method, "fixed")) {
      if (!is.numeric(ref$value) || length(ref$value) != 1L ||
          is.na(ref$value) || !is.finite(ref$value)) {
        stop(
          "For `ref$method = 'fixed'`, `ref$value` must be one finite ",
          "numeric value.",
          call. = FALSE
        )
      }

      if (length(metadata$vars) > 1L) {
        warning(
          "Method-based `ref$method = 'fixed'` applies the same numeric ",
          "reference to every fitted exposure. For multivariable models with ",
          "different exposure units, a named exposure-specific `ref` list is ",
          "usually preferable.",
          call. = FALSE
        )
      }
    }

    reference_values <- vapply(
      metadata$vars,
      function(variable) {
        x <- data[[variable]]

        switch(
          ref_method,
          median = stats::median(x, na.rm = TRUE),
          percentile = as.numeric(stats::quantile(
            x,
            probs = ref$value,
            na.rm = TRUE,
            names = FALSE
          )),
          fixed = as.numeric(ref$value)
        )
      },
      numeric(1)
    )
  } else {
    if (is.null(names(ref)) || anyNA(names(ref)) ||
        any(!nzchar(names(ref))) || anyDuplicated(names(ref))) {
      stop(
        "Exposure-specific `ref` must be a named list with unique, ",
        "non-empty fitted exposure names.",
        call. = FALSE
      )
    }

    unknown_ref <- setdiff(names(ref), metadata$vars)
    missing_ref <- setdiff(metadata$vars, names(ref))

    if (length(unknown_ref)) {
      stop(
        "`ref` contains variable(s) not found in the fitted model: ",
        paste(unknown_ref, collapse = ", "), ".",
        call. = FALSE
      )
    }

    if (length(missing_ref)) {
      stop(
        "Exposure-specific `ref` must contain one value for every fitted ",
        "exposure. Missing: ", paste(missing_ref, collapse = ", "), ".",
        call. = FALSE
      )
    }

    ref <- ref[metadata$vars]

    valid <- vapply(
      ref,
      function(x) {
        is.numeric(x) && length(x) == 1L &&
          !is.na(x) && is.finite(x)
      },
      logical(1)
    )

    if (any(!valid)) {
      stop(
        "Each exposure-specific reference value must be one finite numeric ",
        "value.",
        call. = FALSE
      )
    }

    reference_values <- vapply(ref, as.numeric, numeric(1))
  }

  names(reference_values) <- metadata$vars

  for (variable in metadata$vars) {
    handle_extrapolation(
      reference_values[[variable]],
      variable,
      "reference"
    )
  }

  reference_profiles <- lapply(
    metadata$vars,
    function(variable) {
      rep(
        reference_values[[variable]],
        metadata$history_length
      )
    }
  )
  names(reference_profiles) <- metadata$vars

  baseline_design <- .build_newdata_basis(
    fit = fit,
    profiles = reference_profiles,
    extrapolation = extrapolation
  )

  link_object <- .epix_link_object(metadata$link)

  central_parameters <- .extract_central_parameters(fit)

  baseline_eta_point <- as.numeric(
    .predict_point_population(
      fit = fit,
      newdata = baseline_design,
      type = "link",
      warn_fallback = FALSE
    )
  )[1L]

  baseline_response_point <- as.numeric(
    link_object$linkinv(baseline_eta_point)
  )

  if (!is.finite(baseline_response_point)) {
    stop(
      "The joint reference profile produced a non-finite baseline expected ",
      "response.",
      call. = FALSE
    )
  }

  # One joint draw matrix is generated ONCE and reused for every requested
  # variable. Thus sample s has the same meaning throughout the entire result.
  parameter_draws <- NULL
  baseline_eta_draws <- NULL
  baseline_response_draws <- NULL
  baseline_summary <- NULL

  if (uncertainty) {
    parameter_draws <- .extract_parameter_draws(
      fit = fit,
      n_samples = n_samples
    )

    baseline_eta_draws <- as.numeric(
      .predict_draws_population(
        fit = fit,
        newdata = baseline_design,
        n_samples = n_samples,
        type = "link",
        parameter_draws = parameter_draws
      )[, 1L]
    )

    baseline_response_draws <- as.numeric(
      link_object$linkinv(baseline_eta_draws)
    )

    if (any(!is.finite(baseline_response_draws))) {
      stop(
        "Draw-by-draw baseline prediction produced non-finite values.",
        call. = FALSE
      )
    }

    baseline_summary <- .summarise_prediction_draws(
      draws = matrix(baseline_response_draws, ncol = 1L),
      probs = interval_probs,
      value_name = "baseline"
    )
  }

  # CROSSPRED-BASED EFFECT DESIGN

  crosspred_eta_matrix <- function(
    basis,
    coefficients,
    at_values,
    center_value
  ) {
    if (!requireNamespace("dlnm", quietly = TRUE)) {
      stop(
        "Package 'dlnm' is required by `summarise_effects()`.",
        call. = FALSE
      )
    }

    coefficients <- as.numeric(coefficients)
    p <- ncol(basis)

    if (length(coefficients) != p) {
      stop(
        "Cross-basis coefficient dimension mismatch: basis has ", p,
        " column(s), but ", length(coefficients),
        " coefficient(s) were supplied.",
        call. = FALSE
      )
    }

    V0 <- matrix(0, nrow = p, ncol = p)

    prediction <- dlnm::crosspred(
      basis = basis,
      coef = coefficients,
      vcov = V0,
      model.link = metadata$link,
      at = at_values,
      cen = center_value,
      bylag = 1,
      cumul = FALSE
    )

    eta <- as.matrix(prediction$matfit)

    if (!nrow(eta) || !ncol(eta) || any(!is.finite(eta))) {
      stop(
        "`dlnm::crosspred()` did not return a finite lag-specific association ",
        "matrix.",
        call. = FALSE
      )
    }

    lag_index <- suppressWarnings(
      as.integer(gsub("^lag", "", colnames(eta)))
    )

    if (is.null(colnames(eta)) || anyNA(lag_index)) {
      lag_index <- seq.int(0L, ncol(eta) - 1L)
    }

    lag_order <- order(lag_index)
    list(
      eta = eta[, lag_order, drop = FALSE],
      lag = lag_index[lag_order]
    )
  }

  transform_lag_matrix <- function(
    eta_matrix,
    lag_index,
    at_values,
    variable
  ) {
    if (identical(scale, "lag")) {
      grid <- expand.grid(
        value = at_values,
        lag = lag_index,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      grid$scale <- "lag"
      grid <- grid[, c("lag", "scale", "value"), drop = FALSE]

      return(list(
        matrix = eta_matrix,
        grid = grid
      ))
    }

    if (incremental) {
      cumulative <- t(apply(eta_matrix, 1L, cumsum))

      grid <- expand.grid(
        value = at_values,
        lag = lag_index,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      grid$period <- paste0("0-", grid$lag)
      grid$scale <- "period"
      grid <- grid[, c("lag", "period", "scale", "value"), drop = FALSE]

      return(list(
        matrix = cumulative,
        grid = grid
      ))
    }

    if (any(lag_periods$lag_end > max(lag_index))) {
      stop(
        "`lag_periods` exceeds the fitted maximum lag for variable '",
        variable, "'. Maximum available lag: ", max(lag_index), ".",
        call. = FALSE
      )
    }

    period_matrix <- vapply(
      seq_len(nrow(lag_periods)),
      function(i) {
        selected <- which(
          lag_index >= lag_periods$lag_start[i] &
            lag_index <= lag_periods$lag_end[i]
        )

        if (!length(selected)) {
          stop(
            "No fitted lags fall inside period '",
            lag_periods$period[i], "'.",
            call. = FALSE
          )
        }

        rowSums(eta_matrix[, selected, drop = FALSE])
      },
      numeric(nrow(eta_matrix))
    )

    if (is.null(dim(period_matrix))) {
      period_matrix <- matrix(
        period_matrix,
        nrow = nrow(eta_matrix),
        ncol = nrow(lag_periods)
      )
    }

    grid <- expand.grid(
      value = at_values,
      period = lag_periods$period,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    grid$scale <- "period"
    grid <- grid[, c("period", "scale", "value"), drop = FALSE]

    list(
      matrix = period_matrix,
      grid = grid
    )
  }

  # Builds the exact centered prediction design used by dlnm::crosspred().
  # Because matfit is linear in the supplied cross-basis coefficients, running
  # crosspred with unit coefficient vectors recovers the columns of that design.
  # Full parameter draws can then be propagated with one matrix multiplication.
  build_effect_design <- function(
    basis,
    at_values,
    center_value,
    lag_index_reference,
    transformed_grid,
    variable
  ) {
    p <- ncol(basis)
    design_columns <- vector("list", p)

    for (j in seq_len(p)) {
      unit <- numeric(p)
      unit[j] <- 1

      unit_result <- crosspred_eta_matrix(
        basis = basis,
        coefficients = unit,
        at_values = at_values,
        center_value = center_value
      )

      if (!identical(unit_result$lag, lag_index_reference)) {
        stop(
          "Internal lag-order mismatch while constructing the DLNM effect ",
          "design for variable '", variable, "'.",
          call. = FALSE
        )
      }

      transformed <- transform_lag_matrix(
        eta_matrix = unit_result$eta,
        lag_index = unit_result$lag,
        at_values = at_values,
        variable = variable
      )

      if (nrow(transformed$grid) != nrow(transformed_grid) ||
          !identical(names(transformed$grid), names(transformed_grid))) {
        stop(
          "Internal grid mismatch while constructing the DLNM effect design.",
          call. = FALSE
        )
      }

      design_columns[[j]] <- as.vector(transformed$matrix)
    }

    A <- do.call(cbind, design_columns)
    storage.mode(A) <- "double"

    if (any(!is.finite(A))) {
      stop(
        "The centered DLNM effect design contains non-finite values.",
        call. = FALSE
      )
    }

    A
  }

  # DIAGNOSTICS

  compute_one_diagnostic <- function(df) {
    lag_values <- sort(unique(df$lag))
    rows <- vector("list", length(lag_values))

    for (i in seq_along(lag_values)) {
      current_lag <- lag_values[i]
      values <- df$eta[df$lag == current_lag]
      values <- values[is.finite(values)]

      if (!length(values)) {
        rows[[i]] <- data.frame(
          lag = current_lag,
          mean_effect = NA_real_,
          max_effect = NA_real_,
          max_abs_effect = NA_real_,
          absmean_effect = NA_real_,
          eta_lag = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      max_pos <- which.max(abs(values))

      rows[[i]] <- data.frame(
        lag = current_lag,
        mean_effect = mean(values),
        max_effect = values[max_pos],
        max_abs_effect = max(abs(values)),
        absmean_effect = mean(abs(values)),
        eta_lag = sum(values),
        stringsAsFactors = FALSE
      )
    }

    out <- do.call(rbind, rows)

    out$score_max <- out$max_abs_effect
    out$score_mean <- abs(out$mean_effect)
    out$score_absmean <- out$absmean_effect

    out$rank_max <- rank(-out$score_max, ties.method = "min", na.last = "keep")
    out$rank_mean <- rank(-out$score_mean, ties.method = "min", na.last = "keep")
    out$rank_absmean <- rank(
      -out$score_absmean,
      ties.method = "min",
      na.last = "keep"
    )

    absolute_denominator <- sum(abs(out$eta_lag), na.rm = TRUE)
    signed_denominator <- sum(out$eta_lag, na.rm = TRUE)

    tolerance <- sqrt(.Machine$double.eps) *
      max(1, absolute_denominator)

    signed_stable <- is.finite(signed_denominator) &&
      abs(signed_denominator) > tolerance

    out$contribution_absolute <- if (
      is.finite(absolute_denominator) &&
      absolute_denominator > 0
    ) {
      abs(out$eta_lag) / absolute_denominator
    } else {
      NA_real_
    }

    out$contribution_absolute_percent <- 100 * out$contribution_absolute

    out$contribution_signed <- if (signed_stable) {
      out$eta_lag / signed_denominator
    } else {
      NA_real_
    }

    out$contribution_signed_percent <- 100 * out$contribution_signed
    out$signed_contribution_available <- signed_stable

    out
  }

  compute_lag_diagnostics <- function(diag_data) {
    if (!all(c("var", "lag", "eta") %in% names(diag_data))) {
      stop(
        "Internal diagnostics require `var`, `lag`, and `eta`.",
        call. = FALSE
      )
    }

    has_samples <- "sample" %in% names(diag_data)

    if (!has_samples) {
      split_var <- split(diag_data, diag_data$var, drop = TRUE)
      result <- lapply(names(split_var), function(variable) {
        current <- compute_one_diagnostic(split_var[[variable]])
        current$var <- variable
        current[, c("var", setdiff(names(current), "var")), drop = FALSE]
      })

      out <- do.call(rbind, result)
      rownames(out) <- NULL
      out <- out[order(out$var, out$rank_absmean, out$lag), , drop = FALSE]
      return(out)
    }

    keys <- unique(diag_data[, c("sample", "var"), drop = FALSE])
    per_sample <- vector("list", nrow(keys))

    for (i in seq_len(nrow(keys))) {
      current_data <- diag_data[
        diag_data$sample == keys$sample[i] &
          diag_data$var == keys$var[i],
        ,
        drop = FALSE
      ]

      current <- compute_one_diagnostic(current_data)
      current$sample <- keys$sample[i]
      current$var <- keys$var[i]
      per_sample[[i]] <- current[, c(
        "sample", "var",
        setdiff(names(current), c("sample", "var"))
      ), drop = FALSE]
    }

    samples <- do.call(rbind, per_sample)
    rownames(samples) <- NULL

    numeric_metrics <- c(
      "mean_effect", "max_effect", "max_abs_effect", "absmean_effect",
      "eta_lag", "score_max", "score_mean", "score_absmean",
      "rank_max", "rank_mean", "rank_absmean",
      "contribution_absolute", "contribution_absolute_percent",
      "contribution_signed", "contribution_signed_percent"
    )

    groups <- unique(samples[, c("var", "lag"), drop = FALSE])
    summary_rows <- vector("list", nrow(groups))

    for (i in seq_len(nrow(groups))) {
      current <- samples[
        samples$var == groups$var[i] &
          samples$lag == groups$lag[i],
        ,
        drop = FALSE
      ]

      row <- data.frame(
        var = groups$var[i],
        lag = groups$lag[i],
        stringsAsFactors = FALSE
      )

      for (metric in numeric_metrics) {
        x <- current[[metric]]
        q <- safe_quantile(x)

        row[[metric]] <- stats::median(x, na.rm = TRUE)
        row[[paste0(metric, "_sd")]] <- safe_sd(x)
        row[[paste0(metric, "_lower")]] <- q[1]
        row[[paste0(metric, "_upper")]] <- q[2]
      }

      row$signed_contribution_available <- all(
        current$signed_contribution_available,
        na.rm = TRUE
      )

      summary_rows[[i]] <- row
    }

    out <- do.call(rbind, summary_rows)
    rownames(out) <- NULL
    out <- out[order(out$var, out$rank_absmean, out$lag), , drop = FALSE]
    out
  }

  # ONE VARIABLE

  summarise_one_variable <- function(variable) {
    at_values <- get_at_values(variable)
    center_value <- as.numeric(reference_values[[variable]])
    basis <- get_effect_basis(variable)

    cb_names <- .epix_cb_cols_for_var(metadata$cb_cols, variable)

    if (length(cb_names) != ncol(basis)) {
      stop(
        "Cross-basis dimension mismatch for variable '", variable,
        "': fitted metadata contain ", length(cb_names),
        " coefficient(s), but the effect basis contains ",
        ncol(basis), " column(s).",
        call. = FALSE
      )
    }

    missing_central <- setdiff(cb_names, names(central_parameters))
    if (length(missing_central)) {
      stop(
        "Central parameter vector is missing cross-basis coefficient(s) for '",
        variable, "': ", paste(missing_central, collapse = ", "), ".",
        call. = FALSE
      )
    }

    point_lag <- crosspred_eta_matrix(
      basis = basis,
      coefficients = central_parameters[cb_names],
      at_values = at_values,
      center_value = center_value
    )

    max_lag_expected <- metadata$max_lag
    if (max(point_lag$lag) != max_lag_expected) {
      stop(
        "DLNM lag mismatch for variable '", variable,
        "': crosspred returned maximum lag ", max(point_lag$lag),
        " but fitted metadata specify ", max_lag_expected, ".",
        call. = FALSE
      )
    }

    point_transformed <- transform_lag_matrix(
      eta_matrix = point_lag$eta,
      lag_index = point_lag$lag,
      at_values = at_values,
      variable = variable
    )

    grid <- point_transformed$grid
    grid$var <- variable
    grid <- grid[, c(
      "var",
      setdiff(names(grid), "var")
    ), drop = FALSE]

    eta_point <- as.vector(point_transformed$matrix)
    effect_point <- transform_effect(eta_point)

    predicted_point <- as.numeric(
      link_object$linkinv(baseline_eta_point + eta_point)
    )
    delta_point <- predicted_point - baseline_response_point

    if (any(!is.finite(predicted_point)) || any(!is.finite(delta_point))) {
      stop(
        "Response-scale point effects produced non-finite values for variable '",
        variable, "'.",
        call. = FALSE
      )
    }

    point_output <- grid
    point_output$eta <- eta_point
    point_output$effect <- as.numeric(effect_point)
    point_output$baseline <- baseline_response_point
    point_output$predicted <- predicted_point
    point_output$delta <- delta_point

    if (!uncertainty) {
      diagnostic_source <- if (diagnostics) {
        point_output[, intersect(
          c("var", "lag", "value", "eta"),
          names(point_output)
        ), drop = FALSE]
      } else {
        NULL
      }

      return(list(
        effects = point_output,
        diagnostics_source = diagnostic_source
      ))
    }

    # Full joint parameter uncertainty

    missing_draws <- setdiff(cb_names, colnames(parameter_draws))
    if (length(missing_draws)) {
      stop(
        "Parameter-draw matrix is missing cross-basis coefficient(s) for '",
        variable, "': ", paste(missing_draws, collapse = ", "), ".",
        call. = FALSE
      )
    }

    A <- build_effect_design(
      basis = basis,
      at_values = at_values,
      center_value = center_value,
      lag_index_reference = point_lag$lag,
      transformed_grid = point_transformed$grid,
      variable = variable
    )

    colnames(A) <- cb_names

    eta_draws <- parameter_draws[, cb_names, drop = FALSE] %*% t(A)
    storage.mode(eta_draws) <- "double"

    if (any(!is.finite(eta_draws))) {
      stop(
        "DLNM effect draws contain non-finite values for variable '",
        variable, "'.",
        call. = FALSE
      )
    }

    effect_draws <- transform_effect(eta_draws)

    target_eta_draws <- sweep(
      eta_draws,
      MARGIN = 1L,
      STATS = baseline_eta_draws,
      FUN = "+"
    )

    predicted_draws <- matrix(
      as.numeric(link_object$linkinv(as.vector(target_eta_draws))),
      nrow = nrow(target_eta_draws),
      ncol = ncol(target_eta_draws),
      byrow = FALSE
    )

    delta_draws <- sweep(
      predicted_draws,
      MARGIN = 1L,
      STATS = baseline_response_draws,
      FUN = "-"
    )

    if (any(!is.finite(effect_draws)) ||
        any(!is.finite(predicted_draws)) ||
        any(!is.finite(delta_draws))) {
      stop(
        "Draw-by-draw effect transformation produced non-finite values for ",
        "variable '", variable, "'.",
        call. = FALSE
      )
    }

    diagnostic_source <- NULL

    if (identical(output, "samples")) {
      n_effect_rows <- nrow(grid)
      n_draws <- nrow(eta_draws)

      long_grid <- grid[
        rep(seq_len(n_effect_rows), times = n_draws),
        ,
        drop = FALSE
      ]
      rownames(long_grid) <- NULL

      sample_output <- long_grid
      sample_output$sample <- rep(seq_len(n_draws), each = n_effect_rows)
      sample_output <- sample_output[, c(
        "sample",
        setdiff(names(sample_output), "sample")
      ), drop = FALSE]

      sample_output$eta <- as.vector(t(eta_draws))
      sample_output$effect <- as.vector(t(effect_draws))
      sample_output$baseline <- rep(
        baseline_response_draws,
        each = n_effect_rows
      )
      sample_output$predicted <- as.vector(t(predicted_draws))
      sample_output$delta <- as.vector(t(delta_draws))

      if (diagnostics) {
        diagnostic_source <- sample_output[, intersect(
          c("sample", "var", "lag", "value", "eta"),
          names(sample_output)
        ), drop = FALSE]
      }

      return(list(
        effects = sample_output,
        diagnostics_source = diagnostic_source
      ))
    }

    # Summary output is calculated directly from the draw matrices after all
    # transformations. No fixed baseline is inserted.
    eta_summary <- .summarise_prediction_draws(
      draws = eta_draws,
      probs = interval_probs,
      value_name = "eta"
    )
    effect_summary <- .summarise_prediction_draws(
      draws = effect_draws,
      probs = interval_probs,
      value_name = "effect"
    )
    predicted_summary <- .summarise_prediction_draws(
      draws = predicted_draws,
      probs = interval_probs,
      value_name = "predicted"
    )
    delta_summary <- .summarise_prediction_draws(
      draws = delta_draws,
      probs = interval_probs,
      value_name = "delta"
    )

    summary_output <- grid

    summary_output <- cbind(
      summary_output,
      eta_summary[, setdiff(names(eta_summary), "row_id"), drop = FALSE],
      effect_summary[, setdiff(names(effect_summary), "row_id"), drop = FALSE]
    )

    baseline_values <- baseline_summary[
      rep(1L, nrow(summary_output)),
      setdiff(names(baseline_summary), "row_id"),
      drop = FALSE
    ]
    rownames(baseline_values) <- NULL

    summary_output <- cbind(
      summary_output,
      baseline_values,
      predicted_summary[
        ,
        setdiff(names(predicted_summary), "row_id"),
        drop = FALSE
      ],
      delta_summary[
        ,
        setdiff(names(delta_summary), "row_id"),
        drop = FALSE
      ]
    )

    if (diagnostics) {
      n_effect_rows <- nrow(grid)
      n_draws <- nrow(eta_draws)

      diagnostic_source <- grid[
        rep(seq_len(n_effect_rows), times = n_draws),
        ,
        drop = FALSE
      ]
      rownames(diagnostic_source) <- NULL

      diagnostic_source$sample <- rep(
        seq_len(n_draws),
        each = n_effect_rows
      )
      diagnostic_source$eta <- as.vector(t(eta_draws))

      diagnostic_source <- diagnostic_source[, intersect(
        c("sample", "var", "lag", "value", "eta"),
        names(diagnostic_source)
      ), drop = FALSE]
    }

    list(
      effects = as.data.frame(summary_output),
      diagnostics_source = diagnostic_source
    )
  }

  # RUN REQUESTED VARIABLES

  variable_results <- lapply(variables, summarise_one_variable)

  output_data <- do.call(
    rbind,
    lapply(variable_results, `[[`, "effects")
  )
  rownames(output_data) <- NULL

  # Metadata on the output make the interpretation recoverable downstream.
  attr(output_data, "epiexposure_reference") <- reference_values
  attr(output_data, "epiexposure_reference_type") <- if (ref_is_method) {
    paste0("method:", ref_method)
  } else {
    "exposure_specific"
  }
  attr(output_data, "epiexposure_effect_measure") <- effect_measure
  attr(output_data, "epiexposure_link") <- metadata$link
  attr(output_data, "epiexposure_prediction_level") <- "population"
  attr(output_data, "epiexposure_prediction_estimand") <- "expected_response"
  attr(output_data, "epiexposure_uncertainty") <- uncertainty
  attr(output_data, "epiexposure_scale") <- scale
  attr(output_data, "epiexposure_max_lag") <- metadata$max_lag
  attr(output_data, "epiexposure_history_length") <- metadata$history_length
  attr(output_data, "epiexposure_history_contract") <- metadata$history_contract

  if (uncertainty) {
    attr(output_data, "epiexposure_n_samples") <- n_samples
    attr(output_data, "epiexposure_interval_probs") <- interval_probs
    attr(output_data, "epiexposure_summary_center") <- if (
      identical(output, "summary")
    ) {
      "median"
    } else {
      "samples"
    }
  }

  if (!diagnostics) {
    return(output_data)
  }

  diagnostics_source <- do.call(
    rbind,
    lapply(variable_results, `[[`, "diagnostics_source")
  )
  rownames(diagnostics_source) <- NULL

  diagnostics_data <- compute_lag_diagnostics(diagnostics_source)

  result <- structure(
    list(
      effects = output_data,
      diagnostics = diagnostics_data
    ),
    class = c("epiexposure_effects", "list")
  )

  attr(result, "epiexposure_reference") <- reference_values
  attr(result, "epiexposure_effect_measure") <- effect_measure
  attr(result, "epiexposure_link") <- metadata$link
  attr(result, "epiexposure_prediction_level") <- "population"
  attr(result, "epiexposure_prediction_estimand") <- "expected_response"
  attr(result, "epiexposure_uncertainty") <- uncertainty
  attr(result, "epiexposure_max_lag") <- metadata$max_lag
  attr(result, "epiexposure_history_length") <- metadata$history_length
  attr(result, "epiexposure_history_contract") <- metadata$history_contract

  result
}
