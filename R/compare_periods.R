#' Compare period-specific DLNM effects
#'
#' Compares period-specific DLNM associations returned by
#' `summarise_effects(scale = "period")` against one selected reference period.
#'
#' The comparison preserves the EpiExposure effect contract:
#'
#' - `eta` is the fundamental additive DLNM contrast on the fitted link scale;
#' - `effect` is the user-requested representation of `eta`
#'   (`"linear"`, `"exponentiated"`, or `"percent"`);
#' - `predicted` is the population-level expected outcome for that period;
#' - `delta = predicted - baseline` is the absolute response-scale change from
#'   the joint exposure reference profile.
#'
#' Period comparisons are always made at the same exposure variable and the
#' same exposure evaluation value.
#'
#' @param period_df Output from the current
#'   `summarise_effects(scale = "period")`. A data frame is expected. If an
#'   `"epiexposure_effects"` object produced with diagnostics is supplied, its
#'   `$effects` component is used. The object must carry the current
#'   EpiExposure temporal metadata: one common fitted `max_lag`, exact
#'   `history_length = max_lag + 1`, and the exact-history contract.
#'
#'   Current EpiExposure period output must contain `var`, `value`, `period`,
#'   `scale`, `eta`, `effect`, `baseline`, `predicted`, and `delta`.
#'
#'   If `summarise_effects(..., uncertainty = TRUE)` was used,
#'   `period_df` must be the **sample-level** output produced with
#'   `output = "samples"`. Comparing already-summarized medians and interval
#'   bounds would not preserve covariance among period effects.
#' @param period_ref Reference period. Either:
#'
#'   - one positive integer selecting the period by its defined order; or
#'   - one exact non-empty character period label, for example `"W1"`.
#'
#'   Numeric references do not require period labels to end in numbers. If
#'   `periods` is supplied, its row order defines the numeric period order.
#'   Otherwise the first-occurrence order in `period_df$period` is used.
#' @param periods Optional period-definition data frame, typically returned by
#'   `define_periods()`, containing `period`, `lag_start`, and `lag_end`.
#'
#'   When supplied, period labels are validated against `period_df`, every
#'   period boundary must lie within the inherited common fitted lag range
#'   `0:max_lag`, and `lag_start`, `lag_end`, and `n_lags` are appended to the
#'   comparison output. The number of lags is
#'
#'   `n_lags = lag_end - lag_start + 1`.
#'
#'   Supplying `periods` is recommended when the compared periods have unequal
#'   widths because cumulative DLNM effects are sums across included lags.
#' @param output Character. `"summary"` or `"samples"`.
#'
#'   For deterministic input (`uncertainty = FALSE`), only `"summary"` is
#'   available and contains the direct period comparisons.
#'
#'   For sample-level uncertain input, `"samples"` returns draw-by-draw
#'   comparisons and `"summary"` returns medians, SDs, and empirical intervals
#'   calculated **after** the draw-by-draw comparison.
#' @param interval_probs Optional numeric vector of length two defining the
#'   empirical interval for uncertain comparisons. If `NULL`, the interval
#'   stored by `summarise_effects()` is reused; otherwise the default
#'   `c(0.025, 0.975)` is used when no stored interval is available.
#' @param eps Positive finite numerical tolerance. It is used only for
#'   consistency checks of values that should be equal under the EpiExposure
#'   contract; it is not used to alter effects or to force a denominator away
#'   from zero.
#'
#' @return A data frame containing the original period-specific quantities and
#'   the corresponding reference-period quantities. Principal comparison
#'   columns are:
#'
#'   \describe{
#'     \item{`eta_diff`}{Difference in additive DLNM contrasts:
#'
#'       `eta - ref_eta`.
#'
#'       This is the canonical comparison on the fitted link scale.}
#'     \item{`diff`}{Arithmetic difference in the displayed `effect`:
#'
#'       `effect - ref_effect`.
#'
#'       Its units depend on the `effect_measure` used by
#'       `summarise_effects()`.}
#'     \item{`ratio`}{For `log` links, the ratio of the two period-specific
#'       response ratios. For `logit` links, the ratio of the two
#'       period-specific odds ratios:
#'
#'       `exp(eta - ref_eta)`.
#'
#'       For other links, no general multiplicative interpretation is defined,
#'       so `ratio` and `ratio_percent` are `NA`.}
#'     \item{`ratio_percent`}{`100 * (ratio - 1)` when `ratio` is defined.}
#'     \item{`predicted_diff`}{Difference between the period-specific
#'       population expected responses:
#'
#'       `predicted - ref_predicted`.}
#'     \item{`delta_diff`}{Difference between period-specific absolute
#'       response-scale changes:
#'
#'       `delta - ref_delta`.
#'
#'       Because the joint `baseline` is identical across compared periods,
#'       `delta_diff` and `predicted_diff` should be numerically equal.}
#'   }
#'
#'   With uncertain sample-level input and `output = "summary"`, each numeric
#'   quantity is summarized by its median plus `_sd`, `_lower`, and `_upper`
#'   columns.
#'
#'   Output attributes record the reference period, fitted link,
#'   effect measure, prediction contract, uncertainty contract, ratio
#'   interpretation, common fitted `max_lag`, exact history length, and the
#'   inherited EpiExposure exact-history contract.
#'
#' @details
#' ## Exact common lag/history contract
#'
#' `compare_periods()` does not receive the original long-format exposure
#' histories, so it cannot recount time rows itself. Instead, it requires the
#' temporal metadata propagated by the current `summarise_effects()` output:
#'
#' \deqn{
#'   history\_length = max\_lag + 1.
#' }
#'
#' The stored contract must state that all fitted exposure variables used the
#' same exact lag window. This prevents period comparisons from silently mixing
#' effect objects created under incompatible temporal histories.
#'
#' When `periods` is supplied, its lag boundaries are additionally checked
#' against the inherited common fitted `max_lag`. `compare_periods()` does not
#' require a manually supplied `periods` object to cover every fitted lag,
#' because analysts may intentionally compare a subset of defined periods; it
#' only requires every supplied period to remain inside the fitted lag window.
#'
#' ## What is being compared?
#'
#' `summarise_effects(scale = "period")` first sums lag-specific DLNM
#' contributions on the additive linear-predictor scale within each requested
#' period. `compare_periods()` then compares those period-specific cumulative
#' contrasts at the **same exposure value**.
#'
#' If period \eqn{j} has contrast \eqn{\eta_j(x)} and the selected reference
#' period has contrast \eqn{\eta_r(x)}, the fundamental comparison is
#'
#' \deqn{\Delta\eta_{j:r}(x) = \eta_j(x) - \eta_r(x).}
#'
#' Under a log link,
#'
#' \deqn{\exp\{\Delta\eta_{j:r}(x)\}
#'       = RR_j(x) / RR_r(x),}
#'
#' while under a logit link,
#'
#' \deqn{\exp\{\Delta\eta_{j:r}(x)\}
#'       = OR_j(x) / OR_r(x).}
#'
#' Therefore the returned multiplicative `ratio` is derived from `eta_diff`,
#' not from blindly dividing whichever transformed `effect` happens to be
#' displayed.
#'
#' This distinction matters for `effect_measure = "percent"`. For example,
#' a +20 percent relative effect and a +10 percent relative effect do **not**
#' imply a period ratio of `20 / 10 = 2`. Their underlying multiplicative
#' effects are 1.20 and 1.10, and the appropriate relative ratio is
#' `1.20 / 1.10`.
#'
#' ## Response-scale comparison
#'
#' `predicted_diff` compares expected outcomes under the same joint reference
#' profile while changing the focal exposure over different lag periods. For
#' Beta and Binomial outcomes, multiplying `predicted_diff` by 100 gives the
#' difference in percentage points between the two period-specific expected
#' proportions/probabilities.
#'
#' Since `delta = predicted - baseline` and the same baseline is used for every
#' period comparison at a given parameter draw, `delta_diff` must equal
#' `predicted_diff` up to numerical tolerance. The function validates this
#' identity.
#'
#' ## Unequal period widths
#'
#' Period effects are cumulative sums over lags. Consequently, comparing a
#' 21-lag period with a 10-lag period compares the cumulative associations as
#' defined; the contrast can reflect both the lag-response pattern and the
#' different number of included lags. `compare_periods()` does not divide
#' cumulative effects by period width. Use period definitions with the
#' scientific interpretation intended by the analysis.
#'
#' ## Uncertainty
#'
#' When `summarise_effects()` uses uncertainty, all requested effects are based
#' on the same joint parameter draw. `compare_periods()` preserves this
#' covariance by matching the target and reference periods **within the same
#' draw**, calculating `eta_diff`, `diff`, `ratio`, `predicted_diff`, and
#' `delta_diff` draw by draw, and only then summarizing those comparisons.
#'
#' Already summarized uncertain output is therefore rejected explicitly.
#' Subtracting medians or combining marginal interval endpoints would not
#' recover the distribution of a difference or ratio.
#'
#' @export
compare_periods <- function(
    period_df,
    period_ref = 1,
    periods = NULL,
    output = c("summary", "samples"),
    interval_probs = NULL,
    eps = 1e-10
) {

  # ==========================================================================
  # SMALL VALIDATORS
  # ==========================================================================

  output <- match.arg(output)

  valid_scalar_string <- function(x) {
    is.character(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      nzchar(x)
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

  nearly_equal <- function(a, b, tolerance) {
    if (length(a) != length(b)) {
      return(rep(FALSE, max(length(a), length(b))))
    }

    abs(a - b) <= tolerance *
      pmax(1, abs(a), abs(b))
  }

  summarize_vector <- function(x, probs) {
    x <- as.numeric(x)

    if (!length(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      return(c(
        value = NA_real_,
        sd = NA_real_,
        lower = NA_real_,
        upper = NA_real_
      ))
    }

    c(
      value = stats::median(x),
      sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
      lower = stats::quantile(
        x,
        probs = probs[1L],
        names = FALSE,
        type = 7
      ),
      upper = stats::quantile(
        x,
        probs = probs[2L],
        names = FALSE,
        type = 7
      )
    )
  }

  # ==========================================================================
  # ACCEPT THE DIAGNOSTIC WRAPPER EXPLICITLY
  # ==========================================================================

  if (!is.data.frame(period_df) &&
      is.list(period_df) &&
      is.data.frame(period_df$effects)) {
    period_df <- period_df$effects
  }

  if (!is.data.frame(period_df) ||
      !nrow(period_df)) {
    stop(
      "`period_df` must be a non-empty data.frame returned by ",
      "`summarise_effects(scale = 'period')`.",
      call. = FALSE
    )
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

  tolerance <- max(
    as.numeric(eps),
    sqrt(.Machine$double.eps)
  )

  # ==========================================================================
  # STRICT SUMMARISE_EFFECTS METADATA CONTRACT
  # ==========================================================================

  scale_attr <- attr(
    period_df,
    "epiexposure_scale",
    exact = TRUE
  )
  effect_measure <- attr(
    period_df,
    "epiexposure_effect_measure",
    exact = TRUE
  )
  link_name <- attr(
    period_df,
    "epiexposure_link",
    exact = TRUE
  )
  prediction_level <- attr(
    period_df,
    "epiexposure_prediction_level",
    exact = TRUE
  )
  prediction_estimand <- attr(
    period_df,
    "epiexposure_prediction_estimand",
    exact = TRUE
  )
  uncertainty <- attr(
    period_df,
    "epiexposure_uncertainty",
    exact = TRUE
  )
  max_lag_metadata <- attr(
    period_df,
    "epiexposure_max_lag",
    exact = TRUE
  )
  history_length_metadata <- attr(
    period_df,
    "epiexposure_history_length",
    exact = TRUE
  )
  history_contract_metadata <- attr(
    period_df,
    "epiexposure_history_contract",
    exact = TRUE
  )

  if (is.null(scale_attr) ||
      !identical(scale_attr, "period")) {
    stop(
      "`period_df` must come from `summarise_effects(scale = 'period')` under ",
      "the current EpiExposure contract.",
      call. = FALSE
    )
  }

  if (!valid_scalar_string(effect_measure) ||
      !effect_measure %in%
      c("linear", "exponentiated", "percent")) {
    stop(
      "`period_df` is missing valid `epiexposure_effect_measure` metadata.",
      call. = FALSE
    )
  }

  if (!valid_scalar_string(link_name)) {
    stop(
      "`period_df` is missing valid fitted-link metadata.",
      call. = FALSE
    )
  }
  link_name <- tolower(link_name)

  if (!identical(prediction_level, "population")) {
    stop(
      "`compare_periods()` requires population/fixed-component effects.",
      call. = FALSE
    )
  }

  if (!identical(
    prediction_estimand,
    "expected_response"
  )) {
    stop(
      "`compare_periods()` requires the expected-response estimand.",
      call. = FALSE
    )
  }

  if (!is.logical(uncertainty) ||
      length(uncertainty) != 1L ||
      is.na(uncertainty)) {
    stop(
      "`period_df` is missing valid EpiExposure uncertainty metadata.",
      call. = FALSE
    )
  }

  if (!valid_integer_scalar(max_lag_metadata) ||
      max_lag_metadata < 0L) {
    stop(
      "`period_df` is missing valid common fitted `epiexposure_max_lag` ",
      "metadata from the current `summarise_effects()`.",
      call. = FALSE
    )
  }
  max_lag_metadata <- as.integer(max_lag_metadata)

  if (!valid_integer_scalar(history_length_metadata) ||
      history_length_metadata < 1L) {
    stop(
      "`period_df` is missing valid `epiexposure_history_length` metadata.",
      call. = FALSE
    )
  }
  history_length_metadata <- as.integer(history_length_metadata)

  if (!identical(
    history_length_metadata,
    max_lag_metadata + 1L
  )) {
    stop(
      "`period_df` temporal metadata are inconsistent: history length must ",
      "equal `max_lag + 1`. Expected ",
      max_lag_metadata + 1L,
      " but found ",
      history_length_metadata,
      ".",
      call. = FALSE
    )
  }

  if (!identical(
    history_contract_metadata,
    "all_fitted_exposures_same_exact_max_lag_plus_one"
  )) {
    stop(
      "`compare_periods()` requires period effects produced under the current ",
      "EpiExposure exact-history contract.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # REQUIRED COLUMNS
  # ==========================================================================

  required_columns <- c(
    "var",
    "value",
    "period",
    "scale",
    "eta",
    "effect",
    "baseline",
    "predicted",
    "delta"
  )

  missing_columns <- setdiff(
    required_columns,
    names(period_df)
  )

  if (length(missing_columns)) {
    stop(
      "`period_df` is missing required current EpiExposure column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (anyNA(period_df$var) ||
      any(!nzchar(as.character(period_df$var)))) {
    stop(
      "`period_df$var` must contain non-missing exposure names.",
      call. = FALSE
    )
  }

  if (anyNA(period_df$period) ||
      any(!nzchar(as.character(period_df$period)))) {
    stop(
      "`period_df$period` must contain non-missing, non-empty labels.",
      call. = FALSE
    )
  }

  if (!is.numeric(period_df$value) ||
      anyNA(period_df$value) ||
      any(!is.finite(period_df$value))) {
    stop(
      "`period_df$value` must contain only finite numeric exposure values.",
      call. = FALSE
    )
  }

  numeric_effect_columns <- c(
    "eta",
    "effect",
    "baseline",
    "predicted",
    "delta"
  )

  for (column in numeric_effect_columns) {
    x <- period_df[[column]]

    if (!is.numeric(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "`period_df$",
        column,
        "` must contain only finite numeric values.",
        call. = FALSE
      )
    }
  }

  if (anyNA(period_df$scale) ||
      any(as.character(period_df$scale) != "period")) {
    stop(
      "Every row of `period_df` must have `scale = 'period'`.",
      call. = FALSE
    )
  }

  period_df$var <- as.character(
    period_df$var
  )
  period_df$period <- as.character(
    period_df$period
  )

  has_samples <- "sample" %in%
    names(period_df)

  if (isTRUE(uncertainty) &&
      !has_samples) {
    stop(
      "`period_df` contains summarized uncertain effects. Correct uncertainty ",
      "for period comparisons requires the joint sample-level output. Rerun ",
      "`summarise_effects(..., uncertainty = TRUE, output = 'samples')` and ",
      "pass that result to `compare_periods()`.",
      call. = FALSE
    )
  }

  if (!isTRUE(uncertainty) &&
      has_samples) {
    stop(
      "`period_df` contains a `sample` column but its uncertainty metadata are ",
      "FALSE. The effect object is internally inconsistent.",
      call. = FALSE
    )
  }

  if (!isTRUE(uncertainty) &&
      identical(output, "samples")) {
    stop(
      "`output = 'samples'` is available only when `period_df` contains ",
      "sample-level uncertainty draws.",
      call. = FALSE
    )
  }

  if (has_samples) {
    if (!is.numeric(period_df$sample) ||
        anyNA(period_df$sample) ||
        any(!is.finite(period_df$sample)) ||
        any(period_df$sample < 1) ||
        any(period_df$sample !=
            as.integer(period_df$sample))) {
      stop(
        "`period_df$sample` must contain positive finite integer draw IDs.",
        call. = FALSE
      )
    }

    period_df$sample <- as.integer(
      period_df$sample
    )
  }

  # ==========================================================================
  # VALIDATE effect AGAINST eta AND STORED EFFECT MEASURE
  # ==========================================================================

  expected_effect <- switch(
    effect_measure,
    linear = period_df$eta,
    exponentiated = exp(period_df$eta),
    percent = 100 * (
      exp(period_df$eta) - 1
    )
  )

  if (any(!is.finite(expected_effect))) {
    stop(
      "The stored `eta` values imply non-finite transformed effects.",
      call. = FALSE
    )
  }

  if (any(!nearly_equal(
    period_df$effect,
    expected_effect,
    tolerance = tolerance
  ))) {
    stop(
      "`period_df$effect` is inconsistent with `eta` and the stored ",
      "`effect_measure = '",
      effect_measure,
      "'`. Do not mix effect outputs produced under different effect measures.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # OPTIONAL PERIOD-DEFINITION CONTRACT
  # ==========================================================================

  period_metadata <- NULL

  if (!is.null(periods)) {
    if (!is.data.frame(periods) ||
        !nrow(periods)) {
      stop(
        "`periods` must be NULL or a non-empty data.frame.",
        call. = FALSE
      )
    }

    required_period_columns <- c(
      "period",
      "lag_start",
      "lag_end"
    )

    missing_period_columns <- setdiff(
      required_period_columns,
      names(periods)
    )

    if (length(missing_period_columns)) {
      stop(
        "`periods` is missing required column(s): ",
        paste(
          missing_period_columns,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    periods_local <- periods[
      ,
      required_period_columns,
      drop = FALSE
    ]

    if (anyNA(periods_local$period) ||
        any(!nzchar(
          as.character(
            periods_local$period
          )
        )) ||
        anyDuplicated(
          as.character(
            periods_local$period
          )
        )) {
      stop(
        "`periods$period` must contain unique non-empty labels.",
        call. = FALSE
      )
    }

    periods_local$period <-
      as.character(
        periods_local$period
      )

    for (column in c(
      "lag_start",
      "lag_end"
    )) {
      x <- periods_local[[column]]

      if (!is.numeric(x) ||
          anyNA(x) ||
          any(!is.finite(x)) ||
          any(x < 0) ||
          any(x != as.integer(x))) {
        stop(
          "`periods$",
          column,
          "` must contain non-negative finite integers.",
          call. = FALSE
        )
      }

      periods_local[[column]] <-
        as.integer(x)
    }

    if (any(
      periods_local$lag_start >
      periods_local$lag_end
    )) {
      stop(
        "Every row of `periods` must satisfy `lag_start <= lag_end`.",
        call. = FALSE
      )
    }

    if (any(
      periods_local$lag_start > max_lag_metadata |
      periods_local$lag_end > max_lag_metadata
    )) {
      stop(
        "Every supplied period must lie within the common fitted lag range ",
        "0:", max_lag_metadata, ".",
        call. = FALSE
      )
    }

    periods_max_lag <- attr(
      periods,
      "epiexposure_max_lag",
      exact = TRUE
    )

    if (!is.null(periods_max_lag)) {
      if (!valid_integer_scalar(periods_max_lag) ||
          as.integer(periods_max_lag) != max_lag_metadata) {
        stop(
          "`periods` maximum-lag metadata do not match the fitted period-effect ",
          "object. Expected max_lag = ", max_lag_metadata, ".",
          call. = FALSE
        )
      }
    }

    periods_n_lags <- attr(
      periods,
      "epiexposure_n_lags",
      exact = TRUE
    )

    if (!is.null(periods_n_lags)) {
      if (!valid_integer_scalar(periods_n_lags) ||
          as.integer(periods_n_lags) != history_length_metadata) {
        stop(
          "`periods` lag-count metadata do not match the inherited exact ",
          "history length. Expected ", history_length_metadata,
          " lags (`max_lag + 1`).",
          call. = FALSE
        )
      }
    }

    missing_definitions <- setdiff(
      unique(period_df$period),
      periods_local$period
    )

    if (length(missing_definitions)) {
      stop(
        "`periods` does not define period label(s) present in `period_df`: ",
        paste(
          missing_definitions,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    period_metadata <- periods_local
    period_metadata$n_lags <-
      period_metadata$lag_end -
      period_metadata$lag_start + 1L

    period_order <- period_metadata$period

  } else {
    period_order <- unique(
      period_df$period
    )
  }

  # ==========================================================================
  # RESOLVE REFERENCE PERIOD
  # ==========================================================================

  if (is.numeric(period_ref)) {
    if (!valid_integer_scalar(period_ref) ||
        period_ref < 1L) {
      stop(
        "Numeric `period_ref` must be one positive integer period index.",
        call. = FALSE
      )
    }

    period_index <- as.integer(
      period_ref
    )

    if (period_index >
        length(period_order)) {
      stop(
        "Numeric `period_ref = ",
        period_index,
        "` exceeds the number of available ordered periods (",
        length(period_order),
        ").",
        call. = FALSE
      )
    }

    reference_period <-
      period_order[
        period_index
      ]

  } else if (valid_scalar_string(
    period_ref
  )) {
    reference_period <-
      as.character(period_ref)

  } else {
    stop(
      "`period_ref` must be one positive integer index or one exact character ",
      "period label.",
      call. = FALSE
    )
  }

  if (!reference_period %in%
      unique(period_df$period)) {
    stop(
      "Reference period '",
      reference_period,
      "' was not found in `period_df$period`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # STRICT UNIQUENESS AND SAMPLE COMPLETENESS
  # ==========================================================================

  identity_columns <- c(
    "var",
    "value",
    "period"
  )

  if (has_samples) {
    identity_columns <- c(
      "sample",
      identity_columns
    )
  }

  if (anyDuplicated(
    period_df[
      ,
      identity_columns,
      drop = FALSE
    ]
  )) {
    stop(
      "`period_df` contains duplicate rows for the same ",
      paste(
        identity_columns,
        collapse = " + "
      ),
      " combination.",
      call. = FALSE
    )
  }

  if (has_samples) {
    expected_samples <- sort(
      unique(period_df$sample)
    )
    n_expected_samples <- length(
      expected_samples
    )

    group_key <- interaction(
      period_df$var,
      period_df$value,
      period_df$period,
      drop = TRUE,
      lex.order = TRUE
    )

    sample_sets <- split(
      period_df$sample,
      group_key
    )

    incomplete <- vapply(
      sample_sets,
      function(x) {
        !identical(
          sort(unique(x)),
          expected_samples
        )
      },
      logical(1)
    )

    if (any(incomplete)) {
      stop(
        "At least one variable/value/period combination is missing one or more ",
        "parameter draws. Period uncertainty cannot be compared on different ",
        "sample supports.",
        call. = FALSE
      )
    }

    stored_n_samples <- attr(
      period_df,
      "epiexposure_n_samples",
      exact = TRUE
    )

    if (!is.null(stored_n_samples) &&
        (!valid_integer_scalar(
          stored_n_samples
        ) ||
        as.integer(stored_n_samples) !=
        n_expected_samples)) {
      stop(
        "`period_df` sample IDs disagree with stored `epiexposure_n_samples`.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # MATCH EACH ROW TO THE REFERENCE PERIOD
  # ==========================================================================

  match_columns <- c(
    "var",
    "value"
  )

  if (has_samples) {
    match_columns <- c(
      "sample",
      match_columns
    )
  }

  reference_rows <- period_df[
    period_df$period ==
      reference_period,
    c(
      match_columns,
      "eta",
      "effect",
      "baseline",
      "predicted",
      "delta"
    ),
    drop = FALSE
  ]

  if (!nrow(reference_rows)) {
    stop(
      "No rows were found for reference period '",
      reference_period,
      "'.",
      call. = FALSE
    )
  }

  if (anyDuplicated(
    reference_rows[
      ,
      match_columns,
      drop = FALSE
    ]
  )) {
    stop(
      "Reference period '",
      reference_period,
      "' is not unique within the comparison keys.",
      call. = FALSE
    )
  }

  names(reference_rows)[
    names(reference_rows) == "eta"
  ] <- "ref_eta"
  names(reference_rows)[
    names(reference_rows) == "effect"
  ] <- "ref_effect"
  names(reference_rows)[
    names(reference_rows) == "baseline"
  ] <- "ref_baseline"
  names(reference_rows)[
    names(reference_rows) == "predicted"
  ] <- "ref_predicted"
  names(reference_rows)[
    names(reference_rows) == "delta"
  ] <- "ref_delta"

  working <- period_df
  working$.epix_row_id <- seq_len(
    nrow(working)
  )

  joined <- merge(
    working,
    reference_rows,
    by = match_columns,
    all.x = TRUE,
    sort = FALSE
  )

  joined <- joined[
    order(joined$.epix_row_id),
    ,
    drop = FALSE
  ]
  rownames(joined) <- NULL

  reference_value_columns <- c(
    "ref_eta",
    "ref_effect",
    "ref_baseline",
    "ref_predicted",
    "ref_delta"
  )

  if (anyNA(
    joined[
      ,
      reference_value_columns,
      drop = FALSE
    ]
  )) {
    stop(
      "The selected reference period is missing for at least one variable, ",
      "exposure value",
      if (has_samples) ", or parameter draw" else "",
      ". Every comparison row must have a matched reference-period row.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # BASELINE CONSISTENCY
  # ==========================================================================

  baseline_equal <- nearly_equal(
    joined$baseline,
    joined$ref_baseline,
    tolerance = tolerance
  )

  if (any(!baseline_equal)) {
    stop(
      "Compared periods do not share the same joint response-scale baseline. ",
      "This usually means that outputs with different exposure references or ",
      "different fitted-model contracts were combined before calling ",
      "`compare_periods()`.",
      call. = FALSE
    )
  }

  # Under a common baseline:
  # (predicted_j - baseline) - (predicted_r - baseline)
  # == predicted_j - predicted_r.
  target_delta_identity <-
    (joined$predicted -
       joined$ref_predicted) -
    (joined$delta -
       joined$ref_delta)

  if (any(
    abs(target_delta_identity) >
    tolerance * pmax(
      1,
      abs(joined$predicted),
      abs(joined$ref_predicted),
      abs(joined$delta),
      abs(joined$ref_delta)
    )
  )) {
    stop(
      "`predicted` and `delta` columns are inconsistent with the shared ",
      "baseline contract.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # DRAW-BY-DRAW / DETERMINISTIC COMPARISONS
  # ==========================================================================

  joined$reference_period <-
    reference_period

  joined$eta_diff <-
    joined$eta -
    joined$ref_eta

  joined$diff <-
    joined$effect -
    joined$ref_effect

  joined$predicted_diff <-
    joined$predicted -
    joined$ref_predicted

  joined$delta_diff <-
    joined$delta -
    joined$ref_delta

  multiplicative_ratio_defined <-
    link_name %in% c(
      "log",
      "logit"
    )

  if (multiplicative_ratio_defined) {
    joined$ratio <- exp(
      joined$eta_diff
    )

    if (any(!is.finite(joined$ratio)) ||
        any(joined$ratio <= 0)) {
      stop(
        "Period-relative multiplicative ratios are non-finite. The difference ",
        "in linear-predictor effects is too extreme for stable exponentiation.",
        call. = FALSE
      )
    }

    joined$ratio_percent <-
      100 * (
        joined$ratio - 1
      )

  } else {
    joined$ratio <- NA_real_
    joined$ratio_percent <- NA_real_
  }

  joined$.epix_row_id <- NULL

  # ==========================================================================
  # ADD OPTIONAL PERIOD WIDTHS
  # ==========================================================================

  if (!is.null(period_metadata)) {
    period_match <- match(
      joined$period,
      period_metadata$period
    )

    joined$lag_start <-
      period_metadata$lag_start[
        period_match
      ]
    joined$lag_end <-
      period_metadata$lag_end[
        period_match
      ]
    joined$n_lags <-
      period_metadata$n_lags[
        period_match
      ]
  }

  # ==========================================================================
  # ORDER DETERMINISTIC / SAMPLE OUTPUT
  # ==========================================================================

  period_rank <- match(
    joined$period,
    period_order
  )

  order_args <- list(
    joined$var,
    joined$value,
    period_rank
  )

  if (has_samples) {
    order_args <- c(
      list(joined$sample),
      order_args
    )
  }

  joined <- joined[
    do.call(
      order,
      order_args
    ),
    ,
    drop = FALSE
  ]
  rownames(joined) <- NULL

  # ==========================================================================
  # SAMPLE-LEVEL RETURN
  # ==========================================================================

  if (has_samples &&
      identical(output, "samples")) {

    attr(
      joined,
      "epiexposure_reference_period"
    ) <- reference_period

    attr(
      joined,
      "epiexposure_effect_measure"
    ) <- effect_measure

    attr(
      joined,
      "epiexposure_link"
    ) <- link_name

    attr(
      joined,
      "epiexposure_prediction_level"
    ) <- prediction_level

    attr(
      joined,
      "epiexposure_prediction_estimand"
    ) <- prediction_estimand

    attr(
      joined,
      "epiexposure_uncertainty"
    ) <- TRUE

    attr(
      joined,
      "epiexposure_comparison_output"
    ) <- "samples"

    attr(
      joined,
      "epiexposure_ratio_defined"
    ) <- multiplicative_ratio_defined

    attr(
      joined,
      "epiexposure_ratio_interpretation"
    ) <- if (
      identical(link_name, "log")
    ) {
      "ratio_of_period_specific_response_ratios"
    } else if (
      identical(link_name, "logit")
    ) {
      "ratio_of_period_specific_odds_ratios"
    } else {
      "not_defined_for_this_link"
    }

    attr(
      joined,
      "epiexposure_period_comparison_contract"
    ) <- "same_value_same_draw_reference_period"

    attr(
      joined,
      "epiexposure_max_lag"
    ) <- max_lag_metadata

    attr(
      joined,
      "epiexposure_history_length"
    ) <- history_length_metadata

    attr(
      joined,
      "epiexposure_history_contract"
    ) <- history_contract_metadata

    return(joined)
  }

  # ==========================================================================
  # DETERMINISTIC SUMMARY RETURN
  # ==========================================================================

  if (!has_samples) {
    attr(
      joined,
      "epiexposure_reference_period"
    ) <- reference_period

    attr(
      joined,
      "epiexposure_effect_measure"
    ) <- effect_measure

    attr(
      joined,
      "epiexposure_link"
    ) <- link_name

    attr(
      joined,
      "epiexposure_prediction_level"
    ) <- prediction_level

    attr(
      joined,
      "epiexposure_prediction_estimand"
    ) <- prediction_estimand

    attr(
      joined,
      "epiexposure_uncertainty"
    ) <- FALSE

    attr(
      joined,
      "epiexposure_comparison_output"
    ) <- "summary"

    attr(
      joined,
      "epiexposure_ratio_defined"
    ) <- multiplicative_ratio_defined

    attr(
      joined,
      "epiexposure_ratio_interpretation"
    ) <- if (
      identical(link_name, "log")
    ) {
      "ratio_of_period_specific_response_ratios"
    } else if (
      identical(link_name, "logit")
    ) {
      "ratio_of_period_specific_odds_ratios"
    } else {
      "not_defined_for_this_link"
    }

    attr(
      joined,
      "epiexposure_period_comparison_contract"
    ) <- "same_value_reference_period"

    attr(
      joined,
      "epiexposure_max_lag"
    ) <- max_lag_metadata

    attr(
      joined,
      "epiexposure_history_length"
    ) <- history_length_metadata

    attr(
      joined,
      "epiexposure_history_contract"
    ) <- history_contract_metadata

    return(joined)
  }

  # ==========================================================================
  # SUMMARISE DRAW-BY-DRAW COMPARISONS
  # ==========================================================================

  if (is.null(interval_probs)) {
    interval_probs <- attr(
      period_df,
      "epiexposure_interval_probs",
      exact = TRUE
    )

    if (is.null(interval_probs)) {
      interval_probs <- c(
        0.025,
        0.975
      )
    }
  }

  interval_probs <-
    validate_interval_probs(
      interval_probs
    )

  summary_keys <- c(
    "var",
    "value",
    "period"
  )

  group_id <- interaction(
    joined$var,
    joined$value,
    joined$period,
    drop = TRUE,
    lex.order = TRUE
  )

  rows_by_group <- split(
    seq_len(nrow(joined)),
    group_id
  )

  quantities_to_summarize <- c(
    "eta",
    "ref_eta",
    "eta_diff",
    "effect",
    "ref_effect",
    "diff",
    "baseline",
    "predicted",
    "ref_predicted",
    "predicted_diff",
    "delta",
    "ref_delta",
    "delta_diff"
  )

  if (multiplicative_ratio_defined) {
    quantities_to_summarize <- c(
      quantities_to_summarize,
      "ratio",
      "ratio_percent"
    )
  }

  summary_rows <- vector(
    "list",
    length(rows_by_group)
  )

  row_index <- 1L

  for (idx in rows_by_group) {
    current <- joined[
      idx,
      ,
      drop = FALSE
    ]

    base <- current[
      1L,
      summary_keys,
      drop = FALSE
    ]

    base$scale <- "period"
    base$reference_period <-
      reference_period

    if (!is.null(period_metadata)) {
      base$lag_start <-
        current$lag_start[1L]
      base$lag_end <-
        current$lag_end[1L]
      base$n_lags <-
        current$n_lags[1L]
    }

    for (quantity in quantities_to_summarize) {
      summary_values <- summarize_vector(
        current[[quantity]],
        probs = interval_probs
      )

      base[[quantity]] <-
        unname(
          summary_values["value"]
        )
      base[[paste0(
        quantity,
        "_sd"
      )]] <-
        unname(
          summary_values["sd"]
        )
      base[[paste0(
        quantity,
        "_lower"
      )]] <-
        unname(
          summary_values["lower"]
        )
      base[[paste0(
        quantity,
        "_upper"
      )]] <-
        unname(
          summary_values["upper"]
        )
    }

    if (!multiplicative_ratio_defined) {
      base$ratio <- NA_real_
      base$ratio_sd <- NA_real_
      base$ratio_lower <- NA_real_
      base$ratio_upper <- NA_real_

      base$ratio_percent <- NA_real_
      base$ratio_percent_sd <- NA_real_
      base$ratio_percent_lower <- NA_real_
      base$ratio_percent_upper <- NA_real_
    }

    summary_rows[[row_index]] <-
      base
    row_index <- row_index + 1L
  }

  out <- do.call(
    rbind,
    summary_rows
  )
  rownames(out) <- NULL

  out_period_rank <- match(
    out$period,
    period_order
  )

  out <- out[
    order(
      out$var,
      out$value,
      out_period_rank
    ),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL

  attr(
    out,
    "epiexposure_reference_period"
  ) <- reference_period

  attr(
    out,
    "epiexposure_effect_measure"
  ) <- effect_measure

  attr(
    out,
    "epiexposure_link"
  ) <- link_name

  attr(
    out,
    "epiexposure_prediction_level"
  ) <- prediction_level

  attr(
    out,
    "epiexposure_prediction_estimand"
  ) <- prediction_estimand

  attr(
    out,
    "epiexposure_uncertainty"
  ) <- TRUE

  attr(
    out,
    "epiexposure_n_samples"
  ) <- length(
    unique(period_df$sample)
  )

  attr(
    out,
    "epiexposure_interval_probs"
  ) <- interval_probs

  attr(
    out,
    "epiexposure_summary_center"
  ) <- "median"

  attr(
    out,
    "epiexposure_comparison_output"
  ) <- "summary"

  attr(
    out,
    "epiexposure_ratio_defined"
  ) <- multiplicative_ratio_defined

  attr(
    out,
    "epiexposure_ratio_interpretation"
  ) <- if (
    identical(link_name, "log")
  ) {
    "ratio_of_period_specific_response_ratios"
  } else if (
    identical(link_name, "logit")
  ) {
    "ratio_of_period_specific_odds_ratios"
  } else {
    "not_defined_for_this_link"
  }

  attr(
    out,
    "epiexposure_period_comparison_contract"
  ) <- "same_value_same_draw_reference_period"

  attr(
    out,
    "epiexposure_max_lag"
  ) <- max_lag_metadata

  attr(
    out,
    "epiexposure_history_length"
  ) <- history_length_metadata

  attr(
    out,
    "epiexposure_history_contract"
  ) <- history_contract_metadata

  out
}
