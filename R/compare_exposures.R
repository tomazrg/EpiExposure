#' Compare chronological exposure profiles
#'
#' Compares two or more exposure histories, simulated scenarios, or
#' epidemic-specific exposure profiles while preserving the EpiExposure
#' chronological profile convention.
#'
#' Each profile must be ordered from the **oldest/earliest** exposure
#' observation to the **most recent** observation. Consequently, for a profile
#' of length `n`, the final element corresponds to retrospective lag 0 and the
#' first element corresponds to positional lag `n - 1`.
#'
#' `compare_exposures()` is descriptive. It compares exposure histories
#' themselves and does not estimate DLNM effects, disease risks, or model
#' predictions.
#'
#' @param exposures Named list of numeric exposure profiles. Every element must
#'   represent the same exposure variable, in the same units, over the same
#'   chronological positions, and all profiles must have equal length.
#'
#'   Profiles must be stored from the earliest/oldest observation to the most
#'   recent observation. A single numeric vector is still accepted for backward
#'   compatibility, but at least two profiles or aggregated groups are required
#'   for an actual comparison.
#' @param exposure1 Legacy numeric exposure vector. It must be supplied together
#'   with `exposure2` when `exposures = NULL`.
#' @param exposure2 Legacy numeric exposure vector. It must be supplied together
#'   with `exposure1` when `exposures = NULL`.
#'
#'   If `exposures` is supplied, `exposure1` and `exposure2` must both be
#'   `NULL`; legacy arguments are not silently ignored.
#' @param time Optional numeric vector of chronological time positions. It must
#'   have the same length as every exposure profile, contain unique finite
#'   values, be strictly increasing, and be equally spaced.
#'
#'   Equal spacing is required because EpiExposure profiles represent complete
#'   regular exposure histories and several global metrics treat each temporal
#'   position equally.
#'
#'   If `NULL`, `0:(n - 1)` is used.
#' @param group Optional group membership for the supplied profiles. When
#'   provided, profiles are first aggregated within each group at every time
#'   position using `agg_fun`, and pairwise comparisons are made between the
#'   resulting group-level profiles.
#'
#'   An unnamed vector is matched by profile position. A named vector is safer:
#'   its names must match `names(exposures)` exactly and are used to align group
#'   membership to profiles before aggregation.
#' @param agg_fun Character. `"median"` or `"mean"`. Used only when `group` is
#'   supplied. The default is `"median"`.
#' @param mode Character. `"global"` or `"timewise"`.
#'
#'   \describe{
#'     \item{`"global"`}{Returns one row per pairwise profile comparison with
#'       distance, magnitude, correlation, peak-timing, center-of-mass, and
#'       high-exposure overlap summaries.}
#'     \item{`"timewise"`}{Returns one row per chronological position and
#'       pairwise comparison, retaining the complete exposure histories.}
#'   }
#' @param cor_method Character correlation method used in `mode = "global"`.
#'   Available options are:
#'
#'   - `"pearson"` (default): linear association between the original exposure
#'     values;
#'   - `"spearman"`: monotonic association calculated from exposure ranks;
#'   - `"kendall"`: Kendall's rank correlation (tau), also measuring monotonic
#'     association and often useful for small samples or many non-linear
#'     monotonic relationships.
#'
#'   The selected method changes only the `corr` summary. All distance,
#'   magnitude, timing, center-of-mass, and overlap metrics are unchanged.
#'
#'   Correlation is returned as `NA` when either profile is constant because
#'   none of the supported correlation coefficients is defined in that case.
#' @param q Finite numeric scalar between 0 and 1. Each profile's own `q`
#'   quantile defines its high-exposure positions for the Jaccard-style
#'   `overlap_above_q` metric. The default is `0.8`.
#' @param eps Positive finite numeric tolerance. It is used only to decide when
#'   a time-wise ratio denominator or a non-negative profile's total
#'   center-of-mass weight is numerically indistinguishable from zero.
#'
#'   `eps` is **not added to exposure values** and therefore never modifies a
#'   profile before comparison.
#'
#' @return A data frame.
#'
#'   With `mode = "timewise"`, the output contains:
#'
#'   - `exposure1`, `exposure2`: compared profile/group names;
#'   - `position`: chronological position `1, ..., n`;
#'   - `time`: chronological time supplied by the user;
#'   - `lag`: positional retrospective lag, with lag 0 at the most recent
#'     observation;
#'   - `value1`, `value2`;
#'   - `diff = value1 - value2`;
#'   - `abs_diff = abs(diff)`;
#'   - `ratio = value1 / value2` when `abs(value2) > eps`, otherwise `NA`;
#'   - `ratio_defined`: whether that numerical ratio was computed.
#'
#'   With `mode = "global"`, the output contains one row per pair and includes:
#'
#'   - discrete distances `L1`, `L2`, and `Linf`;
#'   - `MAE` and `RMSE`, which are length-normalized profile discrepancies;
#'   - individual and between-profile sums, means, and SDs;
#'   - the selected Pearson, Spearman, or Kendall correlation when both
#'     profiles vary;
#'   - peak values and peak timing;
#'   - explicit information about tied maxima;
#'   - temporal centers of mass when the corresponding profile is non-negative
#'     and has positive total mass;
#'   - profile-specific high-exposure thresholds and Jaccard overlap.
#'
#'   Output attributes record chronological order, common time step, comparison
#'   direction, grouping information, and the ratio/center-of-mass contracts.
#'
#' @details
#' ## Comparison direction
#'
#' Pairwise differences are directional:
#'
#' \deqn{difference = exposure1 - exposure2.}
#'
#' Therefore positive `diff`, `diff_sum`, or `diff_mean` values indicate larger
#' exposure in the first named profile than in the second.
#'
#' `timing_shift_peak` and `center_of_mass_shift` follow the same direction:
#'
#' \deqn{time_1 - time_2.}
#'
#' Positive values therefore mean that the corresponding timing summary occurs
#' later in `exposure1`.
#'
#' ## Discrete distance metrics
#'
#' For paired profile values \eqn{x_{1t}} and \eqn{x_{2t}},
#'
#' \deqn{L1 = \sum_t |x_{1t} - x_{2t}|,}
#'
#' \deqn{L2 = \sqrt{\sum_t (x_{1t} - x_{2t})^2},}
#'
#' and
#'
#' \deqn{Linf = \max_t |x_{1t} - x_{2t}|.}
#'
#' `L1` and `L2` increase with the number of represented time positions. `MAE`
#' and `RMSE` are therefore also returned when comparisons across equal-unit
#' histories of different lengths are needed.
#'
#' These distances retain the units/scaling of the exposure profile and should
#' only be compared across profiles of the same exposure variable and units.
#'
#' ## Time-wise ratios
#'
#' Earlier EpiExposure code calculated
#'
#' \deqn{x_1 / (x_2 + eps),}
#'
#' which changed every denominator, including perfectly valid non-zero values.
#' The current implementation never alters observed exposure values. It returns
#'
#' \deqn{x_1/x_2}
#'
#' only when `abs(x2) > eps`; otherwise the ratio is `NA`.
#'
#' A numerical ratio is not automatically scientifically meaningful. It should
#' only be interpreted as a relative exposure contrast when the exposure is
#' measured on a ratio scale with a meaningful zero. For example, ratios of
#' temperatures expressed in degrees Celsius generally should not be interpreted
#' as relative temperature effects.
#'
#' ## Correlation method
#'
#' `cor_method` controls only the global `corr` column:
#'
#' - Pearson measures linear co-variation of the original profile values;
#' - Spearman is Pearson correlation of the ranks and measures monotonic
#'   association;
#' - Kendall returns Kendall's tau, based on concordant and discordant pairs.
#'
#' Correlation describes similarity in temporal pattern, not similarity in
#' absolute exposure magnitude. For example, two profiles can have correlation
#' close to 1 while one is systematically much larger than the other. Distance
#' and magnitude metrics should therefore be interpreted alongside `corr`.
#'
#' No correlation is causal and no p-value is computed by this descriptive
#' function.
#'
#' ## Peak timing and ties
#'
#' `peak_time1` and `peak_time2` retain the first chronological occurrence of
#' the maximum for backward compatibility. The function additionally returns
#' first and last peak times, the midpoint of their temporal range, the number
#' of maximum positions, and indicators of tied maxima. Thus a plateau or
#' repeated equal maximum is never hidden by `which.max()`.
#'
#' `timing_shift_peak` uses the first peak times for backward compatibility.
#' `timing_shift_peak_midpoint` compares the midpoints of the full peak-time
#' ranges and is usually more informative when maxima are tied.
#'
#' ## Temporal center of mass
#'
#' For a non-negative profile with positive total exposure, the temporal center
#' of mass is
#'
#' \deqn{COM = \frac{\sum_t time_t x_t}{\sum_t x_t}.}
#'
#' Earlier versions silently replaced negative exposure values by zero before
#' computing COM. That changes the supplied profile and is no longer done.
#'
#' If a profile contains any negative value, or its non-negative total mass is
#' not greater than `eps`, its COM is returned as `NA` and `com_defined` is
#' `FALSE`. Even for non-negative data, COM is scientifically interpretable only
#' when zero and exposure magnitude have a meaningful weighting interpretation.
#'
#' ## High-exposure overlap
#'
#' Each profile obtains its own threshold:
#'
#' \deqn{Q_j(q) = quantile(x_j, q).}
#'
#' High-exposure positions satisfy \eqn{x_{jt} >= Q_j(q)}. The returned
#' `overlap_above_q` is the Jaccard index
#'
#' \deqn{
#'   \frac{|H_1 \cap H_2|}{|H_1 \cup H_2|}.
#' }
#'
#' This measures whether the two profiles experience their **own relatively high
#' exposure** at the same times. It is not an overlap above one common absolute
#' exposure threshold. Ties at the quantile threshold can make the number of
#' high positions larger than exactly `(1 - q) * n`; counts and thresholds are
#' returned explicitly.
#'
#' ## Group aggregation
#'
#' With `group`, aggregation occurs independently at every chronological
#' position before any pairwise metric is calculated. The resulting comparison
#' therefore describes representative group-level profiles, not a distribution
#' of pairwise epidemic-level differences. No uncertainty interval or
#' inferential test is produced by this function.
#'
#' @export
compare_exposures <- function(
    exposures = NULL,
    exposure1 = NULL,
    exposure2 = NULL,
    time = NULL,
    group = NULL,
    agg_fun = c("median", "mean"),
    mode = c("global", "timewise"),
    cor_method = c("pearson", "spearman", "kendall"),
    q = 0.8,
    eps = 1e-12
) {

  # ==========================================================================
  # ARGUMENT MATCHING AND SMALL VALIDATORS
  # ==========================================================================

  mode <- match.arg(mode)
  agg_fun <- match.arg(agg_fun)
  cor_method <- match.arg(cor_method)

  valid_scalar_numeric <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x)
  }

  nearly_equal <- function(a, b, tolerance) {
    abs(a - b) <= tolerance *
      max(1, abs(a), abs(b))
  }

  if (!valid_scalar_numeric(q) ||
      q < 0 ||
      q > 1) {
    stop(
      "`q` must be one finite numeric value between 0 and 1.",
      call. = FALSE
    )
  }

  if (!valid_scalar_numeric(eps) ||
      eps <= 0) {
    stop(
      "`eps` must be one positive finite numeric value.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # BACKWARD-COMPATIBILITY INPUT
  # ==========================================================================

  if (is.null(exposures)) {
    if (is.null(exposure1) ||
        is.null(exposure2)) {
      stop(
        "Provide either `exposures` or both legacy arguments `exposure1` and ",
        "`exposure2`.",
        call. = FALSE
      )
    }

    exposures <- list(
      exposure1 = exposure1,
      exposure2 = exposure2
    )

  } else {
    if (!is.null(exposure1) ||
        !is.null(exposure2)) {
      stop(
        "When `exposures` is supplied, legacy arguments `exposure1` and ",
        "`exposure2` must both be NULL.",
        call. = FALSE
      )
    }

    if (is.numeric(exposures)) {
      exposures <- list(
        exposure = exposures
      )
    }
  }

  # ==========================================================================
  # VALIDATE EXPOSURE PROFILES
  # ==========================================================================

  if (!is.list(exposures) ||
      !length(exposures)) {
    stop(
      "`exposures` must be a non-empty named list of numeric vectors.",
      call. = FALSE
    )
  }

  exposure_names <- names(exposures)

  if (is.null(exposure_names) ||
      anyNA(exposure_names) ||
      any(!nzchar(exposure_names)) ||
      anyDuplicated(exposure_names)) {
    stop(
      "All exposure profiles must have unique, non-empty names.",
      call. = FALSE
    )
  }

  profile_lengths <- vapply(
    exposures,
    length,
    integer(1)
  )

  if (any(profile_lengths < 2L)) {
    stop(
      "Every exposure profile must contain at least two chronological ",
      "observations.",
      call. = FALSE
    )
  }

  if (length(unique(profile_lengths)) != 1L) {
    stop(
      "All exposure profiles must have the same length and represent the same ",
      "chronological positions.",
      call. = FALSE
    )
  }

  n_time <- unname(
    profile_lengths[1L]
  )

  for (profile_name in exposure_names) {
    x <- exposures[[profile_name]]

    if (!is.numeric(x)) {
      stop(
        "Exposure profile '",
        profile_name,
        "' must be numeric.",
        call. = FALSE
      )
    }

    if (anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "Exposure profile '",
        profile_name,
        "' must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    exposures[[profile_name]] <-
      as.numeric(x)
  }

  # ==========================================================================
  # CHRONOLOGICAL TIME CONTRACT
  # ==========================================================================

  if (is.null(time)) {
    time <- seq.int(
      from = 0,
      to = n_time - 1L
    )
  }

  if (!is.numeric(time) ||
      length(time) != n_time ||
      anyNA(time) ||
      any(!is.finite(time))) {
    stop(
      "`time` must contain one finite numeric value for every chronological ",
      "profile position.",
      call. = FALSE
    )
  }

  time <- as.numeric(time)

  if (anyDuplicated(time)) {
    stop(
      "`time` must contain unique chronological positions.",
      call. = FALSE
    )
  }

  if (is.unsorted(
    time,
    strictly = TRUE
  )) {
    stop(
      "`time` must increase strictly from the earliest/oldest exposure ",
      "observation to the most recent observation.",
      call. = FALSE
    )
  }

  time_differences <- diff(time)

  if (any(time_differences <= 0) ||
      any(!is.finite(time_differences))) {
    stop(
      "`time` must increase by a positive finite interval.",
      call. = FALSE
    )
  }

  time_step <- time_differences[1L]
  time_tolerance <- sqrt(
    .Machine$double.eps
  )

  if (length(time_differences) > 1L &&
      any(
        abs(
          time_differences -
          time_step
        ) >
        time_tolerance *
        pmax(
          1,
          abs(time_differences),
          abs(time_step)
        )
      )) {
    stop(
      "`time` must be equally spaced. EpiExposure exposure profiles represent ",
      "complete regular temporal histories.",
      call. = FALSE
    )
  }

  lag_position <- rev(
    seq.int(
      from = 0L,
      to = n_time - 1L
    )
  )

  chronological_position <-
    seq_len(n_time)

  # ==========================================================================
  # OPTIONAL GROUP ALIGNMENT AND AGGREGATION
  # ==========================================================================

  original_n_profiles <-
    length(exposures)

  original_profile_names <-
    names(exposures)

  grouped <- !is.null(group)
  group_sizes <- NULL
  group_members <- NULL

  if (grouped) {
    if (is.factor(group)) {
      group <- as.character(group)
    }

    group_names <- names(group)

    if (!is.null(group_names)) {
      if (anyNA(group_names) ||
          any(!nzchar(group_names)) ||
          anyDuplicated(group_names) ||
          !setequal(
            group_names,
            original_profile_names
          ) ||
          length(group) !=
          original_n_profiles) {
        stop(
          "When `group` is named, its names must match `names(exposures)` ",
          "exactly and uniquely.",
          call. = FALSE
        )
      }

      group <- group[
        original_profile_names
      ]
    } else if (
      length(group) !=
      original_n_profiles
    ) {
      stop(
        "Unnamed `group` must have the same length as `exposures`.",
        call. = FALSE
      )
    }

    if (anyNA(group)) {
      stop(
        "`group` cannot contain missing values.",
        call. = FALSE
      )
    }

    group <- as.character(group)

    if (any(!nzchar(group))) {
      stop(
        "`group` cannot contain empty labels.",
        call. = FALSE
      )
    }

    group_levels <- unique(group)

    if (length(group_levels) < 2L) {
      stop(
        "`group` must define at least two groups for pairwise comparison.",
        call. = FALSE
      )
    }

    exposure_matrix <- do.call(
      rbind,
      exposures
    )

    if (!is.matrix(exposure_matrix) ||
        nrow(exposure_matrix) !=
        original_n_profiles ||
        ncol(exposure_matrix) !=
        n_time) {
      stop(
        "Internal error while assembling exposure profiles for group ",
        "aggregation.",
        call. = FALSE
      )
    }

    aggregate_function <- if (
      identical(
        agg_fun,
        "median"
      )
    ) {
      stats::median
    } else {
      base::mean
    }

    grouped_exposures <-
      vector(
        "list",
        length(group_levels)
      )
    names(grouped_exposures) <-
      group_levels

    group_sizes <-
      stats::setNames(
        integer(
          length(group_levels)
        ),
        group_levels
      )

    group_members <-
      stats::setNames(
        vector(
          "list",
          length(group_levels)
        ),
        group_levels
      )

    for (group_level in group_levels) {
      selected <- group ==
        group_level

      group_sizes[[group_level]] <-
        sum(selected)

      group_members[[group_level]] <-
        original_profile_names[
          selected
        ]

      grouped_exposures[[group_level]] <-
        apply(
          exposure_matrix[
            selected,
            ,
            drop = FALSE
          ],
          MARGIN = 2,
          FUN = aggregate_function
        )

      if (anyNA(
        grouped_exposures[[
          group_level
        ]]
      ) ||
      any(!is.finite(
        grouped_exposures[[
          group_level
        ]]
      ))) {
        stop(
          "Group aggregation produced non-finite values for group '",
          group_level,
          "'.",
          call. = FALSE
        )
      }
    }

    exposures <- grouped_exposures
  }

  if (length(exposures) < 2L) {
    stop(
      "At least two exposure profiles or aggregated groups are required for ",
      "comparison.",
      call. = FALSE
    )
  }

  comparison_names <- names(exposures)

  pairwise_combinations <-
    utils::combn(
      comparison_names,
      2,
      simplify = FALSE
    )

  # ==========================================================================
  # OUTPUT METADATA HELPER
  # ==========================================================================

  attach_metadata <- function(out) {
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
      "epiexposure_time_step"
    ) <- as.numeric(time_step)

    attr(
      out,
      "epiexposure_n_time"
    ) <- n_time

    attr(
      out,
      "epiexposure_compare_mode"
    ) <- mode

    attr(
      out,
      "epiexposure_correlation_method"
    ) <- cor_method

    attr(
      out,
      "epiexposure_difference_direction"
    ) <- "exposure1_minus_exposure2"

    attr(
      out,
      "epiexposure_ratio_contract"
    ) <- "x1_over_x2_when_abs_x2_gt_eps"

    attr(
      out,
      "epiexposure_center_of_mass_contract"
    ) <- "nonnegative_profiles_only_no_clipping"

    attr(
      out,
      "epiexposure_high_overlap_contract"
    ) <- "profile_specific_quantile_jaccard"

    attr(
      out,
      "epiexposure_grouped"
    ) <- grouped

    attr(
      out,
      "epiexposure_aggregation"
    ) <- if (grouped) {
      agg_fun
    } else {
      NULL
    }

    attr(
      out,
      "epiexposure_group_sizes"
    ) <- group_sizes

    attr(
      out,
      "epiexposure_group_members"
    ) <- group_members

    attr(
      out,
      "epiexposure_descriptive_only"
    ) <- TRUE

    out
  }

  # ==========================================================================
  # TIME-WISE MODE
  # ==========================================================================

  if (identical(
    mode,
    "timewise"
  )) {

    timewise_rows <- lapply(
      pairwise_combinations,
      function(comparison) {

        name1 <- comparison[1L]
        name2 <- comparison[2L]

        x1 <- exposures[[name1]]
        x2 <- exposures[[name2]]

        denominator_defined <-
          abs(x2) > eps

        ratio <- rep(
          NA_real_,
          n_time
        )

        ratio[
          denominator_defined
        ] <-
          x1[
            denominator_defined
          ] /
          x2[
            denominator_defined
          ]

        if (any(
          !is.finite(
            ratio[
              denominator_defined
            ]
          )
        )) {
          stop(
            "A finite time-wise ratio could not be calculated for comparison '",
            name1,
            "' versus '",
            name2,
            "'.",
            call. = FALSE
          )
        }

        data.frame(
          exposure1 = name1,
          exposure2 = name2,
          position = chronological_position,
          time = time,
          lag = lag_position,
          value1 = x1,
          value2 = x2,
          diff = x1 - x2,
          abs_diff = abs(
            x1 - x2
          ),
          ratio = ratio,
          ratio_defined =
            denominator_defined,
          stringsAsFactors = FALSE
        )
      }
    )

    out <- do.call(
      rbind,
      timewise_rows
    )
    rownames(out) <- NULL

    return(
      attach_metadata(out)
    )
  }

  # ==========================================================================
  # GLOBAL-METRIC HELPERS
  # ==========================================================================

  peak_summary <- function(x) {
    peak_value <- max(x)
    peak_indices <- which(
      x == peak_value
    )

    if (!length(peak_indices)) {
      stop(
        "Internal error while identifying profile maximum.",
        call. = FALSE
      )
    }

    peak_times <- time[
      peak_indices
    ]

    peak_lags <- lag_position[
      peak_indices
    ]

    list(
      value = peak_value,
      n = length(peak_indices),
      tied = length(peak_indices) > 1L,
      first_time = peak_times[1L],
      last_time = peak_times[
        length(peak_times)
      ],
      midpoint_time = mean(
        range(peak_times)
      ),
      first_lag = peak_lags[1L],
      last_lag = peak_lags[
        length(peak_lags)
      ],
      midpoint_lag = mean(
        range(peak_lags)
      )
    )
  }

  center_of_mass <- function(x) {
    if (any(x < 0)) {
      return(
        list(
          value = NA_real_,
          defined = FALSE,
          reason = "negative_values"
        )
      )
    }

    total_mass <- sum(x)

    if (!is.finite(total_mass) ||
        total_mass <= eps) {
      return(
        list(
          value = NA_real_,
          defined = FALSE,
          reason = "zero_or_near_zero_total_mass"
        )
      )
    }

    com_value <- sum(
      time * x
    ) / total_mass

    if (!is.finite(com_value)) {
      stop(
        "Internal center-of-mass calculation produced a non-finite value.",
        call. = FALSE
      )
    }

    list(
      value = com_value,
      defined = TRUE,
      reason = "defined"
    )
  }

  high_overlap <- function(
    x1,
    x2
  ) {

    threshold1 <-
      as.numeric(
        stats::quantile(
          x1,
          probs = q,
          names = FALSE,
          type = 7
        )
      )

    threshold2 <-
      as.numeric(
        stats::quantile(
          x2,
          probs = q,
          names = FALSE,
          type = 7
        )
      )

    if (!is.finite(threshold1) ||
        !is.finite(threshold2)) {
      stop(
        "Could not calculate finite high-exposure thresholds.",
        call. = FALSE
      )
    }

    high1 <- x1 >=
      threshold1
    high2 <- x2 >=
      threshold2

    n_high1 <- sum(high1)
    n_high2 <- sum(high2)
    n_intersection <- sum(
      high1 & high2
    )
    n_union <- sum(
      high1 | high2
    )

    overlap <- if (
      n_union == 0L
    ) {
      NA_real_
    } else {
      n_intersection /
        n_union
    }

    list(
      threshold1 = threshold1,
      threshold2 = threshold2,
      n_high1 = n_high1,
      n_high2 = n_high2,
      n_intersection = n_intersection,
      n_union = n_union,
      overlap = overlap
    )
  }

  # ==========================================================================
  # GLOBAL MODE
  # ==========================================================================

  global_rows <- lapply(
    pairwise_combinations,
    function(comparison) {

      name1 <- comparison[1L]
      name2 <- comparison[2L]

      x1 <- exposures[[name1]]
      x2 <- exposures[[name2]]

      difference <- x1 - x2
      absolute_difference <-
        abs(difference)

      squared_difference <-
        difference^2

      # ----------------------------------------------------------------------
      # DISTANCE / DISCREPANCY
      # ----------------------------------------------------------------------

      L1 <- sum(
        absolute_difference
      )

      L2 <- sqrt(
        sum(
          squared_difference
        )
      )

      Linf <- max(
        absolute_difference
      )

      MAE <- mean(
        absolute_difference
      )

      RMSE <- sqrt(
        mean(
          squared_difference
        )
      )

      # ----------------------------------------------------------------------
      # PROFILE MAGNITUDE / VARIABILITY
      # ----------------------------------------------------------------------

      sum1 <- sum(x1)
      sum2 <- sum(x2)

      mean1 <- mean(x1)
      mean2 <- mean(x2)

      sd1 <- stats::sd(x1)
      sd2 <- stats::sd(x2)

      if (any(!is.finite(
        c(
          sum1,
          sum2,
          mean1,
          mean2,
          sd1,
          sd2
        )
      ))) {
        stop(
          "Internal profile summary produced a non-finite value.",
          call. = FALSE
        )
      }

      # ----------------------------------------------------------------------
      # CORRELATION
      # ----------------------------------------------------------------------

      # A correlation coefficient is undefined when either profile is constant.
      # Using uniqueness rather than SD makes the check equally appropriate for
      # Pearson, Spearman, and Kendall.
      varying1 <- length(unique(x1)) > 1L
      varying2 <- length(unique(x2)) > 1L

      corr <- if (
        !varying1 ||
        !varying2
      ) {
        NA_real_
      } else {
        as.numeric(
          stats::cor(
            x1,
            x2,
            method = cor_method
          )
        )
      }

      if (!is.na(corr) &&
          !is.finite(corr)) {
        stop(
          "Correlation calculation using method '",
          cor_method,
          "' returned a non-finite value.",
          call. = FALSE
        )
      }

      if (!is.na(corr) &&
          (corr < -1 - 1e-12 ||
           corr > 1 + 1e-12)) {
        stop(
          "Correlation calculation using method '",
          cor_method,
          "' returned a value outside [-1, 1].",
          call. = FALSE
        )
      }

      # Protect only against negligible floating-point excursions.
      if (!is.na(corr)) {
        corr <- min(
          max(corr, -1),
          1
        )
      }

      # ----------------------------------------------------------------------
      # PEAK TIMING
      # ----------------------------------------------------------------------

      peak1 <- peak_summary(x1)
      peak2 <- peak_summary(x2)

      timing_shift_peak <-
        peak1$first_time -
        peak2$first_time

      timing_shift_peak_midpoint <-
        peak1$midpoint_time -
        peak2$midpoint_time

      # ----------------------------------------------------------------------
      # TEMPORAL CENTER OF MASS
      # ----------------------------------------------------------------------

      com1 <- center_of_mass(x1)
      com2 <- center_of_mass(x2)

      center_of_mass_shift <- if (
        com1$defined &&
        com2$defined
      ) {
        com1$value -
          com2$value
      } else {
        NA_real_
      }

      # ----------------------------------------------------------------------
      # HIGH-EXPOSURE OVERLAP
      # ----------------------------------------------------------------------

      overlap <- high_overlap(
        x1,
        x2
      )

      # ----------------------------------------------------------------------
      # GROUP SIZES
      # ----------------------------------------------------------------------

      n_profiles1 <- if (grouped) {
        unname(
          group_sizes[[name1]]
        )
      } else {
        1L
      }

      n_profiles2 <- if (grouped) {
        unname(
          group_sizes[[name2]]
        )
      } else {
        1L
      }

      # ----------------------------------------------------------------------
      # OUTPUT ROW
      # ----------------------------------------------------------------------

      data.frame(
        exposure1 = name1,
        exposure2 = name2,
        n_profiles1 = as.integer(
          n_profiles1
        ),
        n_profiles2 = as.integer(
          n_profiles2
        ),
        n_time = n_time,

        L1 = L1,
        L2 = L2,
        Linf = Linf,
        MAE = MAE,
        RMSE = RMSE,

        sum1 = sum1,
        sum2 = sum2,
        diff_sum = sum1 - sum2,

        mean1 = mean1,
        mean2 = mean2,
        diff_mean = mean1 - mean2,

        sd1 = sd1,
        sd2 = sd2,
        diff_sd = sd1 - sd2,

        corr = corr,
        cor_method = cor_method,
        corr_defined =
          varying1 && varying2,

        peak_value1 = peak1$value,
        peak_value2 = peak2$value,

        peak_time1 = peak1$first_time,
        peak_time2 = peak2$first_time,
        timing_shift_peak =
          timing_shift_peak,

        peak_time_first1 =
          peak1$first_time,
        peak_time_last1 =
          peak1$last_time,
        peak_time_midpoint1 =
          peak1$midpoint_time,
        n_peak1 = peak1$n,
        peak_tied1 = peak1$tied,

        peak_time_first2 =
          peak2$first_time,
        peak_time_last2 =
          peak2$last_time,
        peak_time_midpoint2 =
          peak2$midpoint_time,
        n_peak2 = peak2$n,
        peak_tied2 = peak2$tied,

        timing_shift_peak_midpoint =
          timing_shift_peak_midpoint,

        peak_lag_first1 =
          peak1$first_lag,
        peak_lag_last1 =
          peak1$last_lag,
        peak_lag_midpoint1 =
          peak1$midpoint_lag,

        peak_lag_first2 =
          peak2$first_lag,
        peak_lag_last2 =
          peak2$last_lag,
        peak_lag_midpoint2 =
          peak2$midpoint_lag,

        com1 = com1$value,
        com2 = com2$value,
        com_defined1 =
          com1$defined,
        com_defined2 =
          com2$defined,
        com_reason1 =
          com1$reason,
        com_reason2 =
          com2$reason,
        center_of_mass_shift =
          center_of_mass_shift,

        threshold1 =
          overlap$threshold1,
        threshold2 =
          overlap$threshold2,
        n_high1 =
          overlap$n_high1,
        n_high2 =
          overlap$n_high2,
        n_high_intersection =
          overlap$n_intersection,
        n_high_union =
          overlap$n_union,
        overlap_above_q =
          overlap$overlap,
        q = q,

        stringsAsFactors = FALSE
      )
    }
  )

  out <- do.call(
    rbind,
    global_rows
  )
  rownames(out) <- NULL

  attach_metadata(out)
}
