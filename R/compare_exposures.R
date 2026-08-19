#' Compare exposure profiles through time
#'
#' Compares exposure profiles between two or more scenarios, epidemics,
#' or simulated exposure histories. The function supports both global
#' aggregated comparisons and time-wise comparisons, enabling detailed
#' temporal analysis of exposure differences.
#'
#' Exposure profiles must be supplied in chronological order, from the
#' earliest observation to the most recent observation. The first element
#' of each exposure vector therefore corresponds to the earliest time
#' position, and the last element corresponds to the most recent time
#' position.
#'
#' The function also supports grouping of exposure profiles, such as multiple
#' epidemics belonging to the same treatment or scenario. When grouping is
#' requested, profiles are aggregated at each time position using either
#' the median or mean before pairwise comparisons are performed.
#'
#' @param exposures Named list of numeric vectors. Each element represents
#'   an exposure profile supplied in chronological order, from the earliest
#'   to the most recent observation. All vectors must have the same length.
#'   Alternatively, a single numeric vector is accepted and converted to
#'   a named list.
#'
#' @param exposure1 Legacy numeric vector used for backward compatibility.
#'   Must be supplied together with `exposure2` when `exposures` is `NULL`.
#'
#' @param exposure2 Legacy numeric vector used for backward compatibility.
#'   Must be supplied together with `exposure1` when `exposures` is `NULL`.
#'
#' @param time Optional numeric vector defining chronological time positions.
#'   The vector must have the same length as each exposure profile and must
#'   be ordered from the earliest to the most recent observation. If `NULL`,
#'   the sequence `0:(n - 1)` is used.
#'
#' @param group Optional vector defining group membership for each exposure
#'   profile. Its length must equal the number of profiles in `exposures`.
#'   When provided, profiles are first aggregated within each group and
#'   comparisons are performed between the resulting group-level profiles.
#'
#' @param agg_fun Aggregation function used when `group` is provided.
#'   Available options are `"median"` and `"mean"`. The median is the default
#'   and is generally more robust to extreme exposure profiles.
#'
#' @param mode Type of comparison:
#'   \describe{
#'     \item{"global"}{
#'       Returns one row per pairwise comparison, including distance,
#'       correlation, magnitude, timing, center-of-mass, and overlap metrics.
#'     }
#'     \item{"timewise"}{
#'       Returns one row per chronological time position and pairwise
#'       comparison, retaining the full temporal resolution of the profiles.
#'     }
#'   }
#'
#' @param q Numeric scalar between 0 and 1. Quantile threshold used to define
#'   high-exposure time positions when calculating overlap. The default is
#'   `0.8`.
#'
#' @param eps Positive numeric scalar used to avoid division by zero in ratio
#'   and center-of-mass calculations. The default is `1e-12`.
#'
#' @return A data.frame.
#'
#' If `mode = "global"`, the output contains one row per pairwise comparison
#' with:
#' \itemize{
#'   \item `exposure1` and `exposure2`: names of the compared profiles;
#'   \item `n_time`: number of time positions;
#'   \item `L1`, `L2`, and `Linf`: distance metrics;
#'   \item `diff_sum`, `diff_mean`, and `diff_sd`: differences in profile
#'     summary statistics;
#'   \item `corr`: Pearson correlation between profiles;
#'   \item `peak_time1` and `peak_time2`: chronological positions of maximum
#'     exposure;
#'   \item `timing_shift_peak`: difference between peak times;
#'   \item `com1` and `com2`: temporal centers of mass;
#'   \item `center_of_mass_shift`: difference between centers of mass;
#'   \item `overlap_above_q`: overlap between high-exposure time positions;
#'   \item `q`: quantile used to define high exposure.
#' }
#'
#' If `mode = "timewise"`, the output contains one row per time position
#' and pairwise comparison with:
#' \itemize{
#'   \item `exposure1` and `exposure2`;
#'   \item `time`;
#'   \item `value1` and `value2`;
#'   \item `diff`;
#'   \item `abs_diff`;
#'   \item `ratio`.
#' }
#'
#' @details
#' This function operates exclusively on chronological exposure profiles.
#' No reversal is applied internally. Users must supply each profile from
#' the earliest observation to the most recent observation.
#'
#' Global comparisons summarize overall differences and similarities across
#' the complete exposure history. Time-wise comparisons preserve the full
#' chronological resolution and are suitable for identifying and plotting
#' when differences occurred.
#'
#' When observational datasets contain many epidemic-specific profiles,
#' `group` can be used to create representative group-level profiles before
#' comparison. Aggregation is performed independently at each time position.
#'
#' The ratio at each time position is calculated as:
#'
#' \deqn{
#'   \frac{x_1}{x_2 + \epsilon}
#' }
#'
#' where \eqn{\epsilon} is controlled by `eps`.
#'
#' The temporal center of mass is calculated after replacing negative
#' exposure values with zero. Consequently, this metric is primarily
#' meaningful for non-negative exposure variables.
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
    q = 0.8,
    eps = 1e-12
) {

  # ============================================================
  # MATCH ARGUMENTS
  # ============================================================

  mode <- match.arg(mode)
  agg_fun <- match.arg(agg_fun)

  # ============================================================
  # BASIC ARGUMENT VALIDATION
  # ============================================================

  if (
    !is.numeric(q) ||
    length(q) != 1L ||
    !is.finite(q) ||
    q < 0 ||
    q > 1
  ) {

    stop(
      "`q` must be a single finite numeric value between 0 and 1."
    )
  }

  if (
    !is.numeric(eps) ||
    length(eps) != 1L ||
    !is.finite(eps) ||
    eps <= 0
  ) {

    stop(
      "`eps` must be a single positive finite numeric value."
    )
  }

  # ============================================================
  # BACKWARD COMPATIBILITY
  # ============================================================

  if (is.null(exposures)) {

    if (
      is.null(exposure1) ||
      is.null(exposure2)
    ) {

      stop(
        paste0(
          "Provide either `exposures` or both `exposure1` and ",
          "`exposure2`."
        )
      )
    }

    exposures <- list(
      exposure1 = exposure1,
      exposure2 = exposure2
    )

  } else if (is.numeric(exposures)) {

    exposures <- list(
      exposure = exposures
    )
  }

  # ============================================================
  # VALIDATE EXPOSURE LIST
  # ============================================================

  if (!is.list(exposures)) {

    stop(
      "`exposures` must be a named list of numeric vectors."
    )
  }

  if (length(exposures) == 0L) {

    stop(
      "`exposures` cannot be empty."
    )
  }

  if (
    is.null(names(exposures)) ||
    anyNA(names(exposures)) ||
    any(names(exposures) == "")
  ) {

    stop(
      "All exposure profiles in `exposures` must have names."
    )
  }

  if (anyDuplicated(names(exposures))) {

    stop(
      "Exposure-profile names must be unique."
    )
  }

  lengths <- vapply(
    exposures,
    length,
    integer(1)
  )

  if (length(unique(lengths)) != 1L) {

    stop(
      "All exposure profiles must have the same length."
    )
  }

  n <- unname(lengths[1])

  if (n < 2L) {

    stop(
      "Exposure profiles must have length greater than or equal to 2."
    )
  }

  for (nm in names(exposures)) {

    x <- exposures[[nm]]

    if (!is.numeric(x)) {

      stop(
        paste0(
          "Exposure profile '",
          nm,
          "' must be numeric."
        )
      )
    }

    if (any(!is.finite(x))) {

      stop(
        paste0(
          "Exposure profile '",
          nm,
          "' must contain only finite values."
        )
      )
    }

    exposures[[nm]] <- as.numeric(x)
  }

  # ============================================================
  # CHRONOLOGICAL TIME
  # ============================================================

  if (is.null(time)) {

    time <- 0:(n - 1L)
  }

  if (!is.numeric(time)) {

    stop(
      "`time` must be a numeric vector."
    )
  }

  if (length(time) != n) {

    stop(
      "`time` must have the same length as each exposure profile."
    )
  }

  if (any(!is.finite(time))) {

    stop(
      "`time` must contain only finite values."
    )
  }

  if (anyDuplicated(time)) {

    stop(
      "`time` must contain unique chronological positions."
    )
  }

  time <- as.numeric(time)

  if (is.unsorted(time, strictly = TRUE)) {

    stop(
      paste0(
        "`time` must be supplied in strictly increasing chronological ",
        "order, from the earliest to the most recent observation."
      )
    )
  }

  # ============================================================
  # GROUP AGGREGATION
  # ============================================================

  if (!is.null(group)) {

    if (length(group) != length(exposures)) {

      stop(
        "`group` must have the same length as the number of exposure profiles."
      )
    }

    if (anyNA(group)) {

      stop(
        "`group` cannot contain missing values."
      )
    }

    group <- as.character(group)

    group_levels <- unique(group)

    if (length(group_levels) < 2L) {

      stop(
        "`group` must define at least two groups for comparison."
      )
    }

    mat <- do.call(
      rbind,
      exposures
    )

    agg_fun_use <- if (
      identical(agg_fun, "median")
    ) {
      stats::median
    } else {
      base::mean
    }

    exposures <- lapply(
      group_levels,
      function(g) {

        rows <- group == g

        apply(
          mat[rows, , drop = FALSE],
          MARGIN = 2,
          FUN = agg_fun_use
        )
      }
    )

    names(exposures) <- group_levels
  }

  # ============================================================
  # REQUIRE AT LEAST TWO PROFILES AFTER AGGREGATION
  # ============================================================

  if (length(exposures) < 2L) {

    stop(
      paste0(
        "At least two exposure profiles or groups are required ",
        "for comparison."
      )
    )
  }

  combinations <- utils::combn(
    names(exposures),
    2,
    simplify = FALSE
  )

  # ============================================================
  # TIME-WISE MODE
  # ============================================================

  if (identical(mode, "timewise")) {

    out_list <- lapply(
      combinations,
      function(comparison) {

        name1 <- comparison[1]
        name2 <- comparison[2]

        x1 <- exposures[[name1]]
        x2 <- exposures[[name2]]

        denominator <- x2 + eps

        ratio <- x1 / denominator

        data.frame(
          exposure1 = name1,
          exposure2 = name2,
          time = time,
          value1 = x1,
          value2 = x2,
          diff = x1 - x2,
          abs_diff = abs(x1 - x2),
          ratio = ratio,
          stringsAsFactors = FALSE
        )
      }
    )

    out <- do.call(
      rbind,
      out_list
    )

    rownames(out) <- NULL

    return(out)
  }

  # ============================================================
  # GLOBAL MODE
  # ============================================================

  out_list <- lapply(
    combinations,
    function(comparison) {

      name1 <- comparison[1]
      name2 <- comparison[2]

      x1 <- exposures[[name1]]
      x2 <- exposures[[name2]]

      difference <- x1 - x2
      absolute_difference <- abs(difference)

      # --------------------------------------------------------
      # DISTANCE METRICS
      # --------------------------------------------------------

      L1 <- sum(
        absolute_difference
      )

      L2 <- sqrt(
        sum(difference^2)
      )

      Linf <- max(
        absolute_difference
      )

      # --------------------------------------------------------
      # PROFILE SUMMARIES
      # --------------------------------------------------------

      sum1 <- sum(x1)
      sum2 <- sum(x2)

      mean1 <- mean(x1)
      mean2 <- mean(x2)

      sd1 <- stats::sd(x1)
      sd2 <- stats::sd(x2)

      # --------------------------------------------------------
      # CORRELATION
      # --------------------------------------------------------

      corr <- if (
        stats::sd(x1) == 0 ||
        stats::sd(x2) == 0
      ) {

        NA_real_

      } else {

        suppressWarnings(
          stats::cor(
            x1,
            x2
          )
        )
      }

      # --------------------------------------------------------
      # PEAK TIMING
      # --------------------------------------------------------

      peak_time1 <- time[
        which.max(x1)
      ]

      peak_time2 <- time[
        which.max(x2)
      ]

      timing_shift_peak <- (
        peak_time1 -
          peak_time2
      )

      # --------------------------------------------------------
      # TEMPORAL CENTER OF MASS
      # --------------------------------------------------------

      weight1 <- pmax(
        x1,
        0
      )

      weight2 <- pmax(
        x2,
        0
      )

      center_of_mass1 <- if (
        sum(weight1) <= eps
      ) {

        NA_real_

      } else {

        sum(
          time * weight1
        ) / sum(weight1)
      }

      center_of_mass2 <- if (
        sum(weight2) <= eps
      ) {

        NA_real_

      } else {

        sum(
          time * weight2
        ) / sum(weight2)
      }

      center_of_mass_shift <- (
        center_of_mass1 -
          center_of_mass2
      )

      # --------------------------------------------------------
      # HIGH-EXPOSURE OVERLAP
      # --------------------------------------------------------

      threshold1 <- stats::quantile(
        x1,
        probs = q,
        names = FALSE,
        na.rm = TRUE
      )

      threshold2 <- stats::quantile(
        x2,
        probs = q,
        names = FALSE,
        na.rm = TRUE
      )

      high1 <- x1 >= threshold1
      high2 <- x2 >= threshold2

      union_high <- sum(
        high1 | high2
      )

      overlap_above_q <- if (
        union_high == 0
      ) {

        NA_real_

      } else {

        sum(
          high1 & high2
        ) / union_high
      }

      # --------------------------------------------------------
      # OUTPUT
      # --------------------------------------------------------

      data.frame(
        exposure1 = name1,
        exposure2 = name2,
        n_time = n,
        L1 = L1,
        L2 = L2,
        Linf = Linf,
        diff_sum = sum1 - sum2,
        diff_mean = mean1 - mean2,
        diff_sd = sd1 - sd2,
        corr = corr,
        peak_time1 = peak_time1,
        peak_time2 = peak_time2,
        timing_shift_peak = timing_shift_peak,
        com1 = center_of_mass1,
        com2 = center_of_mass2,
        center_of_mass_shift = center_of_mass_shift,
        overlap_above_q = overlap_above_q,
        q = q,
        stringsAsFactors = FALSE
      )
    }
  )

  out <- do.call(
    rbind,
    out_list
  )

  rownames(out) <- NULL

  return(out)
}
