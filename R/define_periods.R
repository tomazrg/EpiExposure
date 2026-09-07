#' Define epidemiological lag periods
#'
#' Creates contiguous, non-overlapping lag periods for period-specific DLNM
#' summaries and scenario definitions.
#'
#' Periods are defined on the EpiExposure retrospective lag scale:
#'
#' - lag 0 = the most recent exposure observation;
#' - increasing lag = progressively older exposure observations;
#' - `max_lag` = the oldest exposure observation included in the fitted history.
#'
#' @param max_lag Non-negative integer maximum lag. The standard form is a
#'   scalar such as `85`. For compatibility with lag-range metadata,
#'   `c(0, 85)` is also accepted. EpiExposure v1 requires lag histories to start
#'   at 0.
#' @param cuts Optional numeric vector of **inclusive period-end lags**.
#'
#'   For example, with `max_lag = 20`:
#'
#'   ```
#'   cuts = c(7, 14)
#'   ```
#'
#'   creates:
#'
#'   ```
#'   0--7
#'   8--14
#'   15--20
#'   ```
#'
#'   Every cut must therefore be a unique finite integer strictly between
#'   `0` and `max_lag`. Cut order in the input is not meaningful; valid cuts
#'   are arranged in increasing lag order before periods are constructed.
#'
#'   `NULL` or `numeric(0)` requests one cumulative period from lag 0 through
#'   `max_lag`.
#' @param prefix Non-empty character scalar used to label periods. With the
#'   default `prefix = "W"`, periods are labeled `"W1"`, `"W2"`, and so on.
#'
#' @return A data frame with exactly three columns:
#'
#'   - `period`: unique period label;
#'   - `lag_start`: first lag included in the period;
#'   - `lag_end`: last lag included in the period.
#'
#'   Period bounds are integer and inclusive. Rows are ordered from the most
#'   recent lag interval to the oldest.
#'
#'   The returned object also stores:
#'
#'   - `attr(x, "epiexposure_max_lag")`: validated maximum lag;
#'   - `attr(x, "epiexposure_cuts")`: validated internal cut points;
#'   - `attr(x, "epiexposure_n_lags")`: total number of represented lags,
#'     `max_lag + 1`;
#'   - `attr(x, "epiexposure_period_contract")`:
#'     `"contiguous_inclusive_retrospective_lags"`;
#'   - `attr(x, "epiexposure_lag_scale")`:
#'     `"lag_0_most_recent"`;
#'   - `attr(x, "epiexposure_full_coverage")`: `TRUE`.
#'
#' @details
#' ## Inclusive cut semantics
#'
#' `cuts` are not arbitrary break coordinates and are not the first lag of the
#' following period. Each cut is the **inclusive final lag of the current
#' period**.
#'
#' Thus:
#'
#' ```
#' define_periods(
#'   max_lag = 20,
#'   cuts = c(7, 14)
#' )
#' ```
#'
#' returns:
#'
#' ```
#'   period lag_start lag_end
#'   W1             0       7
#'   W2             8      14
#'   W3            15      20
#' ```
#'
#' Every lag from `0` through `max_lag` appears exactly once. There are no gaps
#' and no overlaps.
#'
#' ## Single cumulative period
#'
#' When `cuts = NULL` or `cuts = numeric(0)`, the full lag history is represented
#' by one period:
#'
#' ```
#'   period lag_start lag_end
#'   W1             0 max_lag
#' ```
#'
#' This also supports `max_lag = 0`, in which case the single valid period is
#' `0--0`.
#'
#' ## Validation philosophy
#'
#' Earlier versions silently converted cut points with `as.integer()` and
#' discarded cuts outside the fitted lag range. The current EpiExposure
#' contract does not alter or ignore invalid period definitions silently.
#'
#' Non-integer, duplicated, non-finite, boundary (`0` or `max_lag`), or
#' out-of-range cut points produce explicit errors. This prevents downstream
#' functions such as `summarise_effects()` and scenario simulation from using a
#' period definition different from the one requested by the analyst.
#'
#' ## Relationship with `summarise_effects()`
#'
#' For period-specific effects, `summarise_effects()` treats both
#' `lag_start` and `lag_end` as inclusive. Therefore a period `8--14`
#' represents exactly the lags:
#'
#' ```
#' 8, 9, 10, 11, 12, 13, 14
#' ```
#'
#' `define_periods()` constructs exactly that contract and orders periods by
#' increasing retrospective lag.
#'
#' @export
define_periods <- function(
    max_lag,
    cuts = NULL,
    prefix = "W"
) {

  # ==========================================================================
  # SMALL VALIDATORS
  # ==========================================================================

  is_integerish <- function(x) {
    is.numeric(x) &&
      length(x) >= 1L &&
      !anyNA(x) &&
      all(is.finite(x)) &&
      all(x == as.integer(x))
  }

  valid_prefix <- function(x) {
    is.character(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      nzchar(trimws(x))
  }

  # ==========================================================================
  # VALIDATE max_lag
  # ==========================================================================

  if (!is.numeric(max_lag) ||
      !length(max_lag) ||
      length(max_lag) > 2L ||
      anyNA(max_lag) ||
      any(!is.finite(max_lag)) ||
      any(max_lag < 0) ||
      any(max_lag != as.integer(max_lag))) {
    stop(
      "`max_lag` must be one non-negative integer or the lag range `c(0, L)`.",
      call. = FALSE
    )
  }

  max_lag <- as.integer(max_lag)

  if (length(max_lag) == 2L) {
    if (max_lag[1L] != 0L) {
      stop(
        "When `max_lag` is supplied as a two-value lag range, it must start ",
        "at 0, for example `c(0, 85)`.",
        call. = FALSE
      )
    }

    if (max_lag[2L] < max_lag[1L]) {
      stop(
        "The two-value `max_lag` range must be ordered as `c(0, L)`.",
        call. = FALSE
      )
    }

    maximum_lag <- max_lag[2L]
  } else {
    maximum_lag <- max_lag[1L]
  }

  # ==========================================================================
  # VALIDATE prefix
  # ==========================================================================

  if (!valid_prefix(prefix)) {
    stop(
      "`prefix` must be one non-empty character string.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # VALIDATE cuts
  # ==========================================================================

  if (is.null(cuts)) {
    cuts_valid <- integer(0)

  } else {
    if (!is.numeric(cuts)) {
      stop(
        "`cuts` must be NULL or a numeric vector of inclusive period-end lags.",
        call. = FALSE
      )
    }

    if (!length(cuts)) {
      cuts_valid <- integer(0)

    } else {
      if (anyNA(cuts) ||
          any(!is.finite(cuts))) {
        stop(
          "`cuts` must contain only finite, non-missing values.",
          call. = FALSE
        )
      }

      if (any(cuts != as.integer(cuts))) {
        invalid <- unique(
          cuts[cuts != as.integer(cuts)]
        )

        stop(
          "`cuts` must contain integer lag values. Non-integer value(s): ",
          paste(
            utils::head(invalid, 10L),
            collapse = ", "
          ),
          if (length(invalid) > 10L) ", ..." else "",
          ".",
          call. = FALSE
        )
      }

      cuts_integer <- as.integer(cuts)

      if (anyDuplicated(cuts_integer)) {
        duplicated_values <- unique(
          cuts_integer[
            duplicated(cuts_integer) |
              duplicated(cuts_integer, fromLast = TRUE)
          ]
        )

        stop(
          "`cuts` must contain unique lag values. Duplicated value(s): ",
          paste(
            duplicated_values,
            collapse = ", "
          ),
          ".",
          call. = FALSE
        )
      }

      if (maximum_lag == 0L) {
        stop(
          "`cuts` cannot be supplied when `max_lag = 0`; the only valid lag ",
          "period is 0--0.",
          call. = FALSE
        )
      }

      invalid_boundary <- cuts_integer <= 0L |
        cuts_integer >= maximum_lag

      if (any(invalid_boundary)) {
        invalid <- cuts_integer[
          invalid_boundary
        ]

        stop(
          "Every cut must be strictly between 0 and `max_lag` because cuts are ",
          "inclusive period-end lags. Invalid value(s): ",
          paste(
            invalid,
            collapse = ", "
          ),
          ". For max_lag = ",
          maximum_lag,
          ", valid cuts are integers from 1 through ",
          maximum_lag - 1L,
          ".",
          call. = FALSE
        )
      }

      cuts_valid <- sort(cuts_integer)
    }
  }

  # ==========================================================================
  # CONSTRUCT PERIODS
  # ==========================================================================

  if (!length(cuts_valid)) {
    starts <- 0L
    ends <- maximum_lag

  } else {
    starts <- c(
      0L,
      cuts_valid + 1L
    )

    ends <- c(
      cuts_valid,
      maximum_lag
    )
  }

  if (length(starts) != length(ends)) {
    stop(
      "Internal error: mismatched lag-period boundaries.",
      call. = FALSE
    )
  }

  if (any(starts < 0L) ||
      any(ends > maximum_lag) ||
      any(starts > ends)) {
    stop(
      "Internal error: generated lag-period boundaries are invalid.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # STRICT COVERAGE / CONTIGUITY CHECK
  # ==========================================================================

  if (starts[1L] != 0L ||
      ends[length(ends)] != maximum_lag) {
    stop(
      "Internal error: generated periods do not cover the complete fitted lag ",
      "range.",
      call. = FALSE
    )
  }

  if (length(starts) > 1L) {
    expected_starts <- ends[
      seq_len(length(ends) - 1L)
    ] + 1L

    actual_starts <- starts[
      2:length(starts)
    ]

    if (!identical(
      as.integer(actual_starts),
      as.integer(expected_starts)
    )) {
      stop(
        "Internal error: generated periods contain a gap or overlap.",
        call. = FALSE
      )
    }
  }

  represented_lags <- unlist(
    Map(
      function(start, end) {
        seq.int(start, end)
      },
      starts,
      ends
    ),
    use.names = FALSE
  )

  expected_lags <- seq.int(
    0L,
    maximum_lag
  )

  if (!identical(
    as.integer(represented_lags),
    as.integer(expected_lags)
  )) {
    stop(
      "Internal error: each lag from 0 through `max_lag` must occur exactly ",
      "once across the generated periods.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # OUTPUT
  # ==========================================================================

  out <- data.frame(
    period = paste0(
      prefix,
      seq_along(starts)
    ),
    lag_start = as.integer(starts),
    lag_end = as.integer(ends),
    stringsAsFactors = FALSE
  )

  if (anyNA(out$period) ||
      any(!nzchar(out$period)) ||
      anyDuplicated(out$period)) {
    stop(
      "Internal error: generated period labels are invalid.",
      call. = FALSE
    )
  }

  attr(
    out,
    "epiexposure_max_lag"
  ) <- maximum_lag

  attr(
    out,
    "epiexposure_cuts"
  ) <- cuts_valid

  attr(
    out,
    "epiexposure_n_lags"
  ) <- maximum_lag + 1L

  attr(
    out,
    "epiexposure_period_contract"
  ) <- "contiguous_inclusive_retrospective_lags"

  attr(
    out,
    "epiexposure_lag_scale"
  ) <- "lag_0_most_recent"

  attr(
    out,
    "epiexposure_full_coverage"
  ) <- TRUE

  attr(
    out,
    "epiexposure_period_order"
  ) <- "increasing_lag"

  out
}
