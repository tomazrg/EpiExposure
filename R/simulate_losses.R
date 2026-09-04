#' Simulate yield and economic losses
#'
#' Simulates predicted yield, relative yield, yield loss, and economic loss
#' from a user-specified response variable. Optional lower and upper response
#' bounds are propagated through all derived loss metrics.
#'
#' Each row of the original data is retained and expanded across combinations
#' of attainable yield, price, and optional parameter simulations.
#'
#' @param data A data.frame containing the response variable and any additional
#'   scenario variables.
#' @param y Character string giving the name of the numeric response column
#'   used to calculate yield and economic losses.
#' @param lower Optional character string giving the column containing the
#'   lower bound of `y`. If supplied, `upper` must also be supplied.
#' @param upper Optional character string giving the column containing the
#'   upper bound of `y`. If supplied, `lower` must also be supplied.
#' @param slope Numeric. Regression slope relating `y` to yield loss.
#'   A single value uses a fixed slope. Two values are interpreted as the
#'   minimum and maximum of a uniform distribution. Alternatively, a vector
#'   of length `n` can contain user-supplied parameter draws.
#' @param intercept Numeric. Yield-model intercept. A single value uses a
#'   fixed intercept. Two values are interpreted as the minimum and maximum
#'   of a uniform distribution. Alternatively, a vector of length `n`
#'   can contain user-supplied parameter draws.
#' @param attainable_yield Numeric vector of attainable yield values.
#' @param price Numeric vector of commodity prices.
#' @param n Number of parameter simulations. Default is 1.
#' @param random_sd Standard deviation of the normally distributed random
#'   effect used in the yield model. Default is 0 (no random effect).
#' @param y_multiplier Numeric multiplier applied internally to `y`, `lower`,
#'   and `upper` before loss calculations. Default is 1. For example, use
#'   `y_multiplier = 100` when model predictions are proportions (0-1) but
#'   the yield-loss coefficients were estimated using percentages (0-100).
#' @param seed Optional random seed.
#'
#' @return A data.frame containing all original columns from `data` plus
#'   simulation identifiers, parameter draws, transformed response values,
#'   and loss metrics. When `lower` and `upper` are supplied, lower and upper
#'   bounds are returned for every derived metric.
#'
#' @details
#' The response used internally is:
#'
#'   response = y * y_multiplier
#'
#' Predicted yield is:
#'
#'   intercept - slope * response - random_effect
#'
#' and is constrained to be non-negative.
#'
#' Relative yield proportion is:
#'
#'   predicted_yield / intercept
#'
#' Relative yield percentage is:
#'
#'   relative_yield_proportion * 100
#'
#' Yield-loss proportion is:
#'
#'   (slope / intercept) * response
#'
#' and is constrained between 0 and 1.
#'
#' Yield-loss percentage is:
#'
#'   yield_loss_proportion * 100
#'
#' Absolute yield loss is:
#'
#'   attainable_yield * yield_loss_proportion
#'
#' Economic loss is:
#'
#'   (yield_loss / 1000) * price
#'
#' The division by 1000 is fixed internally, assuming yield is expressed in
#' kg per unit area and price is expressed per metric ton.
#'
#' When response bounds are supplied, the function evaluates each metric at
#' both response bounds and then uses the minimum and maximum resulting values
#' as the metric-specific lower and upper bounds. This is important because
#' increasing disease response increases loss metrics but decreases predicted
#' and relative yield.
#'
#' @export
simulate_losses <- function(
    data,
    y,
    lower = NULL,
    upper = NULL,
    slope = 100,
    intercept = 11142.94,
    attainable_yield = seq(4000, 12000, by = 50),
    price = seq(100, 300, by = 10),
    n = 1,
    random_sd = 0,
    y_multiplier = 1,
    seed = NULL
) {

  # ============================================================
  # 1. Validate input data
  # ============================================================

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  if (!is.character(y) || length(y) != 1L || is.na(y) || !nzchar(y)) {
    stop("`y` must be a single character string giving the response column.")
  }

  if (!y %in% names(data)) {
    stop("Column `", y, "` was not found in `data`.")
  }

  if (!is.numeric(data[[y]])) {
    stop("The response column `", y, "` must be numeric.")
  }

  # lower and upper are optional, but must be supplied together
  has_bounds <- !is.null(lower) || !is.null(upper)

  if (has_bounds) {

    if (is.null(lower) || is.null(upper)) {
      stop("`lower` and `upper` must be supplied together.")
    }

    if (!is.character(lower) ||
        length(lower) != 1L ||
        is.na(lower) ||
        !nzchar(lower)) {
      stop("`lower` must be NULL or a single character string.")
    }

    if (!is.character(upper) ||
        length(upper) != 1L ||
        is.na(upper) ||
        !nzchar(upper)) {
      stop("`upper` must be NULL or a single character string.")
    }

    if (!lower %in% names(data)) {
      stop("Column `", lower, "` was not found in `data`.")
    }

    if (!upper %in% names(data)) {
      stop("Column `", upper, "` was not found in `data`.")
    }

    if (!is.numeric(data[[lower]])) {
      stop("The lower-bound column `", lower, "` must be numeric.")
    }

    if (!is.numeric(data[[upper]])) {
      stop("The upper-bound column `", upper, "` must be numeric.")
    }

    finite_bounds <- is.finite(data[[lower]]) & is.finite(data[[upper]])

    if (any(
      data[[lower]][finite_bounds] >
      data[[upper]][finite_bounds]
    )) {
      stop(
        "Values in `lower` must be less than or equal to values in `upper`."
      )
    }

    finite_all <- is.finite(data[[lower]]) &
      is.finite(data[[y]]) &
      is.finite(data[[upper]])

    if (any(
      data[[y]][finite_all] < data[[lower]][finite_all] |
      data[[y]][finite_all] > data[[upper]][finite_all]
    )) {
      warning(
        "Some values of `y` fall outside the corresponding `lower` and ",
        "`upper` bounds."
      )
    }
  }

  if (!is.numeric(slope) ||
      length(slope) == 0L ||
      any(!is.finite(slope))) {
    stop("`slope` must contain finite numeric values.")
  }

  if (!is.numeric(intercept) ||
      length(intercept) == 0L ||
      any(!is.finite(intercept))) {
    stop("`intercept` must contain finite numeric values.")
  }

  if (!is.numeric(attainable_yield) ||
      length(attainable_yield) == 0L ||
      any(!is.finite(attainable_yield))) {
    stop("`attainable_yield` must contain finite numeric values.")
  }

  if (any(attainable_yield < 0)) {
    stop("`attainable_yield` cannot contain negative values.")
  }

  if (!is.numeric(price) ||
      length(price) == 0L ||
      any(!is.finite(price))) {
    stop("`price` must contain finite numeric values.")
  }

  if (any(price < 0)) {
    stop("`price` cannot contain negative values.")
  }

  if (length(n) != 1L ||
      !is.numeric(n) ||
      !is.finite(n) ||
      n < 1 ||
      n != as.integer(n)) {
    stop("`n` must be a positive integer.")
  }

  n <- as.integer(n)

  if (length(random_sd) != 1L ||
      !is.numeric(random_sd) ||
      !is.finite(random_sd) ||
      random_sd < 0) {
    stop("`random_sd` must be a single non-negative numeric value.")
  }

  if (length(y_multiplier) != 1L ||
      !is.numeric(y_multiplier) ||
      !is.finite(y_multiplier) ||
      y_multiplier <= 0) {
    stop("`y_multiplier` must be a single positive numeric value.")
  }


  # ============================================================
  # 2. Check output-column conflicts
  # ============================================================

  base_output_columns <- c(
    ".loss_row_id",
    ".sim",
    "slope_sim",
    "intercept_sim",
    "random_effect",
    "response_used",
    "predicted_yield",
    "relative_yield_proportion",
    "relative_yield_pct",
    "yield_loss_proportion",
    "yield_loss_pct",
    "attainable_yield",
    "yield_loss",
    "price",
    "economic_loss"
  )

  bound_output_columns <- c(
    "response_lower_used",
    "response_upper_used",
    "predicted_yield_lower",
    "predicted_yield_upper",
    "relative_yield_proportion_lower",
    "relative_yield_proportion_upper",
    "relative_yield_pct_lower",
    "relative_yield_pct_upper",
    "yield_loss_proportion_lower",
    "yield_loss_proportion_upper",
    "yield_loss_pct_lower",
    "yield_loss_pct_upper",
    "yield_loss_lower",
    "yield_loss_upper",
    "economic_loss_lower",
    "economic_loss_upper"
  )

  new_columns <- if (has_bounds) {
    c(base_output_columns, bound_output_columns)
  } else {
    base_output_columns
  }

  conflicts <- intersect(
    names(data),
    new_columns
  )

  if (length(conflicts) > 0L) {
    stop(
      "The following columns already exist in `data`: ",
      paste(conflicts, collapse = ", "),
      ". Rename them before using `simulate_losses()`."
    )
  }


  # ============================================================
  # 3. Helper to generate parameter values
  # ============================================================

  draw_parameter <- function(x, n, name) {

    # Fixed coefficient
    if (length(x) == 1L) {
      return(rep(as.numeric(x), n))
    }

    # Two values = uniform range
    if (length(x) == 2L) {
      return(
        stats::runif(
          n,
          min = min(x),
          max = max(x)
        )
      )
    }

    # User-supplied parameter draws
    if (length(x) == n) {
      return(as.numeric(x))
    }

    stop(
      "`", name, "` must have length 1, 2, or `n` (",
      n, ")."
    )
  }


  # ============================================================
  # 4. Seed without permanently changing the user's RNG state
  # ============================================================

  if (!is.null(seed)) {

    if (length(seed) != 1L ||
        !is.numeric(seed) ||
        !is.finite(seed)) {
      stop("`seed` must be NULL or a single finite numeric value.")
    }

    seed_existed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )

    if (seed_existed) {
      old_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }

    on.exit(
      {
        if (seed_existed) {
          assign(
            ".Random.seed",
            old_seed,
            envir = .GlobalEnv
          )
        } else if (
          exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
          )
        ) {
          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )

    set.seed(seed)
  }


  # ============================================================
  # 5. Simulate model parameters
  # ============================================================

  slope_sim <- draw_parameter(
    slope,
    n,
    "slope"
  )

  intercept_sim <- draw_parameter(
    intercept,
    n,
    "intercept"
  )

  if (any(intercept_sim <= 0)) {
    stop("`intercept` values must be greater than zero.")
  }

  if (random_sd > 0) {

    random_effect <- stats::rnorm(
      n,
      mean = 0,
      sd = random_sd
    )

  } else {

    random_effect <- rep(
      0,
      n
    )
  }


  # ============================================================
  # 6. Construct simulation grid
  # ============================================================

  simulation_grid <- expand.grid(
    attainable_yield = attainable_yield,
    price = price,
    .sim = seq_len(n),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  simulation_grid$slope_sim <-
    slope_sim[simulation_grid$.sim]

  simulation_grid$intercept_sim <-
    intercept_sim[simulation_grid$.sim]

  simulation_grid$random_effect <-
    random_effect[simulation_grid$.sim]


  # ============================================================
  # 7. Expand original data
  #
  # Every original row is retained and repeated for each
  # combination of attainable yield, price, and simulation.
  # ============================================================

  n_grid <- nrow(simulation_grid)

  original_index <- rep(
    seq_len(nrow(data)),
    each = n_grid
  )

  simulation_index <- rep(
    seq_len(n_grid),
    times = nrow(data)
  )

  out <- data[
    original_index,
    ,
    drop = FALSE
  ]

  out$.loss_row_id <- original_index

  out$.sim <-
    simulation_grid$.sim[simulation_index]

  out$slope_sim <-
    simulation_grid$slope_sim[simulation_index]

  out$intercept_sim <-
    simulation_grid$intercept_sim[simulation_index]

  out$random_effect <-
    simulation_grid$random_effect[simulation_index]

  out$attainable_yield <-
    simulation_grid$attainable_yield[simulation_index]

  out$price <-
    simulation_grid$price[simulation_index]


  # ============================================================
  # 8. Response values used in the loss equations
  # ============================================================

  out$response_used <-
    out[[y]] * y_multiplier

  if (has_bounds) {

    out$response_lower_used <-
      out[[lower]] * y_multiplier

    out$response_upper_used <-
      out[[upper]] * y_multiplier
  }


  # ============================================================
  # 9. Helper to calculate all loss metrics for one response
  # ============================================================

  calculate_metrics <- function(response) {

    # Predicted yield
    predicted_yield <-
      out$intercept_sim -
      out$slope_sim * response -
      out$random_effect

    predicted_yield <-
      pmax(
        predicted_yield,
        0
      )

    # Relative yield proportion
    relative_yield_proportion <-
      predicted_yield /
      out$intercept_sim

    relative_yield_proportion <-
      pmin(
        pmax(
          relative_yield_proportion,
          0
        ),
        1
      )

    # Relative yield percentage
    relative_yield_pct <-
      relative_yield_proportion *
      100

    # Yield-loss proportion
    yield_loss_proportion <-
      (
        out$slope_sim /
          out$intercept_sim
      ) *
      response

    yield_loss_proportion <-
      pmin(
        pmax(
          yield_loss_proportion,
          0
        ),
        1
      )

    # Yield-loss percentage
    yield_loss_pct <-
      yield_loss_proportion *
      100

    # Absolute yield loss
    yield_loss <-
      out$attainable_yield *
      yield_loss_proportion

    # Economic loss
    economic_loss <-
      (
        yield_loss /
          1000
      ) *
      out$price

    list(
      predicted_yield = predicted_yield,
      relative_yield_proportion = relative_yield_proportion,
      relative_yield_pct = relative_yield_pct,
      yield_loss_proportion = yield_loss_proportion,
      yield_loss_pct = yield_loss_pct,
      yield_loss = yield_loss,
      economic_loss = economic_loss
    )
  }


  # ============================================================
  # 10. Central estimates
  # ============================================================

  central_metrics <-
    calculate_metrics(
      out$response_used
    )

  out$predicted_yield <-
    central_metrics$predicted_yield

  out$relative_yield_proportion <-
    central_metrics$relative_yield_proportion

  out$relative_yield_pct <-
    central_metrics$relative_yield_pct

  out$yield_loss_proportion <-
    central_metrics$yield_loss_proportion

  out$yield_loss_pct <-
    central_metrics$yield_loss_pct

  out$yield_loss <-
    central_metrics$yield_loss

  out$economic_loss <-
    central_metrics$economic_loss


  # ============================================================
  # 11. Propagate lower and upper response bounds
  # ============================================================

  if (has_bounds) {

    metrics_at_lower <-
      calculate_metrics(
        out$response_lower_used
      )

    metrics_at_upper <-
      calculate_metrics(
        out$response_upper_used
      )

    # ----------------------------------------------------------
    # IMPORTANT:
    #
    # We use pmin()/pmax() for each derived metric instead of
    # assuming the direction of the relationship.
    #
    # For example:
    # - higher disease response -> higher yield loss
    # - higher disease response -> lower predicted yield
    #
    # This guarantees that every *_lower column is numerically
    # <= its corresponding *_upper column.
    # ----------------------------------------------------------

    out$predicted_yield_lower <-
      pmin(
        metrics_at_lower$predicted_yield,
        metrics_at_upper$predicted_yield
      )

    out$predicted_yield_upper <-
      pmax(
        metrics_at_lower$predicted_yield,
        metrics_at_upper$predicted_yield
      )

    out$relative_yield_proportion_lower <-
      pmin(
        metrics_at_lower$relative_yield_proportion,
        metrics_at_upper$relative_yield_proportion
      )

    out$relative_yield_proportion_upper <-
      pmax(
        metrics_at_lower$relative_yield_proportion,
        metrics_at_upper$relative_yield_proportion
      )

    out$relative_yield_pct_lower <-
      pmin(
        metrics_at_lower$relative_yield_pct,
        metrics_at_upper$relative_yield_pct
      )

    out$relative_yield_pct_upper <-
      pmax(
        metrics_at_lower$relative_yield_pct,
        metrics_at_upper$relative_yield_pct
      )

    out$yield_loss_proportion_lower <-
      pmin(
        metrics_at_lower$yield_loss_proportion,
        metrics_at_upper$yield_loss_proportion
      )

    out$yield_loss_proportion_upper <-
      pmax(
        metrics_at_lower$yield_loss_proportion,
        metrics_at_upper$yield_loss_proportion
      )

    out$yield_loss_pct_lower <-
      pmin(
        metrics_at_lower$yield_loss_pct,
        metrics_at_upper$yield_loss_pct
      )

    out$yield_loss_pct_upper <-
      pmax(
        metrics_at_lower$yield_loss_pct,
        metrics_at_upper$yield_loss_pct
      )

    out$yield_loss_lower <-
      pmin(
        metrics_at_lower$yield_loss,
        metrics_at_upper$yield_loss
      )

    out$yield_loss_upper <-
      pmax(
        metrics_at_lower$yield_loss,
        metrics_at_upper$yield_loss
      )

    out$economic_loss_lower <-
      pmin(
        metrics_at_lower$economic_loss,
        metrics_at_upper$economic_loss
      )

    out$economic_loss_upper <-
      pmax(
        metrics_at_lower$economic_loss,
        metrics_at_upper$economic_loss
      )
  }


  # ============================================================
  # 12. Restore row names
  # ============================================================

  rownames(out) <- NULL


  # ============================================================
  # 13. Attributes documenting simulation
  # ============================================================

  attr(out, "loss_spec") <- list(
    y = y,
    lower = lower,
    upper = upper,
    y_multiplier = y_multiplier,
    n = n,
    slope = slope,
    intercept = intercept,
    random_sd = random_sd,
    attainable_yield = attainable_yield,
    price = price,
    price_divisor = 1000,
    cap_loss = TRUE,
    truncate_yield = TRUE,
    keep_row_id = TRUE
  )


  out
}
