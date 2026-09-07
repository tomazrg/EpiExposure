#' Simulate yield and economic losses
#'
#' Converts a non-negative response variable into yield and economic losses
#' using an externally supplied linear yield-response relationship.
#'
#' Each row of `data` is retained and expanded across all combinations of
#' `attainable_yield`, `price`, and `n` loss-model parameter simulations.
#'
#' `simulate_losses()` does not fit an epidemiological model and does not draw
#' coefficients from an EpiExposure fit. The response supplied in `y` is treated
#' as already available. The `n` simulations refer only to uncertainty or
#' variability explicitly supplied for the external yield-loss relationship.
#'
#' @param data Non-empty data frame containing the response variable and any
#'   additional scenario/prediction identifiers that should be retained.
#' @param y Character scalar naming the numeric response column used in the
#'   yield-loss relationship.
#' @param lower Optional character scalar naming a lower response bound. If
#'   supplied, `upper` must also be supplied.
#' @param upper Optional character scalar naming an upper response bound. If
#'   supplied, `lower` must also be supplied.
#'
#'   Bounds are propagated deterministically through the loss equations for
#'   each loss-model simulation. They are therefore endpoint-propagated bounds,
#'   not automatically a joint confidence/credible interval after combining
#'   response uncertainty with yield-model parameter uncertainty.
#' @param slope Numeric specification for the non-negative yield-loss slope.
#'
#'   - length 1: fixed slope used in every loss simulation;
#'   - length 2: lower and upper limits of an independent Uniform distribution;
#'   - length `n`: user-supplied slope draws, paired by simulation index with
#'     `intercept`.
#'
#'   The response-yield relationship assumes increasing response cannot improve
#'   yield, so all supplied slope values must be non-negative.
#' @param intercept Numeric specification for the positive reference-yield
#'   intercept.
#'
#'   - length 1: fixed intercept;
#'   - length 2: lower and upper limits of an independent Uniform distribution;
#'   - length `n`: user-supplied intercept draws, paired by simulation index
#'     with `slope`.
#'
#'   All supplied intercept values must be strictly positive.
#' @param attainable_yield Unique non-negative numeric vector of attainable
#'   yields in kg/ha, used to scale the disease-attributable proportional
#'   loss into absolute yield loss in kg/ha.
#'
#'   `attainable_yield` is conceptually distinct from the intercept of the
#'   external calibration equation. The intercept defines the reference yield
#'   used to estimate proportional loss, whereas `attainable_yield` defines the
#'   scenario-specific yield to which that proportional loss is applied.
#' @param price Unique non-negative numeric vector of commodity prices in
#'   USD per metric ton (USD/t).
#' @param n Positive integer number of loss-model parameter simulations.
#'   Default is 1.
#'
#'   `n > 1` is useful only when at least one component of the external
#'   yield-loss model varies across simulations: `slope`, `intercept`, or the
#'   Gaussian intercept deviation controlled by `random_sd`.
#'
#'   If `n > 1` is requested but all three components are effectively fixed,
#'   the function returns the requested repeated simulation rows but emits a
#'   warning that they are identical and do not represent uncertainty
#'   propagation.
#'
#'   `n` does not represent EpiExposure model-coefficient draws unless the user
#'   explicitly supplies such external parameter draws through `slope` and
#'   `intercept`.
#' @param random_sd Non-negative numeric scalar controlling an optional
#'   simulation-level Gaussian deviation in the external yield-model intercept.
#'   Default is 0.
#'
#'   For simulation `s`,
#'
#'   \deqn{
#'     Y_{0,s} = intercept_s + b_s,\qquad
#'     b_s \sim N(0, random\_sd^2).
#'   }
#'
#'   The same `b_s` is shared by every row, attainable-yield value, and price
#'   within loss simulation `s`. This preserves matched scenario comparisons
#'   under the same external yield-model realization.
#'
#'   This quantity is **not** a random effect from the EpiExposure
#'   epidemiological fit. If residual/future-yield noise is not scientifically
#'   intended, leave `random_sd = 0`.
#' @param y_multiplier Positive numeric scalar applied to `y`, `lower`, and
#'   `upper` before the yield-loss equations. Default is 1.
#'
#'   For example, use `y_multiplier = 100` when the response is stored as a
#'   proportion from 0 to 1 but the external yield-loss coefficients were
#'   calibrated against response percentages from 0 to 100.
#' @param constraint_action Character. Action when the raw proportional
#'   disease-attributable loss falls outside the physically interpretable range
#'   `[0, 1]`:
#'
#'   - `"warn"` (default): constrain to `[0, 1]` and issue one warning;
#'   - `"error"`: stop before returning results.
#'
#'   Negative raw loss is normally prevented by the non-negative response and
#'   slope validation; values above 1 can occur when the supplied linear
#'   relationship extrapolates beyond complete yield loss.
#' @param seed Optional finite integer random seed. It controls only stochastic
#'   generation requested by this function (Uniform parameter ranges and/or the
#'   Gaussian intercept deviation). The caller's RNG state is restored.
#'
#' @return A data frame containing all original columns plus:
#'
#'   \describe{
#'     \item{`.loss_row_id`}{Original row index in `data`.}
#'     \item{`.sim`}{Loss-model simulation index.}
#'     \item{`slope_sim`, `intercept_sim`}{External yield-model parameters used
#'       in that simulation.}
#'     \item{`random_effect`}{Simulation-level external yield-intercept
#'       deviation. Retained under the historical column name for backward
#'       compatibility.}
#'     \item{`reference_yield`}{Zero-response yield for that loss simulation:
#'       `intercept_sim + random_effect`.}
#'     \item{`response_used`}{Scaled response `y * y_multiplier`.}
#'     \item{`yield_loss_proportion_raw`}{Unconstrained proportional loss
#'       `slope_sim * response_used / reference_yield`.}
#'     \item{`loss_constrained`}{Whether the raw proportional loss required
#'       constraining to `[0, 1]`.}
#'     \item{`yield_loss_proportion`, `yield_loss_pct`}{Constrained
#'       disease-attributable loss proportion and percent.}
#'     \item{`relative_yield_proportion`, `relative_yield_pct`}{Remaining yield
#'       relative to the simulation-specific reference yield. By construction,
#'       `relative_yield_proportion = 1 - yield_loss_proportion`.}
#'     \item{`predicted_yield`}{Yield on the external calibration-model scale:
#'       `reference_yield * relative_yield_proportion`.}
#'     \item{`attainable_yield`}{User-specified attainable yield in kg/ha.}
#'     \item{`attainable_predicted_yield`}{Attainable yield remaining after the
#'       same proportional loss.}
#'     \item{`yield_loss`}{Absolute attributable yield loss in kg/ha:
#'       `attainable_yield * yield_loss_proportion`.}
#'     \item{`price`}{Commodity price in USD per metric ton (USD/t).}
#'     \item{`economic_loss`}{Economic loss in USD/ha:
#'       `(yield_loss / 1000) * price`.}
#'   }
#'
#'   If `lower` and `upper` are supplied, corresponding `_lower` and `_upper`
#'   columns are added for every response-dependent metric.
#'
#'   Attributes document the loss-model specification, unit conversion,
#'   constraint behavior, parameter-simulation contract, and propagated-bound
#'   interpretation.
#'
#' @details
#' ## Loss equations
#'
#' The internally scaled response is
#'
#' \deqn{
#'   R = y \times y\_multiplier.
#' }
#'
#' For loss simulation \eqn{s}, the zero-response reference yield is
#'
#' \deqn{
#'   Y_{0,s} = \alpha_s + b_s,
#' }
#'
#' where \eqn{\alpha_s} is `intercept_sim` and \eqn{b_s} is the optional
#' simulation-level `random_effect`.
#'
#' The raw disease-attributable proportional loss is
#'
#' \deqn{
#'   L_{raw} =
#'   \frac{\beta_s R}{Y_{0,s}},
#' }
#'
#' where \eqn{\beta_s} is `slope_sim`.
#'
#' The physically interpretable loss proportion is
#'
#' \deqn{
#'   L = min\{1, max(0, L_{raw})\}.
#' }
#'
#' Consequently,
#'
#' \deqn{
#'   relative\_yield = 1 - L
#' }
#'
#' and the calibration-scale predicted yield is
#'
#' \deqn{
#'   predicted\_yield =
#'   Y_{0,s}(1-L).
#' }
#'
#' This construction deliberately keeps `predicted_yield`,
#' `relative_yield_proportion`, and `yield_loss_proportion` internally
#' consistent. Earlier code subtracted the simulated intercept deviation from
#' `predicted_yield` but omitted it from the proportional-loss denominator,
#' allowing those quantities to disagree.
#'
#' ## Attainable-yield scaling
#'
#' The proportional loss estimated from the external calibration relationship
#' is applied to each scenario-specific attainable yield:
#'
#' \deqn{
#'   absolute\_loss = attainable\_yield \times L,
#' }
#'
#' \deqn{
#'   attainable\_predicted\_yield =
#'   attainable\_yield \times (1-L).
#' }
#'
#' These two quantities sum exactly to `attainable_yield` up to floating-point
#' tolerance.
#'
#' ## Economic conversion
#'
#' The unit contract of `simulate_losses()` is fixed:
#'
#' - `attainable_yield` and `yield_loss` are in kg/ha;
#' - `price` is in USD per metric ton (USD/t);
#' - `economic_loss` is returned in USD/ha.
#'
#' Because one metric ton equals 1000 kg,
#'
#' \deqn{
#'   economic\_loss =
#'   \left(\frac{yield\_loss}{1000}\right)
#'   \times price.
#' }
#'
#' The division by 1000 is therefore part of the function's unit contract, not
#' a tunable modeling parameter.
#'
#' ## Parameter simulations
#'
#' `n` controls only simulations of the external loss relationship. Scalar
#' parameters are fixed. A numeric vector of length two is interpreted as a
#' Uniform range, preserving the historical EpiExposure interface. Longer
#' vectors of length `n` are treated as user-supplied paired parameter draws.
#'
#' When both `slope` and `intercept` are supplied as two-value ranges, they are
#' drawn independently. This does **not** preserve slope-intercept covariance
#' from an external fitted yield model. If joint parameter uncertainty is
#' available, provide paired user-supplied vectors of length `n` instead.
#'
#' Because a two-value numeric input is reserved for a Uniform range, exactly
#' two externally paired draws cannot be distinguished from a range under this
#' backward-compatible interface. Use `n > 2` for supplied paired draws if that
#' distinction matters.
#'
#' If `n > 1` but `slope`, `intercept`, and the `random_sd` component all
#' produce the same values in every simulation, `.sim` merely repeats identical
#' loss-model calculations. The function warns explicitly in this situation so
#' repeated rows are not mistaken for propagated uncertainty.
#'
#' ## Propagating response bounds
#'
#' When `lower` and `upper` are supplied, every loss simulation evaluates both
#' endpoints using the **same** slope, intercept, random-effect realization,
#' attainable yield, and price. Metric-specific minima and maxima are then
#' reported as `_lower` and `_upper`.
#'
#' These are propagated endpoint bounds. If `lower` and `upper` are marginal
#' prediction intervals and `n > 1` simultaneously represents uncertainty in an
#' external yield model, the resulting endpoint columns are not automatically a
#' calibrated joint probability interval. For fully joint Monte Carlo
#' propagation, supply draw-level response values as rows of `data` and preserve
#' their draw identifiers.
#'
#' ## Relationship to EpiExposure predictions
#'
#' `simulate_losses()` is a downstream transformation. It does not alter the
#' EpiExposure prediction contract. If `y` comes from `predict_outcomes()` or
#' `compare_predictions()`, the user should normally provide a response-scale
#' expected outcome that is on the same scale used to calibrate `slope` and
#' `intercept`.
#'
#' `random_sd` pertains only to the external yield-loss relationship. It does
#' not reactivate group-specific random effects from the epidemiological model.
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
    constraint_action = c("warn", "error"),
    seed = NULL
) {

  # ==========================================================================
  # HELPERS
  # ==========================================================================

  valid_name <- function(x) {
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

  validate_numeric_vector <- function(
    x,
    argument,
    nonnegative = FALSE,
    positive = FALSE,
    unique_required = FALSE
  ) {
    if (!is.numeric(x) ||
        !length(x) ||
        anyNA(x) ||
        any(!is.finite(x))) {
      stop(
        "`", argument,
        "` must contain only finite, non-missing numeric values.",
        call. = FALSE
      )
    }

    x <- as.numeric(x)

    if (nonnegative && any(x < 0)) {
      stop(
        "`", argument, "` cannot contain negative values.",
        call. = FALSE
      )
    }

    if (positive && any(x <= 0)) {
      stop(
        "`", argument, "` must contain only values greater than zero.",
        call. = FALSE
      )
    }

    if (unique_required && anyDuplicated(x)) {
      stop(
        "`", argument, "` must contain unique values. Duplicated values would ",
        "create duplicate simulation-grid rows.",
        call. = FALSE
      )
    }

    x
  }

  nearly_equal <- function(a, b, tolerance) {
    abs(a - b) <=
      tolerance *
      pmax(
        1,
        abs(a),
        abs(b)
      )
  }

  # ==========================================================================
  # 1. GENERAL VALIDATION
  # ==========================================================================

  if (!is.data.frame(data) ||
      !nrow(data)) {
    stop(
      "`data` must be a non-empty data.frame.",
      call. = FALSE
    )
  }

  if (!valid_name(y)) {
    stop(
      "`y` must be one non-empty character column name.",
      call. = FALSE
    )
  }

  if (!y %in% names(data)) {
    stop(
      "Response column '", y, "' was not found in `data`.",
      call. = FALSE
    )
  }

  if (!is.numeric(data[[y]]) ||
      anyNA(data[[y]]) ||
      any(!is.finite(data[[y]]))) {
    stop(
      "Response column '", y,
      "' must contain only finite, non-missing numeric values.",
      call. = FALSE
    )
  }

  if (!valid_integer_scalar(n) ||
      n < 1) {
    stop(
      "`n` must be one positive integer.",
      call. = FALSE
    )
  }
  n <- as.integer(n)

  if (!is.numeric(random_sd) ||
      length(random_sd) != 1L ||
      is.na(random_sd) ||
      !is.finite(random_sd) ||
      random_sd < 0) {
    stop(
      "`random_sd` must be one non-negative finite numeric value.",
      call. = FALSE
    )
  }

  if (!is.numeric(y_multiplier) ||
      length(y_multiplier) != 1L ||
      is.na(y_multiplier) ||
      !is.finite(y_multiplier) ||
      y_multiplier <= 0) {
    stop(
      "`y_multiplier` must be one positive finite numeric value.",
      call. = FALSE
    )
  }

  constraint_action <- match.arg(
    constraint_action
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

  # ==========================================================================
  # 2. RESPONSE BOUNDS
  # ==========================================================================

  has_bounds <- !is.null(lower) ||
    !is.null(upper)

  if (has_bounds) {
    if (is.null(lower) ||
        is.null(upper)) {
      stop(
        "`lower` and `upper` must be supplied together.",
        call. = FALSE
      )
    }

    if (!valid_name(lower)) {
      stop(
        "`lower` must be NULL or one non-empty character column name.",
        call. = FALSE
      )
    }

    if (!valid_name(upper)) {
      stop(
        "`upper` must be NULL or one non-empty character column name.",
        call. = FALSE
      )
    }

    if (identical(lower, upper) ||
        identical(lower, y) ||
        identical(upper, y)) {
      stop(
        "`y`, `lower`, and `upper` must identify distinct columns.",
        call. = FALSE
      )
    }

    if (!lower %in% names(data)) {
      stop(
        "Lower-bound column '", lower,
        "' was not found in `data`.",
        call. = FALSE
      )
    }

    if (!upper %in% names(data)) {
      stop(
        "Upper-bound column '", upper,
        "' was not found in `data`.",
        call. = FALSE
      )
    }

    for (column in c(lower, upper)) {
      if (!is.numeric(data[[column]]) ||
          anyNA(data[[column]]) ||
          any(!is.finite(data[[column]]))) {
        stop(
          "Bound column '", column,
          "' must contain only finite, non-missing numeric values.",
          call. = FALSE
        )
      }
    }

    if (any(
      data[[lower]] >
      data[[upper]]
    )) {
      stop(
        "Every response bound must satisfy `lower <= upper`.",
        call. = FALSE
      )
    }

    outside_bounds <-
      data[[y]] < data[[lower]] |
      data[[y]] > data[[upper]]

    if (any(outside_bounds)) {
      stop(
        "Every central response value in `y` must lie within its supplied ",
        "`lower` and `upper` bounds.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # 3. LOSS-MODEL PARAMETERS AND ECONOMIC GRID
  # ==========================================================================

  slope <- validate_numeric_vector(
    slope,
    "slope",
    nonnegative = TRUE
  )

  intercept <- validate_numeric_vector(
    intercept,
    "intercept",
    positive = TRUE
  )

  attainable_yield <- validate_numeric_vector(
    attainable_yield,
    "attainable_yield",
    nonnegative = TRUE,
    unique_required = TRUE
  )

  price <- validate_numeric_vector(
    price,
    "price",
    nonnegative = TRUE,
    unique_required = TRUE
  )

  valid_parameter_length <- function(
    x,
    parameter
  ) {
    allowed <- unique(
      c(
        1L,
        2L,
        n
      )
    )

    if (!length(x) %in%
        allowed) {
      stop(
        "`", parameter,
        "` must have length 1, 2, or `n` (",
        n, ").",
        call. = FALSE
      )
    }
  }

  valid_parameter_length(
    slope,
    "slope"
  )

  valid_parameter_length(
    intercept,
    "intercept"
  )

  # ==========================================================================
  # 4. RESPONSE SCALE VALIDATION
  # ==========================================================================

  response_used <-
    as.numeric(
      data[[y]]
    ) *
    y_multiplier

  if (anyNA(response_used) ||
      any(!is.finite(response_used))) {
    stop(
      "Scaling `y` by `y_multiplier` produced non-finite values.",
      call. = FALSE
    )
  }

  if (any(response_used < 0)) {
    stop(
      "`simulate_losses()` requires a non-negative response after applying ",
      "`y_multiplier`, because the supplied slope is interpreted as an ",
      "increasing disease/damage-to-yield-loss relationship.",
      call. = FALSE
    )
  }

  response_lower_used <- NULL
  response_upper_used <- NULL

  if (has_bounds) {
    response_lower_used <-
      as.numeric(
        data[[lower]]
      ) *
      y_multiplier

    response_upper_used <-
      as.numeric(
        data[[upper]]
      ) *
      y_multiplier

    if (anyNA(response_lower_used) ||
        anyNA(response_upper_used) ||
        any(!is.finite(response_lower_used)) ||
        any(!is.finite(response_upper_used))) {
      stop(
        "Scaling the response bounds by `y_multiplier` produced non-finite ",
        "values.",
        call. = FALSE
      )
    }

    if (any(response_lower_used < 0) ||
        any(response_upper_used < 0)) {
      stop(
        "`lower` and `upper` must remain non-negative after applying ",
        "`y_multiplier`.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # 5. CHECK OUTPUT-COLUMN CONFLICTS
  # ==========================================================================

  base_output_columns <- c(
    ".loss_row_id",
    ".sim",
    "slope_sim",
    "intercept_sim",
    "random_effect",
    "reference_yield",
    "response_used",
    "yield_loss_proportion_raw",
    "loss_constrained",
    "predicted_yield",
    "relative_yield_proportion",
    "relative_yield_pct",
    "yield_loss_proportion",
    "yield_loss_pct",
    "attainable_yield",
    "attainable_predicted_yield",
    "yield_loss",
    "price",
    "economic_loss"
  )

  bound_output_columns <- c(
    "response_lower_used",
    "response_upper_used",
    "yield_loss_proportion_raw_lower",
    "yield_loss_proportion_raw_upper",
    "loss_constrained_lower",
    "loss_constrained_upper",
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
    "attainable_predicted_yield_lower",
    "attainable_predicted_yield_upper",
    "yield_loss_lower",
    "yield_loss_upper",
    "economic_loss_lower",
    "economic_loss_upper"
  )

  new_columns <- if (has_bounds) {
    c(
      base_output_columns,
      bound_output_columns
    )
  } else {
    base_output_columns
  }

  conflicts <- intersect(
    names(data),
    new_columns
  )

  if (length(conflicts)) {
    stop(
      "The following output columns already exist in `data`: ",
      paste(
        conflicts,
        collapse = ", "
      ),
      ". Rename them before calling `simulate_losses()`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 6. RNG CONTRACT
  # ==========================================================================

  stochastic_parameters <-
    length(slope) == 2L ||
    length(intercept) == 2L ||
    random_sd > 0

  if (!is.null(seed)) {
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

    set.seed(seed)
  }

  # ==========================================================================
  # 7. LOSS-PARAMETER SIMULATIONS
  # ==========================================================================

  draw_parameter <- function(
    x,
    n_use,
    parameter
  ) {
    if (length(x) == 1L) {
      return(
        rep(
          as.numeric(x),
          n_use
        )
      )
    }

    if (length(x) == 2L) {
      return(
        stats::runif(
          n_use,
          min = min(x),
          max = max(x)
        )
      )
    }

    if (length(x) == n_use) {
      return(
        as.numeric(x)
      )
    }

    stop(
      "Internal error while resolving `",
      parameter,
      "`.",
      call. = FALSE
    )
  }

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

  if (anyNA(slope_sim) ||
      any(!is.finite(slope_sim)) ||
      any(slope_sim < 0)) {
    stop(
      "Generated/supplied slope simulations must be finite and non-negative.",
      call. = FALSE
    )
  }

  if (anyNA(intercept_sim) ||
      any(!is.finite(intercept_sim)) ||
      any(intercept_sim <= 0)) {
    stop(
      "Generated/supplied intercept simulations must be finite and positive.",
      call. = FALSE
    )
  }

  random_effect <- if (random_sd > 0) {
    stats::rnorm(
      n,
      mean = 0,
      sd = random_sd
    )
  } else {
    rep(
      0,
      n
    )
  }

  if (anyNA(random_effect) ||
      any(!is.finite(random_effect))) {
    stop(
      "The external yield-intercept deviations are non-finite.",
      call. = FALSE
    )
  }

  reference_yield_sim <-
    intercept_sim +
    random_effect

  if (any(reference_yield_sim <= 0) ||
      any(!is.finite(reference_yield_sim))) {
    bad <- which(
      reference_yield_sim <= 0 |
        !is.finite(reference_yield_sim)
    )

    stop(
      "At least one external yield-model simulation produced a non-positive ",
      "`reference_yield = intercept + random_effect` (simulation(s): ",
      paste(
        utils::head(
          bad,
          10L
        ),
        collapse = ", "
      ),
      if (length(bad) > 10L) ", ..." else "",
      "). Reduce `random_sd` or provide a yield model whose zero-response ",
      "reference yield remains positive.",
      call. = FALSE
    )
  }

  no_loss_parameter_variation <-
    length(unique(slope_sim)) == 1L &&
    length(unique(intercept_sim)) == 1L &&
    length(unique(random_effect)) == 1L

  if (n > 1L &&
      no_loss_parameter_variation) {
    warning(
      "`n > 1` was requested, but there is no variation across simulations in ",
      "`slope`, `intercept`, or the component controlled by `random_sd`. ",
      "All `.sim` replicates therefore use the same external yield-loss model ",
      "parameters and do not represent uncertainty propagation. Use `n = 1` ",
      "unless repeated identical simulation rows are intentionally required.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 8. CONSTRUCT LOSS-SIMULATION GRID
  # ==========================================================================

  simulation_grid <- expand.grid(
    attainable_yield =
      attainable_yield,
    price =
      price,
    .sim =
      seq_len(n),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  simulation_grid$slope_sim <-
    slope_sim[
      simulation_grid$.sim
    ]

  simulation_grid$intercept_sim <-
    intercept_sim[
      simulation_grid$.sim
    ]

  simulation_grid$random_effect <-
    random_effect[
      simulation_grid$.sim
    ]

  simulation_grid$reference_yield <-
    reference_yield_sim[
      simulation_grid$.sim
    ]

  # ==========================================================================
  # 9. EXPAND ORIGINAL DATA
  # ==========================================================================

  n_grid <- nrow(
    simulation_grid
  )

  original_index <- rep(
    seq_len(
      nrow(data)
    ),
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

  out$.loss_row_id <-
    original_index

  out$.sim <-
    simulation_grid$.sim[
      simulation_index
    ]

  out$slope_sim <-
    simulation_grid$slope_sim[
      simulation_index
    ]

  out$intercept_sim <-
    simulation_grid$intercept_sim[
      simulation_index
    ]

  out$random_effect <-
    simulation_grid$random_effect[
      simulation_index
    ]

  out$reference_yield <-
    simulation_grid$reference_yield[
      simulation_index
    ]

  out$attainable_yield <-
    simulation_grid$attainable_yield[
      simulation_index
    ]

  out$price <-
    simulation_grid$price[
      simulation_index
    ]

  out$response_used <-
    response_used[
      original_index
    ]

  if (has_bounds) {
    out$response_lower_used <-
      response_lower_used[
        original_index
      ]

    out$response_upper_used <-
      response_upper_used[
        original_index
      ]
  }

  # ==========================================================================
  # 10. LOSS-METRIC CALCULATOR
  # ==========================================================================

  calculate_metrics <- function(
    response,
    context = "central"
  ) {
    raw_loss_proportion <-
      (
        out$slope_sim *
          response
      ) /
      out$reference_yield

    if (anyNA(raw_loss_proportion) ||
        any(!is.finite(
          raw_loss_proportion
        ))) {
      stop(
        "The raw proportional-loss calculation produced non-finite values for ",
        context,
        " response values.",
        call. = FALSE
      )
    }

    outside_constraints <-
      raw_loss_proportion < 0 |
      raw_loss_proportion > 1

    if (any(outside_constraints) &&
        identical(
          constraint_action,
          "error"
        )) {
      raw_range <- range(
        raw_loss_proportion[
          outside_constraints
        ]
      )

      stop(
        "The external linear yield-loss relationship produced raw ",
        "proportional loss outside [0, 1] for ",
        sum(
          outside_constraints
        ),
        " expanded row(s) under ",
        context,
        " response values. Outside-range values span ",
        format(
          raw_range[1L],
          digits = 6
        ),
        " to ",
        format(
          raw_range[2L],
          digits = 6
        ),
        ". Use a scientifically appropriate response/parameter range or set ",
        "`constraint_action = 'warn'` to apply physical [0, 1] constraints.",
        call. = FALSE
      )
    }

    loss_proportion <- pmin(
      pmax(
        raw_loss_proportion,
        0
      ),
      1
    )

    relative_yield_proportion <-
      1 -
      loss_proportion

    predicted_yield <-
      out$reference_yield *
      relative_yield_proportion

    relative_yield_pct <-
      100 *
      relative_yield_proportion

    yield_loss_pct <-
      100 *
      loss_proportion

    attainable_predicted_yield <-
      out$attainable_yield *
      relative_yield_proportion

    yield_loss <-
      out$attainable_yield *
      loss_proportion

    economic_loss <-
      (
        yield_loss /
          1000
      ) *
      out$price

    if (anyNA(
      c(
        predicted_yield,
        relative_yield_proportion,
        relative_yield_pct,
        loss_proportion,
        yield_loss_pct,
        attainable_predicted_yield,
        yield_loss,
        economic_loss
      )
    ) ||
    any(!is.finite(
      c(
        predicted_yield,
        relative_yield_proportion,
        relative_yield_pct,
        loss_proportion,
        yield_loss_pct,
        attainable_predicted_yield,
        yield_loss,
        economic_loss
      )
    ))) {
      stop(
        "A derived loss metric became non-finite for ",
        context,
        " response values.",
        call. = FALSE
      )
    }

    list(
      yield_loss_proportion_raw =
        raw_loss_proportion,
      loss_constrained =
        outside_constraints,
      predicted_yield =
        predicted_yield,
      relative_yield_proportion =
        relative_yield_proportion,
      relative_yield_pct =
        relative_yield_pct,
      yield_loss_proportion =
        loss_proportion,
      yield_loss_pct =
        yield_loss_pct,
      attainable_predicted_yield =
        attainable_predicted_yield,
      yield_loss =
        yield_loss,
      economic_loss =
        economic_loss
    )
  }

  # ==========================================================================
  # 11. CENTRAL RESPONSE METRICS
  # ==========================================================================

  central_metrics <- calculate_metrics(
    out$response_used,
    context = "central"
  )

  for (metric_name in names(
    central_metrics
  )) {
    out[[metric_name]] <-
      central_metrics[[
        metric_name
      ]]
  }

  # ==========================================================================
  # 12. PROPAGATE RESPONSE BOUNDS
  # ==========================================================================

  lower_metrics <- NULL
  upper_metrics <- NULL

  if (has_bounds) {
    metrics_at_lower <- calculate_metrics(
      out$response_lower_used,
      context = "lower-bound"
    )

    metrics_at_upper <- calculate_metrics(
      out$response_upper_used,
      context = "upper-bound"
    )

    response_dependent_metrics <- c(
      "yield_loss_proportion_raw",
      "predicted_yield",
      "relative_yield_proportion",
      "relative_yield_pct",
      "yield_loss_proportion",
      "yield_loss_pct",
      "attainable_predicted_yield",
      "yield_loss",
      "economic_loss"
    )

    for (metric_name in response_dependent_metrics) {
      lower_values <- pmin(
        metrics_at_lower[[
          metric_name
        ]],
        metrics_at_upper[[
          metric_name
        ]]
      )

      upper_values <- pmax(
        metrics_at_lower[[
          metric_name
        ]],
        metrics_at_upper[[
          metric_name
        ]]
      )

      out[[
        paste0(
          metric_name,
          "_lower"
        )
      ]] <- lower_values

      out[[
        paste0(
          metric_name,
          "_upper"
        )
      ]] <- upper_values
    }

    out$loss_constrained_lower <-
      metrics_at_lower$loss_constrained

    out$loss_constrained_upper <-
      metrics_at_upper$loss_constrained

    lower_metrics <-
      metrics_at_lower
    upper_metrics <-
      metrics_at_upper
  }

  # ==========================================================================
  # 13. CONSISTENCY CHECKS
  # ==========================================================================

  tolerance <- 100 *
    .Machine$double.eps

  complement_check <-
    out$relative_yield_proportion +
    out$yield_loss_proportion

  if (any(
    !nearly_equal(
      complement_check,
      rep(
        1,
        length(
          complement_check
        )
      ),
      tolerance
    )
  )) {
    stop(
      "Internal error: relative yield and proportional yield loss are not ",
      "complements.",
      call. = FALSE
    )
  }

  attainable_check <-
    out$attainable_predicted_yield +
    out$yield_loss

  if (any(
    !nearly_equal(
      attainable_check,
      out$attainable_yield,
      tolerance
    )
  )) {
    stop(
      "Internal error: attainable predicted yield plus absolute loss does not ",
      "equal attainable yield.",
      call. = FALSE
    )
  }

  calibration_check <-
    out$predicted_yield -
    (
      out$reference_yield *
        out$relative_yield_proportion
    )

  if (any(
    abs(calibration_check) >
    tolerance *
    pmax(
      1,
      abs(
        out$predicted_yield
      ),
      abs(
        out$reference_yield
      )
    )
  )) {
    stop(
      "Internal error: calibration-scale predicted yield is inconsistent with ",
      "the proportional-loss calculation.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 14. EXPLICIT CONSTRAINT WARNING
  # ==========================================================================

  central_constrained_count <-
    sum(
      out$loss_constrained
    )

  lower_constrained_count <- if (
    has_bounds
  ) {
    sum(
      out$loss_constrained_lower
    )
  } else {
    0L
  }

  upper_constrained_count <- if (
    has_bounds
  ) {
    sum(
      out$loss_constrained_upper
    )
  } else {
    0L
  }

  any_constraints <-
    central_constrained_count > 0L ||
    lower_constrained_count > 0L ||
    upper_constrained_count > 0L

  if (any_constraints &&
      identical(
        constraint_action,
        "warn"
      )) {
    warning(
      "The external linear yield-loss relationship exceeded the physical ",
      "[0, 1] proportional-loss range and was constrained. Expanded-row ",
      "counts: central = ",
      central_constrained_count,
      if (has_bounds) {
        paste0(
          ", lower bound = ",
          lower_constrained_count,
          ", upper bound = ",
          upper_constrained_count
        )
      } else {
        ""
      },
      ". Inspect `yield_loss_proportion_raw` and the `loss_constrained` ",
      "columns before interpreting extrapolated scenarios.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 15. OUTPUT ORDER AND ATTRIBUTES
  # ==========================================================================

  rownames(out) <- NULL

  attr(
    out,
    "loss_spec"
  ) <- list(
    y = y,
    lower = lower,
    upper = upper,
    y_multiplier =
      y_multiplier,
    n =
      n,
    slope =
      slope,
    intercept =
      intercept,
    random_sd =
      random_sd,
    attainable_yield =
      attainable_yield,
    price =
      price,
    yield_unit =
      "kg_per_ha",
    price_unit =
      "USD_per_metric_ton",
    economic_loss_unit =
      "USD_per_ha",
    kg_per_metric_ton =
      1000,
    constraint_action =
      constraint_action,
    seed =
      seed,
    response_contract =
      "nonnegative_response_on_external_yield_model_scale",
    parameter_contract =
      "fixed_uniform_or_user_supplied_loss_model_draws",
    random_effect_contract =
      "simulation_level_external_yield_intercept_deviation",
    loss_contract =
      "attributable_loss_relative_to_simulation_specific_reference_yield",
    economic_contract =
      "yield_loss_kg_per_ha_divided_by_1000_times_USD_per_metric_ton",
    bounds_contract =
      if (has_bounds) {
        "deterministic_endpoint_propagation_not_joint_probability_interval"
      } else {
        NULL
      },
    cap_loss =
      TRUE,
    truncate_yield =
      TRUE,
    keep_row_id =
      TRUE,
    central_constrained_rows =
      central_constrained_count,
    lower_constrained_rows =
      lower_constrained_count,
    upper_constrained_rows =
      upper_constrained_count,
    parameter_variation_present =
      !no_loss_parameter_variation,
    repeated_identical_simulations =
      n > 1L &&
      no_loss_parameter_variation
  )

  attr(
    out,
    "epiexposure_loss_simulation_contract"
  ) <- "external_yield_loss_parameter_simulations"

  attr(
    out,
    "epiexposure_loss_n"
  ) <- n

  attr(
    out,
    "epiexposure_loss_bounds"
  ) <- if (has_bounds) {
    "endpoint_propagation"
  } else {
    "none"
  }

  attr(
    out,
    "epiexposure_loss_constraint_action"
  ) <- constraint_action

  attr(
    out,
    "epiexposure_loss_parameter_dependence"
  ) <- if (
    length(slope) == 2L ||
    length(intercept) == 2L
  ) {
    "uniform_ranges_sampled_independently_unless_user_supplies_paired_draws"
  } else {
    "user_supplied_or_fixed_parameters_paired_by_simulation_index"
  }

  attr(
    out,
    "epiexposure_loss_random_component"
  ) <- if (random_sd > 0) {
    "external_yield_intercept_deviation"
  } else {
    "none"
  }

  attr(
    out,
    "epiexposure_loss_parameter_variation"
  ) <- !no_loss_parameter_variation

  attr(
    out,
    "epiexposure_loss_repeated_identical_simulations"
  ) <- n > 1L &&
    no_loss_parameter_variation

  out
}
