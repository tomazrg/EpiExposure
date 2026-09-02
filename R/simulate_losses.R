#' Simulate yield and economic losses
#'
#' Simulates predicted yield, relative yield, yield loss, and economic loss
#' from a user-specified response variable. Each row of the original data is
#' retained and expanded across combinations of attainable yield, price, and
#' optional parameter simulations.
#'
#' @param data A data.frame containing the response variable and any additional
#'   scenario variables.
#' @param y Character string giving the name of the numeric response column
#'   used to calculate yield and economic losses.
#' @param slope Numeric. Regression slope relating `y` to yield loss.
#'   A single value uses a fixed slope. Two values are interpreted as the
#'   minimum and maximum of a uniform distribution. Alternatively, a vector
#'   of length `n_sim` can contain user-supplied parameter draws.
#' @param intercept Numeric. Yield-model intercept. A single value uses a
#'   fixed intercept. Two values are interpreted as the minimum and maximum
#'   of a uniform distribution. Alternatively, a vector of length `n_sim`
#'   can contain user-supplied parameter draws.
#' @param attainable_yield Numeric vector of attainable yield values.
#' @param price Numeric vector of commodity prices.
#' @param n_sim Number of parameter simulations. Default is 1.
#' @param random_sd Standard deviation of the normally distributed random
#'   effect used in the yield model. Default is 0 (no random effect).
#' @param seed Optional random seed.
#'
#' @return A data.frame containing all original columns from `data`, plus:
#'   `.loss_row_id`, `.sim`, `slope_sim`, `intercept_sim`, `random_effect`,
#'   `predicted_yield`, `relative_yield_pct`, `yield_loss_pct`,
#'   `attainable_yield`, `yield_loss`, `price`, and `economic_loss`.
#'
#' @details
#' Predicted yield is calculated as:
#'
#'   intercept - slope * y - random_effect
#'
#' Predicted yield is automatically constrained to be non-negative.
#'
#' Percentage yield loss is calculated as:
#'
#'   (slope / intercept) * y * 100
#'
#' and is automatically constrained between 0 and 100%.
#'
#' Absolute yield loss is:
#'
#'   attainable_yield * yield_loss_pct / 100
#'
#' Economic loss is:
#'
#'   (yield_loss / 1000) * price
#'
#' The division by 1000 is fixed internally, assuming yield is expressed in
#' kg per unit area and price is expressed per metric ton.
#'
#' @export
simulate_losses <- function(
    data,
    y,
    slope = 100,
    intercept = 11142.94,
    attainable_yield = seq(4000, 12000, by = 50),
    price = seq(100, 300, by = 10),
    n_sim = 1L,
    random_sd = 0,
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
    stop(
      "Column `", y, "` was not found in `data`."
    )
  }
  
  if (!is.numeric(data[[y]])) {
    stop("The response column `", y, "` must be numeric.")
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
  
  if (length(n_sim) != 1L ||
      !is.numeric(n_sim) ||
      !is.finite(n_sim) ||
      n_sim < 1 ||
      n_sim != as.integer(n_sim)) {
    stop("`n_sim` must be a positive integer.")
  }
  
  n_sim <- as.integer(n_sim)
  
  if (length(random_sd) != 1L ||
      !is.numeric(random_sd) ||
      !is.finite(random_sd) ||
      random_sd < 0) {
    stop("`random_sd` must be a single non-negative numeric value.")
  }
  
  
  # ============================================================
  # 2. Check output-column conflicts
  # ============================================================
  
  new_columns <- c(
    ".loss_row_id",
    ".sim",
    "slope_sim",
    "intercept_sim",
    "random_effect",
    "predicted_yield",
    "relative_yield_pct",
    "yield_loss_pct",
    "attainable_yield",
    "yield_loss",
    "price",
    "economic_loss"
  )
  
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
      "`", name, "` must have length 1, 2, or `n_sim` (",
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
    n_sim,
    "slope"
  )
  
  intercept_sim <- draw_parameter(
    intercept,
    n_sim,
    "intercept"
  )
  
  if (any(intercept_sim <= 0)) {
    stop("`intercept` values must be greater than zero.")
  }
  
  if (random_sd > 0) {
    
    random_effect <- stats::rnorm(
      n_sim,
      mean = 0,
      sd = random_sd
    )
    
  } else {
    
    random_effect <- rep(
      0,
      n_sim
    )
  }
  
  
  # ============================================================
  # 6. Construct simulation grid
  # ============================================================
  
  simulation_grid <- expand.grid(
    attainable_yield = attainable_yield,
    price = price,
    .sim = seq_len(n_sim),
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
  # Every row of the original data is retained and repeated for
  # each combination of attainable yield, price, and simulation.
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
  
  # Always retain the identifier of the original row
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
  # 8. Response variable
  # ============================================================
  
  response <- out[[y]]
  
  
  # ============================================================
  # 9. Predicted yield
  #
  # Original equation:
  #
  # yield = b0 - b1 * y - u_j
  #
  # Predicted yield is automatically truncated at zero.
  # ============================================================
  
  out$predicted_yield <-
    out$intercept_sim -
    out$slope_sim * response -
    out$random_effect
  
  out$predicted_yield <-
    pmax(
      out$predicted_yield,
      0
    )
  
  
  # ============================================================
  # 10. Relative yield
  #
  # Original equation:
  #
  # relative = yield * 100 / b0
  # ============================================================
  
  out$relative_yield_pct <-
    out$predicted_yield *
    100 /
    out$intercept_sim
  
  out$relative_yield_pct <-
    pmin(
      pmax(
        out$relative_yield_pct,
        0
      ),
      100
    )
  
  
  # ============================================================
  # 11. Percentage yield loss
  #
  # Based on:
  #
  # (b1 / b0) * y * 100
  #
  # Yield loss is automatically constrained between 0 and 100%.
  # ============================================================
  
  out$yield_loss_pct <-
    (
      out$slope_sim /
        out$intercept_sim
    ) *
    response *
    100
  
  out$yield_loss_pct <-
    pmin(
      pmax(
        out$yield_loss_pct,
        0
      ),
      100
    )
  
  
  # ============================================================
  # 12. Absolute yield loss
  # ============================================================
  
  out$yield_loss <-
    (
      out$yield_loss_pct /
        100
    ) *
    out$attainable_yield
  
  
  # ============================================================
  # 13. Economic loss
  #
  # The divisor of 1000 is fixed internally:
  #
  # economic_loss = (yield_loss / 1000) * price
  #
  # This assumes yield is expressed in kg per unit area and price
  # is expressed per metric ton.
  # ============================================================
  
  out$economic_loss <-
    (
      out$yield_loss /
        1000
    ) *
    out$price
  
  
  # ============================================================
  # 14. Restore row names
  # ============================================================
  
  rownames(out) <- NULL
  
  
  # ============================================================
  # 15. Attributes documenting simulation
  # ============================================================
  
  attr(out, "loss_spec") <- list(
    y = y,
    n_sim = n_sim,
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
