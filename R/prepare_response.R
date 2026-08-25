#' Prepare response variable for modeling
#'
#' Prepares and validates a response variable according to the selected
#' distribution family. When necessary, the function transforms the response
#' to satisfy the family requirements, such as converting percentages to
#' proportions for Beta regression.
#'
#' The `family` argument may be supplied as a supported character name or as a
#' family object containing a valid `family` field. Family aliases are
#' normalized internally to a canonical name.
#'
#' @param data A data.frame containing the response variable.
#' @param y_var Character scalar naming the response variable in `data`.
#' @param family Distribution family supplied as a supported character name or
#'   family object. Canonical family names are:
#'   - `"beta"`: proportions in the open interval (0, 1)
#'   - `"binomial"`: binary outcomes coded as 0 or 1
#'   - `"poisson"`: non-negative integer counts
#'   - `"negative_binomial"`: overdispersed non-negative integer counts
#'   - `"gaussian"`: continuous numeric values
#'   - `"gamma"`: strictly positive continuous values
#'   - `"ordinal"`: ordered integer categories coded as 1, 2, 3, ...
#'
#'   Supported aliases include `"bernoulli"`, `"normal"`, `"negbin"`,
#'   `"nbinom"`, `"nbinom1"`, `"nbinom2"`, and `"cumulative"`.
#'
#' @return A data.frame identical to `data`, with an added `y_model` column.
#'   The returned data frame contains the following attributes:
#'   - `"y_scale_mult"`: scaling factor applied to the original response
#'   - `"response_family_name"`: canonical family name used for preparation
#'   - `"response_levels"`: ordered response levels for ordinal outcomes,
#'     otherwise `NULL`
#'   - `"response_n_categories"`: number of ordinal categories, otherwise
#'     `NULL`
#'
#' @details
#' For `"beta"`, responses expressed on a 0-100 percentage scale are divided
#' by 100 when at least one finite value exceeds 1. Values are then constrained
#' to the open interval (0, 1) using `epsilon = 1e-5`.
#'
#' For `"binomial"`, values must be exactly 0 or 1. For `"poisson"` and
#' `"negative_binomial"`, values must be non-negative integers. The function
#' does not silently truncate or round invalid binary or count responses.
#'
#' For `"gamma"`, values must be strictly positive. For `"ordinal"`, values
#' must be consecutive integers beginning at 1 and at least two categories must
#' be present.
#'
#' No unsupported family is silently converted to Gaussian.
#'
#' @export
prepare_response <- function(data, y_var, family) {

  # =========================================================
  # BASIC VALIDATION
  # =========================================================

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  if (!is.character(y_var) || length(y_var) != 1L ||
      is.na(y_var) || !nzchar(y_var)) {
    stop("`y_var` must be one non-empty column name.")
  }

  if (!y_var %in% names(data)) {
    stop("Column '", y_var, "' was not found in `data`.")
  }

  # =========================================================
  # RESOLVE CANONICAL FAMILY NAME
  # =========================================================

  resolve_family_name <- function(family_input) {
    family_raw <- NULL

    if (is.character(family_input) && length(family_input) == 1L &&
        !is.na(family_input) && nzchar(family_input)) {
      family_raw <- family_input
    } else if (is.list(family_input) && !is.null(family_input$family) &&
               length(family_input$family) >= 1L) {
      family_raw <- as.character(family_input$family[[1]])
    } else if (inherits(family_input, "family") &&
               !is.null(family_input$family)) {
      family_raw <- as.character(family_input$family[[1]])
    } else {
      stop(
        "Unsupported `family` specification. Provide a supported family name ",
        "or a family object containing a valid `family` field."
      )
    }

    normalized <- tolower(trimws(family_raw))
    normalized <- gsub("[[:space:]-]+", "_", normalized)
    normalized <- gsub("[^a-z0-9_]", "", normalized)

    if (normalized %in% c("beta", "beta_family", "beta_proportion")) {
      return("beta")
    }

    if (normalized %in% c("binomial", "bernoulli")) {
      return("binomial")
    }

    if (normalized == "poisson") {
      return("poisson")
    }

    if (normalized == "gamma") {
      return("gamma")
    }

    if (normalized %in% c("gaussian", "normal")) {
      return("gaussian")
    }

    if (
      normalized %in% c(
        "negbin", "nbinom", "nbinom1", "nbinom2",
        "negative_binomial", "negative_binomial_1",
        "negative_binomial_2"
      ) || grepl("negative.*binomial", normalized)
    ) {
      return("negative_binomial")
    }

    if (normalized %in% c("ordinal", "cumulative")) {
      return("ordinal")
    }

    stop(
      "Unsupported family: '", family_raw, "'. Supported canonical families ",
      "are: beta, binomial, poisson, gamma, gaussian, negative_binomial, ",
      "and ordinal."
    )
  }

  family_name <- resolve_family_name(family)

  # =========================================================
  # RESPONSE VALIDATION
  # =========================================================

  original_y <- data[[y_var]]

  if (is.factor(original_y)) {
    if (family_name != "ordinal") {
      stop(
        "The response variable is a factor. Factor responses are supported ",
        "only for `family = 'ordinal'`."
      )
    }

    if (is.ordered(original_y)) {
      y <- as.numeric(original_y)
    } else {
      stop(
        "For `family = 'ordinal'`, a factor response must be an ordered factor."
      )
    }
  } else {
    if (!is.numeric(original_y) && !is.integer(original_y)) {
      stop("The response variable must be numeric or integer.")
    }
    y <- as.numeric(original_y)
  }

  if (!length(y)) {
    stop("The response variable cannot be empty.")
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop("The response variable must contain only finite, non-missing values.")
  }

  y_scale <- 1
  response_levels <- NULL
  response_n_categories <- NULL

  # =========================================================
  # BETA
  # =========================================================

  if (family_name == "beta") {
    if (any(y < 0)) {
      stop("For `family = 'beta'`, response values cannot be negative.")
    }

    if (max(y) > 1) {
      if (max(y) > 100) {
        stop(
          "For `family = 'beta'`, values greater than 100 cannot be ",
          "interpreted as percentages."
        )
      }

      y <- y / 100
      y_scale <- 100
    }

    epsilon <- 1e-5
    y <- pmin(pmax(y, epsilon), 1 - epsilon)
  }

  # =========================================================
  # BINOMIAL
  # =========================================================

  if (family_name == "binomial") {
    invalid_binary <- !y %in% c(0, 1)

    if (any(invalid_binary)) {
      stop("For `family = 'binomial'`, all response values must be 0 or 1.")
    }

    y <- as.integer(y)
  }

  # =========================================================
  # POISSON AND NEGATIVE BINOMIAL
  # =========================================================

  if (family_name %in% c("poisson", "negative_binomial")) {
    if (any(y < 0)) {
      stop(
        "For count families, all response values must be non-negative."
      )
    }

    non_integer <- abs(y - round(y)) > sqrt(.Machine$double.eps)

    if (any(non_integer)) {
      stop(
        "For count families, all response values must be integers. ",
        "Values are not rounded automatically."
      )
    }

    y <- as.integer(round(y))
  }

  # =========================================================
  # GAUSSIAN
  # =========================================================

  if (family_name == "gaussian") {
    y <- as.numeric(y)
  }

  # =========================================================
  # GAMMA
  # =========================================================

  if (family_name == "gamma") {
    if (any(y <= 0)) {
      stop("For `family = 'gamma'`, all response values must be greater than 0.")
    }

    y <- as.numeric(y)
  }

  # =========================================================
  # ORDINAL
  # =========================================================

  if (family_name == "ordinal") {
    non_integer <- abs(y - round(y)) > sqrt(.Machine$double.eps)

    if (any(non_integer)) {
      stop("For `family = 'ordinal'`, all response values must be integers.")
    }

    y <- as.integer(round(y))

    if (any(y < 1L)) {
      stop(
        "For `family = 'ordinal'`, categories must be positive integers ",
        "beginning at 1."
      )
    }

    response_levels <- sort(unique(y))
    expected_levels <- seq_len(max(response_levels))

    if (!identical(response_levels, expected_levels)) {
      stop(
        "For `family = 'ordinal'`, categories must be consecutive integers ",
        "beginning at 1."
      )
    }

    if (length(response_levels) < 2L) {
      stop("For `family = 'ordinal'`, at least two categories are required.")
    }

    response_n_categories <- length(response_levels)
  }

  # =========================================================
  # OUTPUT
  # =========================================================

  data$y_model <- y
  attr(data, "y_scale_mult") <- y_scale
  attr(data, "response_family_name") <- family_name
  attr(data, "response_levels") <- response_levels
  attr(data, "response_n_categories") <- response_n_categories

  data
}
