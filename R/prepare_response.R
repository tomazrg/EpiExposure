#' Prepare response variable for modeling
#'
#' Prepares and validates a response variable (`y`) according to the selected
#' distribution family. When needed, the function applies transformations to
#' ensure compatibility with the model assumptions (e.g., scaling proportions
#' for Beta regression or enforcing integer counts for Poisson models).
#'
#' If the response variable is already compatible with the selected distribution,
#' values are preserved without modification.
#'
#' @param data A data.frame containing the response variable.
#' @param y_var Character scalar. Name of the response variable in `data`.
#' @param family Character scalar. Distribution to be used in modeling.
#'   Supported options are:
#'   - `"beta"`: proportions in (0,1)
#'   - `"binomial"`: binary outcomes (0/1)
#'   - `"poisson"`: count data (non-negative integers)
#'   - `"negbin"` or `"negative_binomial"`: overdispersed counts
#'   - `"gaussian"`: continuous values
#'   - `"gamma"`: positive continuous values
#'   - `"ordinal"`: integer categorical severity levels (1, 2, 3, ...)
#'
#' @return A data.frame identical to `data`, with an added column:
#'   - `y_model`: processed response variable ready for modeling
#'
#'   Additionally, the function attaches an attribute:
#'   - `"y_scale_mult"`: scaling factor applied to `y` (e.g., 100 if percentages were converted to proportions)
#'
#' @details
#' - For `"beta"`, values are constrained to the open interval (0,1) using a small epsilon.
#' - For `"binomial"`, non-binary values are set to NA.
#' - For count families (`"poisson"`, `"negbin"`), values are rounded and negative values are set to NA.
#' - For `"gamma"`, non-positive values are set to NA.
#' - For `"ordinal"`, values must be integers ≥ 1.
#'
#' @export
prepare_response <- function(data, y_var, family) {

  # =========================================================
  # ✅ CHECKS BÁSICOS (AJUSTE DE ROBUSTEZ)
  # =========================================================
  stopifnot(
    y_var %in% names(data),
    is.character(family),
    length(family) == 1
  )

  y <- data[[y_var]]
  y_scale <- 1

  # -------------------------------
  # Beta (severity, prevalence)
  # -------------------------------
  if (family == "beta") {

    # ✅ AJUSTE: garantir tipo numérico
    y <- as.numeric(y)

    if (max(y, na.rm = TRUE) > 1) {
      y <- y / 100
      y_scale <- 100
    }

    eps <- 1e-5
    y <- pmin(pmax(y, eps), 1 - eps)
  }

  # -------------------------------
  # Binomial (incidence / prevalence)
  # -------------------------------
  if (family == "binomial") {

    y <- as.numeric(y)
    y <- as.integer(y)

    y[!y %in% c(0, 1)] <- NA
  }

  # -------------------------------
  # Poisson (counts)
  # -------------------------------
  if (family == "poisson") {

    y <- as.numeric(y)
    y <- round(y)

    y[y < 0] <- NA
  }

  # -------------------------------
  # Negative binomial (overdispersed counts)
  # -------------------------------
  if (family %in% c("negbin", "negative_binomial")) {

    y <- as.numeric(y)
    y <- round(y)

    y[y < 0] <- NA
  }

  # -------------------------------
  # Gaussian (continuous responses)
  # -------------------------------
  if (family == "gaussian") {

    y <- as.numeric(y)
  }

  # -------------------------------
  # Gamma (rate / intensity > 0)
  # -------------------------------
  if (family == "gamma") {

    y <- as.numeric(y)

    y[y <= 0] <- NA
  }

  # -------------------------------
  # Ordinal severity (cumulative logit)
  # -------------------------------
  if (family == "ordinal") {

    y <- as.numeric(y)
    y <- as.integer(y)

    y[y < 1] <- NA
  }

  # =========================================================
  # ✅ OUTPUT
  # =========================================================
  data$y_model <- y
  attr(data, "y_scale_mult") <- y_scale

  data
}
