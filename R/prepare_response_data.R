#' Prepare response variable for modeling
#'
#' If the response variable is already compatible with the selected distribution,
#' the values are preserved without modification.
#'
#' @param dat Data frame
#' @param y_var Response variable
#' @param family_choice Distribution
#'
#' @return data.frame
#' @export
prepare_response_data <- function(dat, y_var, family_choice) {

  y <- dat[[y_var]]
  y_scale <- 1

  # -------------------------------
  # Beta (severity, prevalence)
  # -------------------------------
  if (family_choice == "beta") {

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
  if (family_choice == "binomial") {
    y <- as.integer(y)
    y[!y %in% c(0, 1)] <- NA
  }

  # -------------------------------
  # Poisson (counts)
  # -------------------------------
  if (family_choice == "poisson") {
    y <- round(y)
    y[y < 0] <- NA
  }

  # -------------------------------
  # Negative binomial (overdispersed counts)
  # -------------------------------
  if (family_choice %in% c("negbin", "negative_binomial")) {
    y <- round(y)
    y[y < 0] <- NA
  }

  # -------------------------------
  # Gaussian (continuous responses)
  # -------------------------------
  if (family_choice == "gaussian") {
    y <- as.numeric(y)
  }

  # -------------------------------
  # Gamma (rate / intensity > 0)
  # -------------------------------
  if (family_choice == "gamma") {
    y[y <= 0] <- NA
  }

  # -------------------------------
  # Ordinal severity (cumulative logit)
  # -------------------------------
  if (family_choice == "ordinal") {
    y <- as.integer(y)
    y[y < 1] <- NA
  }

  dat$y_model <- y
  attr(dat, "y_scale_mult") <- y_scale

  dat
}
