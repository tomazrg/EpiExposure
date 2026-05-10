#' Prepare response variable for modeling
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
  
  if (family_choice == "beta") {
    if (max(y, na.rm = TRUE) > 1) {
      y <- y / 100
      y_scale <- 100
    }
    eps <- 1e-5
    y <- pmin(pmax(y, eps), 1 - eps)
  }
  
  if (family_choice == "gamma") {
    y[y <= 0] <- NA
  }
  
  if (family_choice == "poisson") {
    y <- round(y)
    y[y < 0] <- NA
  }
  
  dat$y_model <- y
  attr(dat, "y_scale_mult") <- y_scale
  dat
}