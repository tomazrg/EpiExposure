#' Check DLNM identifiability issues
#'
#' @param cb_template Crossbasis object
#'
#' @return logical
#' @export
check_identifiability <- function(cb_template) {
  
  X <- as.matrix(cb_template)
  k <- qr(X)$rank
  p <- ncol(X)
  
  if (k < p) {
    warning("Crossbasis is rank-deficient (potential identifiability issues).")
    return(FALSE)
  }
  
  TRUE
}