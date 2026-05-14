#' Check DLNM cross-basis identifiability
#'
#' Evaluates whether a DLNM cross-basis matrix is full rank.
#' Rank deficiency may lead to non-identifiable or unstable model parameters.
#'
#' @param cb_template A crossbasis object
#'
#' @return Logical value indicating whether the cross-basis is identifiable.
#' @export
check_identifiability <- function(cb_template) {

  X <- as.matrix(cb_template)
  p <- ncol(X)
  k <- qr(X)$rank

  if (k < p) {
    warning(
      "Cross-basis is rank-deficient (", k, " < ", p,
      "). Consider reducing df_var, df_lag, or lag_max."
    )
    return(FALSE)
  }

  TRUE
}
