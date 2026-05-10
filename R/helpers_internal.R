get_cb_names <- function(dat, var) {
  grep)
}

make_series_band <- function(value, cen, n_ref, band = NULL) {
  x <- rep(cen, n_ref)
  if (!is.null(band)) {
    idx <- n_ref - band
    x[idx] <- value
  }
  x
}