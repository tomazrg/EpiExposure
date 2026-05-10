#' Build epidemic-level DLNM design matrix
#'
#' @param wx_long Weather data
#' @param cb_templates DLNM templates
#' @param lag_max Maximum lag
#'
#' @return data.frame
#' @export
build_design_matrix <- function(wx_long, cb_templates, lag_max) {
  
  extract_last <- function(x, cb_template) {
    cb <- dlnm::crossbasis(
      x,
      lag    = lag_max,
      argvar = attr(cb_template, "argvar"),
      arglag = attr(cb_template, "arglag")
    )
    as.numeric(cb[length(x), ])
  }
  
  X_list <- purrr::imap(cb_templates, function(cb, v) {
    
    tmp <- wx_long |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(
        cb = list(extract_last(.data[[v]], cb)),
        .groups = "drop"
      )
    
    p  <- length(tmp$cb[[1]])
    nm <- paste0("cb_", v, "_", seq_len(p))
    
    tmp |>
      dplyr::mutate(cb = lapply(cb, setNames, nm)) |>
      tidyr::unnest_wider(cb)
  })
  
  purrr::reduce(X_list, dplyr::left_join, by = "epi_id")
}