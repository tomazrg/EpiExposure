#' Read tabular data from csv, txt or Excel#' path File path
#' @param name File name (used to infer extension)
#'
#' @return data.frame
#' @export
read_any_table <- function(path, name) {
  
  ext <- tolower(tools::file_ext(name))
  
  if (ext %in% c("csv", "txt")) {
    readr::read_delim(path, delim = NULL, show_col_types = FALSE)
  } else if (ext %in% c("xls", "xlsx")) {
    readxl::read_excel(path)
  } else {
    stop("Unsupported file format")
  }
}
#'
