#' Build epidemic-level DLNM design matrix
#'
#' @param wx_long Weather data (long format, must contain epi_id and exposures)
#' @param cb_templates DLNM templates
#' @param lag_max Maximum lag
#' @param include_response Logical. If TRUE, include response variable (y) if available
#'
#' @return data.frame (design matrix at epidemic level)
#' @export
build_design_matrix <- function(wx_long, cb_templates, lag_max,
                                include_response = TRUE) {

  # =========================================================
  # ✅ CHECKS
  # =========================================================

  stopifnot("epi_id" %in% names(wx_long))

  for (v in names(cb_templates)) {
    if (!v %in% names(wx_long)) {
      stop(paste0(
        "Variable '", v, "' not found in data. Available variables: ",
        paste(names(wx_long), collapse = ", ")
      ))
    }
  }

  # helper
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # =========================================================
  # ✅ armazenar templates reais usados
  # =========================================================
  cb_templates_used <- list()

  extract_last_cb_row <- function(x, cb_template, var_name) {

    cb <- dlnm::crossbasis(
      x,
      lag    = lag_max,
      argvar = attr(cb_template, "argvar"),
      arglag = attr(cb_template, "arglag")
    )

    # salvar template real apenas uma vez
    if (is.null(cb_templates_used[[var_name]])) {
      cb_templates_used[[var_name]] <<- cb
    }

    as.numeric(cb[length(x), ])
  }

  # =========================================================
  # ✅ BUILD DESIGN MATRIX (X)
  # =========================================================

  X_list <- purrr::imap(cb_templates, function(cb, v) {

    tmp <- wx_long |>
      dplyr::group_by(epi_id) |>
      dplyr::arrange(dpp, .by_group = TRUE) |>
      dplyr::summarise(
        cb = list(extract_last_cb_row(.data[[v]], cb, v)),
        .groups = "drop"
      )

    p  <- length(tmp$cb[[1]])
    nm <- paste0("cb_", v, "_", seq_len(p))

    tmp |>
      dplyr::mutate(cb = lapply(cb, setNames, nm)) |>
      tidyr::unnest_wider(cb)
  })

  out <- purrr::reduce(X_list, dplyr::left_join, by = "epi_id")

  # =========================================================
  # ✅ INCLUIR RESPOSTA (y) - OPCIONAL
  # =========================================================

  if (isTRUE(include_response) && "y" %in% names(wx_long)) {

    y_df <- wx_long |>
      dplyr::distinct(epi_id, y)

    # 🔴 validação crítica: 1 y por epidemia
    if (any(duplicated(y_df$epi_id))) {
      stop("Multiple 'y' values per epi_id detected. Expected one y per epidemic.")
    }

    # 🔴 validação: garantir correspondência completa
    if (nrow(y_df) != length(unique(wx_long$epi_id))) {
      stop("Mismatch between epi_id and y values.")
    }

    out <- dplyr::left_join(out, y_df, by = "epi_id")
  }

  # =========================================================
  # ✅ guardar templates utilizados
  # =========================================================

  attr(out, "cb_templates") <- cb_templates_used

  return(out)
}
