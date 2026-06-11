#' Build epidemic-level DLNM design matrix
#'
#' @param wx_long Weather data (long format, must contain epi_id and exposures)
#' @param cb_templates DLNM templates
#' @param lag_max Maximum lag
#' @param include_response Logical. If TRUE, include response variable (y) if available
#'
#' @return data.frame (design matrix at epidemic level)
#' @export
build_design <- function(wx_long, cb_templates, lag_max,
                                include_response = TRUE) {

  # =========================================================
  # ✅ AJUSTE 1 — DIMENSÃO TEMPORAL (lag_max)
  # =========================================================
  # 🔵 Garante que lag_max seja escalar (evita c(0, 85))
  lag_max <- as.integer(max(lag_max))

  # =========================================================
  # ✅ CHECKS BÁSICOS
  # =========================================================
  stopifnot("epi_id" %in% names(wx_long))
  stopifnot("dpp" %in% names(wx_long))  # 🔵 necessário para validação temporal

  for (v in names(cb_templates)) {
    if (!v %in% names(wx_long)) {
      stop(paste0(
        "Variable '", v, "' not found in data. Available variables: ",
        paste(names(wx_long), collapse = ", ")
      ))
    }
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # =========================================================
  # ✅ AJUSTE 6 (NOVO) — CHECK AUTOMÁTICO DE COBERTURA TEMPORAL
  # =========================================================
  # 🔵 Garante que cada epidemia tenha dias suficientes para o lag_max
  .check_lag_coverage <- function(dat, lag_max) {

    n_required <- lag_max + 1

    bad_ids <- dat |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(n_days = dplyr::n_distinct(dpp), .groups = "drop") |>
      dplyr::filter(n_days < n_required)

    if (nrow(bad_ids) > 0) {
      stop(
        paste0(
          "Some epidemics do not have enough temporal coverage for lag_max.\n",
          "Required days per epi_id: ", n_required, "\n",
          "Example problematic epi_id: ",
          paste(head(bad_ids$epi_id, 5), collapse = ", ")
        )
      )
    }
  }

  .check_lag_coverage(wx_long, lag_max)

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

    # ✅ salvar template real apenas uma vez
    if (is.null(cb_templates_used[[var_name]])) {
      cb_templates_used[[var_name]] <<- cb
    }

    # =========================================================
    # ✅ AJUSTE 4 — ORDEM DOS COEFICIENTES
    # =========================================================
    # 🔵 Usa estrutura da crossbasis (ordem segura)
    cn <- colnames(cb)

    if (is.null(cn)) {
      cn <- paste0("cb_", var_name, "_", seq_len(ncol(cb)))
    } else {
      cn <- paste0("cb_", var_name, "_", seq_len(length(cn)))
    }

    out <- as.numeric(cb[length(x), ])
    names(out) <- cn

    out
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

    tmp |>
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
