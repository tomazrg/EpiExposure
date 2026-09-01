#' Define DLNM exposure templates with flexible splines
#'
#' Creates DLNM exposure templates from long-format exposure histories.
#'
#' The input dataset must contain:
#' - `epi_id`: epidemic identifier
#' - `time`: time index
#'
#' Exposure histories must be supplied in chronological order,
#' from the earliest observation to the most recent observation
#' prior to disease assessment.
#'
#' Internal conversion to the lag structure required by DLNMs is
#' handled automatically through the cross-basis representation.
#'
#' @param data Long-format exposure dataset.
#' @param vars Character vector of exposure variables.
#' @param max_lag Maximum lag.
#' @param df_var Degrees of freedom for the exposure dimension.
#' @param df_lag Degrees of freedom for the lag dimension.
#' @param fun_var Basis function for exposure ("ns", "bs", "poly", "lin").
#' @param fun_lag Basis function for lag ("ns", "ps", "lin").
#'
#' @return A named list of crossbasis templates with attribute `"spec"`.
#'
#' @export
define_exposure <- function(data, vars,
                            max_lag,
                            df_var = 4,
                            df_lag = 4,
                            fun_var = "ns",
                            fun_lag = "ns") {

  # =========================================================
  #  AJUSTE 1 — DIMENSÃO TEMPORAL (max_lag)
  # =========================================================
  #  Garante que max_lag seja escalar consistente
  max_lag <- as.integer(max(max_lag))

  # =========================================================
  # CHECKS BÁSICOS
  # =========================================================
  stopifnot(
    all(c("epi_id", "time") %in% names(data)),
    fun_var %in% c("ns", "bs", "poly", "lin"),
    fun_lag %in% c("ns", "ps", "lin")
  )

  for (v in vars) {
    stopifnot(v %in% names(data))
  }

  # =========================================================
  # JUSTE 6 — CHECK AUTOMÁTICO DE COBERTURA TEMPORAL
  # =========================================================
  # Evita crossbasis mal definida (mesmo raciocínio da build_design)
  .check_lag_coverage <- function(dat, max_lag) {

    n_required <- max_lag + 1

    bad_ids <- dat |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(n_days = dplyr::n_distinct(time), .groups = "drop") |>
      dplyr::filter(n_days < n_required)

    if (nrow(bad_ids) > 0) {
      stop(
        paste0(
          "Some epidemics do not have enough temporal coverage for max_lag.\n",
          "Required days per epi_id: ", n_required, "\n",
          "Example problematic epi_id: ",
          paste(head(bad_ids$epi_id, 5), collapse = ", ")
        )
      )
    }
  }

  .check_lag_coverage(data, max_lag)

  # =========================================================
  # ooled series
  # =========================================================
  build_pooled <- function(dat, var, sep_n) {

    ids <- unique(dat$epi_id)
    out <- vector("list", length(ids))

    for (i in seq_along(ids)) {
      v <- dat |>
        dplyr::filter(epi_id == ids[i]) |>
        dplyr::arrange(time) |>
        dplyr::pull(.data[[var]])

      out[[i]] <- c(v, rep(NA_real_, sep_n))
    }

    unlist(out)
  }

  cb_templates <- list()

  for (v in vars) {

    x_pool <- build_pooled(data, v, max_lag)

    # =========================================================
    # argvar: exposure
    # =========================================================
    argvar <- switch(
      fun_var,
      ns   = list(fun = "ns",  df = df_var),
      bs   = list(fun = "bs",  df = df_var),
      poly = list(fun = "poly", degree = df_var),
      lin  = list(fun = "lin")
    )

    # 🔵 evita intercepto na dimensão de exposição
    if (!is.null(argvar$fun) && argvar$fun != "lin") {
      argvar$intercept <- FALSE
    }

    # =========================================================
    # arglag: lag
    # =========================================================
    arglag <- switch(
      fun_lag,
      ns  = list(fun = "ns", df = df_lag),
      ps  = list(fun = "ps", df = df_lag),
      lin = list(fun = "lin")
    )

    cb <- dlnm::crossbasis(
      x_pool,
      lag    = max_lag,
      argvar = argvar,
      arglag = arglag
    )

    # =========================================================
    #  AJUSTE 4 — ORDEM DOS COEFICIENTES (ESTRUTURAL)
    # =========================================================
    # 🔵 garante indexação consistente futura (evita ambiguidades)
    attr(cb, "cb_colnames") <- colnames(cb)

    cb_templates[[v]] <- cb
  }

  # =========================================================
  #  criar spec automaticamente
  # =========================================================
  epiexposure_spec <- lapply(cb_templates, function(cb) {

    # 🔵 usa SEMPRE max() para evitar c(0,85)
    lag_val <- attr(cb, "lag")
    lag_val <- if (length(lag_val) > 1) max(lag_val) else lag_val

    list(
      max_lag = lag_val,
      argvar  = attr(cb, "argvar"),
      arglag  = attr(cb, "arglag")
    )
  })

  # =========================================================
  #  anexar spec
  # =========================================================
  attr(cb_templates, "spec") <- epiexposure_spec

  return(cb_templates)
}
