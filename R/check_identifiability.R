#' Check DLNM cross-basis identifiability
#'
#' Evaluates whether a DLNM cross-basis matrix is full rank.
#'
#' @param data Long-format exposure data.
#'
#' Must contain:
#' - epi_id
#' - time
#'
#' Observations must be supplied in chronological order
#' (earliest observation → most recent observation).
#' @param var Exposure variable (e.g. "tmax")
#' @param max_lag Maximum lag
#' @param df_var Degrees of freedom (exposure)
#' @param df_lag Degrees of freedom (lag)
#' @param fun_var Basis ("ns","bs","poly","lin")
#' @param fun_lag Basis ("ns","ps","lin")
#'
#' @return TRUE/FALSE
#' @export
check_identifiability <- function(
    data,
    var,
    max_lag,
    df_var = 4,
    df_lag = 4,
    fun_var = "ns",
    fun_lag = "ns"
) {

  # =========================================================
  # AJUSTE 1 — DIMENSÃO TEMPORAL (max_lag)
  # =========================================================
  # 🔵 Garante consistência (evita c(0,85))
  max_lag <- as.integer(max(max_lag))

  # =========================================================
  # ✅ CHECKS BÁSICOS
  # =========================================================
  stopifnot(
    "epi_id" %in% names(data),
    "time" %in% names(data),
    var %in% names(data),
    fun_var %in% c("ns", "bs", "poly", "lin"),
    fun_lag %in% c("ns", "ps", "lin")
  )

  # =========================================================
  # ✅ AJUSTE 6 — CHECK AUTOMÁTICO DE COBERTURA TEMPORAL
  # =========================================================
  .check_lag_coverage <- function(dat, max_lag) {

    n_required <- max_lag + 1L

    temporal_summary <- dat |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        n_time = dplyr::n_distinct(time),
        has_missing_time = any(is.na(time)),
        .groups = "drop"
      )

    missing_time_ids <- temporal_summary |>
      dplyr::filter(has_missing_time)

    if (nrow(missing_time_ids) > 0) {

      stop(
        paste0(
          "Missing values were detected in `time`.\n",
          "Example problematic epi_id: ",
          paste(
            utils::head(missing_time_ids$epi_id, 5),
            collapse = ", "
          )
        )
      )
    }

    duplicated_time_ids <- temporal_summary |>
      dplyr::filter(n_rows != n_time)

    if (nrow(duplicated_time_ids) > 0) {

      stop(
        paste0(
          "Duplicated `time` values were detected within some epidemics.\n",
          "Each epidemic must contain one observation per time value.\n",
          "Example problematic epi_id: ",
          paste(
            utils::head(duplicated_time_ids$epi_id, 5),
            collapse = ", "
          )
        )
      )
    }

    bad_ids <- temporal_summary |>
      dplyr::filter(n_time < n_required)

    if (nrow(bad_ids) > 0) {

      stop(
        paste0(
          "Some epidemics do not have enough temporal coverage for max_lag.\n",
          "Required observations per epi_id: ",
          n_required,
          "\n",
          "Example problematic epi_id: ",
          paste(
            utils::head(bad_ids$epi_id, 5),
            collapse = ", "
          )
        )
      )
    }

    invisible(TRUE)
  }

  .check_lag_coverage(data, max_lag)

  # ----------------------------------------------------------
  # pooled series
  # ----------------------------------------------------------
  build_pooled_series <- function(dat, var, sep_n) {
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

  SEPARATOR <- max_lag
  x_pool <- build_pooled_series(data, var, SEPARATOR)

  # ----------------------------------------------------------
  # reconstruir crossbasis
  # ----------------------------------------------------------
  argvar <- switch(
    fun_var,
    ns   = list(fun = "ns", df = df_var),
    bs   = list(fun = "bs", df = df_var),
    poly = list(fun = "poly", degree = df_var),
    lin  = list(fun = "lin")
  )

  if (!is.null(argvar$fun) && argvar$fun != "lin") {
    argvar$intercept <- FALSE
  }

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
  # ✅ AJUSTE 4 — ORDEM DOS COEFICIENTES (garantia estrutural)
  # =========================================================
  # 🔵 armazenado para consistência/debug (não altera comportamento)
  cb_colnames <- colnames(cb)

  # ----------------------------------------------------------
  # rank check
  # ----------------------------------------------------------
  X <- as.matrix(cb)

  # =========================================================
  # ✅ AJUSTE CRÍTICO — REMOVER NA PARA RANK CORRETO
  # =========================================================
  # 🔴 evita rank falso devido aos NA separators
  X <- X[complete.cases(X), , drop = FALSE]

  p <- ncol(X)
  k <- qr(X)$rank

  if (k < p) {
    warning(
      paste0(
        "Cross-basis is rank-deficient (", k, " < ", p, "). ",
        "Consider reducing df_var, df_lag, or max_lag."
      )
    )
    return(FALSE)
  }

  return(TRUE)
}
