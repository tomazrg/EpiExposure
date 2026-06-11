#' Check DLNM cross-basis identifiability
#'
#' Evaluates whether a DLNM cross-basis matrix is full rank.
#'
#' @param wx_long Long-format weather data
#' @param var Exposure variable (e.g. "tmax")
#' @param lag_max Maximum lag
#' @param df_var Degrees of freedom (exposure)
#' @param df_lag Degrees of freedom (lag)
#' @param fun_var Basis ("ns","bs","poly","lin")
#' @param fun_lag Basis ("ns","ps","lin")
#'
#' @return TRUE/FALSE
#' @export
check_identifiability <- function(
    wx_long,
    var,
    lag_max,
    df_var = 4,
    df_lag = 4,
    fun_var = "ns",
    fun_lag = "ns"
) {

  # =========================================================
  # ✅ AJUSTE 1 — DIMENSÃO TEMPORAL (lag_max)
  # =========================================================
  # 🔵 Garante consistência (evita c(0,85))
  lag_max <- as.integer(max(lag_max))

  # =========================================================
  # ✅ CHECKS BÁSICOS
  # =========================================================
  stopifnot(
    "epi_id" %in% names(wx_long),
    "dpp" %in% names(wx_long),
    var %in% names(wx_long),
    fun_var %in% c("ns", "bs", "poly", "lin"),
    fun_lag %in% c("ns", "ps", "lin")
  )

  # =========================================================
  # ✅ AJUSTE 6 — CHECK AUTOMÁTICO DE COBERTURA TEMPORAL
  # =========================================================
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

  # ----------------------------------------------------------
  # pooled series
  # ----------------------------------------------------------
  build_pooled_series <- function(dat, var, sep_n) {
    ids <- unique(dat$epi_id)
    out <- vector("list", length(ids))

    for (i in seq_along(ids)) {
      v <- dat |>
        dplyr::filter(epi_id == ids[i]) |>
        dplyr::arrange(dpp) |>
        dplyr::pull(.data[[var]])

      out[[i]] <- c(v, rep(NA_real_, sep_n))
    }

    unlist(out)
  }

  SEPARATOR <- lag_max
  x_pool <- build_pooled_series(wx_long, var, SEPARATOR)

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
    lag    = lag_max,
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
        "Consider reducing df_var, df_lag, or lag_max."
      )
    )
    return(FALSE)
  }

  return(TRUE)
}
