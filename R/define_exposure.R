#' Define DLNM exposure templates with flexible splines
#'
#' @param wx_long Long-format weather data (epi_id, dpp, variables)
#' @param vars Character vector of exposure variables
#' @param lag_max Maximum lag
#' @param df_var Degrees of freedom for exposure dimension
#' @param df_lag Degrees of freedom for lag dimension
#' @param fun_var Basis function for exposure ("ns","bs","poly","lin")
#' @param fun_lag Basis function for lag ("ns","ps","lin")
#'
#' @return Named list of crossbasis templates (with attribute "spec")
#' @export
define_exposure <- function(wx_long, vars,
                            lag_max,
                            df_var = 4,
                            df_lag = 4,
                            fun_var = "ns",
                            fun_lag = "ns") {


  # ---- basic input checks ----
  stopifnot(
    all(c("epi_id", "dpp") %in% names(wx_long)),
    fun_var %in% c("ns", "bs", "poly", "lin"),
    fun_lag %in% c("ns", "ps", "lin")
  )

  # ---- check exposure variables ----
  for (v in vars) {
    stopifnot(v %in% names(wx_long))
  }

  # ----- pooled series  -----
  build_pooled <- function(dat, var, sep_n) {
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

  cb_templates <- list()

  for (v in vars) {

    x_pool <- build_pooled(wx_long, v, lag_max)

    # ----- argvar: exposure -----
    argvar <- switch(
      fun_var,
      ns   = list(fun = "ns",  df = df_var),
      bs   = list(fun = "bs",  df = df_var),
      poly = list(fun = "poly", degree = df_var),
      lin  = list(fun = "lin")
    )

    # evitar intercepto na dimensão da exposição
    if (!is.null(argvar$fun) && argvar$fun != "lin") {
      argvar$intercept <- FALSE
    }

    # ----- arglag: lag -----
    arglag <- switch(
      fun_lag,
      ns  = list(fun = "ns", df = df_lag),
      ps  = list(fun = "ps", df = df_lag),
      lin = list(fun = "lin")
    )

    cb_templates[[v]] <- dlnm::crossbasis(
      x_pool,
      lag    = lag_max,
      argvar = argvar,
      arglag = arglag
    )
  }

  # ----------------------------
  # ✅ NOVO: criar spec automaticamente
  # ----------------------------
  epiexposure_spec <- lapply(cb_templates, function(cb) {
    list(
      lag_max = attr(cb, "lag"),
      argvar  = attr(cb, "argvar"),
      arglag  = attr(cb, "arglag")
    )
  })

  # ----------------------------
  # ✅ NOVO: anexar como atributo (SEM quebrar nada)
  # ----------------------------
  attr(cb_templates, "spec") <- epiexposure_spec

  return(cb_templates)
}
