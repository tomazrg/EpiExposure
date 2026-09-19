#' Summarize an EpiExposure GAMM fit
#'
#' Summarizes a generalized additive mixed model fitted through
#' `fit_epidlnm()` with `model_engine = "gamm"`.
#'
#' The GAMM object contains separate `gam` and `lme` components.
#' By default, the method returns the `gam` summary containing the
#' fixed/population coefficients and their inferential statistics.
#'
#' @param object A fitted EpiExposure GAMM object.
#' @param component Character. Component to summarize:
#'
#'   - `"gam"` returns the fixed/population model summary;
#'   - `"lme"` returns the mixed-model summary;
#'   - `"both"` returns both summaries.
#'
#' @param ... Additional arguments passed to the corresponding
#'   summary method.
#'
#' @return A model summary object. When `component = "both"`, a list
#'   containing the GAM and LME summaries is returned.
#'
#' @export
summary.epiexposure_gamm <- function(
    object,
    component = c(
      "gam",
      "lme",
      "both"
    ),
    ...
) {
  
  component <- match.arg(
    component
  )
  
  if (!is.list(object) ||
      is.null(object$gam) ||
      is.null(object$lme)) {
    stop(
      paste0(
        "Invalid EpiExposure GAMM object: expected both ",
        "`$gam` and `$lme` components."
      ),
      call. = FALSE
    )
  }
  
  if (identical(component, "gam")) {
    return(
      summary(
        object$gam,
        ...
      )
    )
  }
  
  if (identical(component, "lme")) {
    return(
      summary(
        object$lme,
        ...
      )
    )
  }
  
  structure(
    list(
      gam = summary(
        object$gam,
        ...
      ),
      lme = summary(
        object$lme,
        ...
      )
    ),
    class = "summary_epiexposure_gamm"
  )
}


#' Print an EpiExposure GAMM summary
#'
#' @param x An object returned by
#'   `summary.epiexposure_gamm(component = "both")`.
#' @param ... Additional arguments passed to the component print methods.
#'
#' @return The summary object, invisibly.
#'
#' @export
print.summary_epiexposure_gamm <- function(
    x,
    ...
) {
  
  cat(
    "\nEpiExposure GAMM\n"
  )
  
  cat(
    "\nFixed/population component:\n\n"
  )
  
  print(
    x$gam,
    ...
  )
  
  cat(
    "\nMixed-model component:\n\n"
  )
  
  print(
    x$lme,
    ...
  )
  
  invisible(
    x
  )
}