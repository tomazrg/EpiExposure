#' Sensitivity analysis over degrees of freedom
#'
#' @param df_grid Vector of dfs
#' @param fit_fun Function that fits the model
#'
#' @return data.frame
#' @export
df_sensitivity <- function(df_grid, fit_fun) {
  
  purrr::map_dfr(df_grid, function(df) {
    fit <- fit_fun(df)
    data.frame(df = df, AIC = AIC(fit))
  })
}