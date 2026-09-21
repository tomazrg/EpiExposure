#' Simulated epidemiological dataset
#'
#' A simulated dataset containing multiple epidemics and environmental
#' exposure variables for demonstrating the modeling functions available
#' in EpiExposure.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure histories and the following variables:
#' \describe{
#'   \item{epi_id}{Unique identifier for each epidemic.}
#'   \item{time}{Chronological time index within each epidemic history.}
#'   \item{tmean}{Mean temperature.}
#'   \item{rain}{Rainfall.}
#'   \item{wetness}{Leaf wetness or wetness-related exposure.}
#'   \item{y}{Simulated disease response.}
#' }
#'
#' @source Simulated data.
"epi_data"
