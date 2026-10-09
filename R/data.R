#' Simulated beta epidemiological dataset
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles and a continuous plant disease outcome expressed as a proportion.
#' The dataset is intended for demonstrating beta regression and the main
#' modeling, interpretation, simulation, and prediction workflows available
#' in EpiExposure.
#'
#' The response is repeated across the temporal rows belonging to the same
#' epidemic because each complete exposure profile is associated with one final
#' epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique numeric identifier for each epidemic.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{y}{Continuous disease response expressed as a proportion strictly
#'   between 0 and 1 and repeated across the rows of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("epi_data", package = "EpiExposure")
#' head(epi_data)
"epi_data"


#' Simulated binomial epidemiological dataset
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles and a binary plant disease outcome. The dataset is intended for
#' demonstrating binomial models in EpiExposure.
#'
#' The response is repeated across the temporal rows belonging to the same
#' epidemic because each complete exposure profile is associated with one final
#' epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique numeric identifier for each epidemic.}
#'   \item{block}{Experimental block identifier.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{y}{Binary disease outcome coded as 0 or 1 and repeated across
#'   the rows of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("binomial_data", package = "EpiExposure")
#' head(binomial_data)
"binomial_data"


#' Simulated gamma epidemiological dataset
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles and a positive continuous plant disease outcome. The dataset is
#' intended for demonstrating gamma models in EpiExposure.
#'
#' The response is repeated across the temporal rows belonging to the same
#' epidemic because each complete exposure profile is associated with one final
#' epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique numeric identifier for each epidemic.}
#'   \item{block}{Experimental block identifier.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{y}{Positive continuous disease outcome repeated across the rows
#'   of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("gamma_data", package = "EpiExposure")
#' head(gamma_data)
"gamma_data"


#' Simulated gaussian epidemiological dataset
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles and a continuous plant disease outcome. The dataset is intended
#' for demonstrating gaussian models in EpiExposure.
#'
#' The response is repeated across the temporal rows belonging to the same
#' epidemic because each complete exposure profile is associated with one final
#' epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique numeric identifier for each epidemic.}
#'   \item{block}{Experimental block identifier.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{y}{Continuous disease outcome repeated across the rows of each
#'   epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("gaussian_data", package = "EpiExposure")
#' head(gaussian_data)
"gaussian_data"


#' Simulated negative-binomial epidemiological dataset
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles and an overdispersed count disease outcome. The dataset is
#' intended for demonstrating quadratic-variance negative-binomial models
#' in EpiExposure.
#'
#' The response is repeated across the temporal rows belonging to the same
#' epidemic because each complete exposure profile is associated with one final
#' epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique numeric identifier for each epidemic.}
#'   \item{block}{Experimental block identifier.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{y}{Non-negative count disease outcome repeated across the rows
#'   of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("nb2_data", package = "EpiExposure")
#' head(nb2_data)
"nb2_data"


#' Simulated poisson epidemiological dataset
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles and a count plant disease outcome. The dataset is intended for
#' demonstrating poisson models in EpiExposure.
#'
#' The response is repeated across the temporal rows belonging to the same
#' epidemic because each complete exposure profile is associated with one final
#' epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique numeric identifier for each epidemic.}
#'   \item{block}{Experimental block identifier.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{y}{Non-negative count disease outcome repeated across the rows
#'   of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("poisson_data", package = "EpiExposure")
#' head(poisson_data)
"poisson_data"


#' Simulated epidemiological data with independent spatial replicates
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles for fixed field locations observed across independent
#' epidemiological replicates. The dataset is intended for demonstrating
#' distributed exposure-lag models with separate spatial fields across
#' epidemic replicates.
#'
#' The response and location-level descriptors are repeated across the temporal
#' rows belonging to the same epidemic because each complete exposure profile
#' is associated with one final epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique identifier combining field location and
#'   epidemiological replicate.}
#'   \item{location_id}{Identifier for the field location.}
#'   \item{epidemic_replicate}{Identifier for an independent epidemiological
#'   replicate.}
#'   \item{block}{Experimental block identifier within a replicate.}
#'   \item{block_id}{Identifier combining epidemiological replicate and
#'   experimental block.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{x_coord}{First numeric spatial coordinate of the field location.}
#'   \item{y_coord}{Second numeric spatial coordinate of the field location.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{y}{Non-negative count disease outcome repeated across the rows
#'   of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("spatial_replicate", package = "EpiExposure")
#' head(spatial_replicate)
"spatial_replicate"


#' Simulated epidemiological data with independent spatial fields by year
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles for fixed field locations observed across multiple years. The
#' dataset is intended for demonstrating distributed exposure-lag models with
#' separate spatial fields by year.
#'
#' The response and location-level descriptors are repeated across the temporal
#' rows belonging to the same epidemic because each complete exposure profile
#' is associated with one final epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique identifier combining field location and year.}
#'   \item{location_id}{Identifier for the field location.}
#'   \item{year}{Year identifying the independent spatial field.}
#'   \item{block}{Experimental block identifier within a year.}
#'   \item{block_id}{Identifier combining year and experimental block.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{x_coord}{First numeric spatial coordinate of the field location.}
#'   \item{y_coord}{Second numeric spatial coordinate of the field location.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{y}{Non-negative count disease outcome repeated across the rows
#'   of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("spatial_year", package = "EpiExposure")
#' head(spatial_year)
"spatial_year"


#' Simulated spatial poisson epidemiological dataset
#'
#' A simulated longitudinal dataset containing complete environmental exposure
#' profiles and spatial coordinates for demonstrating distributed
#' exposure-lag models with Matérn spatial dependence and a poisson disease
#' outcome.
#'
#' The response, block, and spatial coordinates are repeated across the
#' temporal rows belonging to the same epidemic because each complete exposure
#' profile is associated with one final epidemic-level disease outcome.
#'
#' @format A data frame with rows representing temporal observations within
#'   epidemic exposure profiles and the following variables:
#' \describe{
#'   \item{epi_id}{Unique numeric identifier for each epidemic.}
#'   \item{block}{Experimental block identifier.}
#'   \item{time}{Chronological time index within each epidemic profile.}
#'   \item{x_coord}{First numeric spatial coordinate of the epidemic location.}
#'   \item{y_coord}{Second numeric spatial coordinate of the epidemic location.}
#'   \item{tmean}{Daily mean temperature in degrees Celsius.}
#'   \item{wetness}{Daily leaf wetness duration in hours.}
#'   \item{rain}{Daily rainfall in millimeters.}
#'   \item{y}{Non-negative count disease outcome repeated across the rows
#'   of each epidemic profile.}
#' }
#'
#' @source Simulated for the EpiExposure package.
#'
#' @examples
#' data("st_poisson", package = "EpiExposure")
#' head(st_poisson)
"st_poisson"
