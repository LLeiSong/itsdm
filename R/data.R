#' Boundary of mainland Africa
#'
#' The overall continental boundary of mainland Africa queried from `rnaturalearth`
#' and get processed.
#'
#' @format A `sf` with one rows and 2 fields
#' \describe{
#' \item{name (`character`) The name of the polygon: Africa}
#' \item{area (\code{\link{units}}) The united number of the overall area. This
#' is not a consensus area, but just a calculated area under this resolution.}
#' \item{geometry (\code{\link{sfc}}) The simple polygon feature of the boundary}}
#'
#' @source \code{\link{rnaturalearth}}
#'
'mainland_africa'


#' Occurrence dataset of a virtual species
#'
#' A pseudo presence-only occurrence dataset of a virtual species made
#' by package \code{\link{virtualspecies}}.
#'
#' @format A `data.frame` with 2000 rows and 2 fields
#' \describe{
#' \item{x (`numeric`) The x coordinates of the records in geographic coordinate system}
#' \item{y (`numeric`) The y coordinates of the records}}
#'
#' @details
#' The environmental niche of the virtual species is made by defining its
#' response functions to annual temperature and annual precipitation in mainland Africa.
#' The response function of annual temperature is normal distribution with
#' mean = 25 and standard deviation = 10.
#' The response function of annual precipitation is normal distribution with
#' mean = 1000 and standard deviation = 500.
#' Then the suitability is convert to presence-absence map by logistic conversion
#' with beta = 0.7, alpha = -0.05, and species prevalence = 0.27.
#' Finally 2000 presence-only points are sampled across the whole region.
#'
#'
#' @seealso \code{\link{virtualspecies}}
#'
'occ_virtual_species'