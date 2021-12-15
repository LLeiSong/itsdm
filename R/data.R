#' Boundary of mainland Africa
#'
#' The overall continental boundary of mainland Africa queried from
#' `rnaturalearth` and get processed.
#'
#' @format A \code{\link{sf}} with one rows and 2 fields
#' \describe{
#' \item{name}{(`character`) The name of the polygon: Africa}
#' \item{area}{(\code{\link{units}}) The united number of the overall area in
#' km2. This is not a consensus area, but just a calculated area under this
#' resolution.}
#' \item{geometry}{(\code{\link{sfc}}) The simple polygon feature of the
#' boundary}}
#'
#' @source \code{rnaturalearth}
#'
'mainland_africa'


#' Occurrence dataset of a virtual species
#'
#' A pseudo presence-only occurrence dataset of a virtual species made
#' by package `virtualspecies`.
#'
#' @format A `data.frame` with 300 rows and 2 fields
#' \describe{
#' \item{x}{(`numeric`) The x coordinates of the records in
#' geographic coordinate system}
#' \item{y}{(`numeric`) The y coordinates of the records}}
#'
#' @details
#' The environmental niche of the virtual species is made by defining its
#' response functions to annual temperature and annual precipitation
#' in mainland Africa.
#' The response function of annual temperature is normal distribution with
#' mean = 22 and standard deviation = 5.
#' The response function of annual precipitation is normal distribution with
#' mean = 1000 and standard deviation = 200.
#' Then the suitability is convert to presence-absence map by logistic
#' conversion with beta = 0.7, alpha = -0.05, and species prevalence = 0.27.
#' Finally 300 presence-only points are sampled across the whole region.
#'
#'
#' @source \code{virtualspecies}
#'
'occ_virtual_species'
