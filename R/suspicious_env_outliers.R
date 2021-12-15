#' @title Function to detect suspicious outliers based on environmental variables.
#' @description Run \code{\link{outlier.tree}} to detect suspicious outliers in observations.
#' @param occ (`data.frame`, `sf`, `SpatialPointsDataFrame`)
#' The occurrence dataset for training.
#' There must be column `x` and `y` for coordinates if it is a regular `data.frame`.
#' @param occ_crs (`numeric` or \code{\link{crs}}) The EPSG number or
#' \code{\link{crs}} object of occurrence CRS.
#' The default value is `4326`, which is the geographic coordinate system.
#' @param variables (`RasterStack` or `stars`) The stack of environmental variables.
#' @param rm_outliers (`logical`) The option to remove the suspicious outliers or not.
#' The default is `FALSE`.
#' @param seed (`integer`) The random seed used in the modeling. It should be an
#' integer. The default is `10L`.
#' @param ... Other arguments passed to function \code{\link{outlier.tree}} in
#' package `outliertree`.
#' @return (`EnvironmentalOutlier`) A list that contains
#' \itemize{
#' \item{outliers (\code{\link{sf}}) The \code{\link{sf}} points of outliers}
#' \item{outlier_details (`tibble`) A table of outlier details returned from
#' function \code{\link{outlier.tree}} in package `outliertree`}
#' \item{pts_occ (\code{\link{sf}}) The \code{\link{sf}} points of occurrence.
#' If `rm_outliers` is `TRUE`, outliers are deleted from points of
#' occurrence. If `FALSE`, the full observations are returned.}}
#'
#' @seealso
#' \code{\link{print.EnvironmentalOutlier}}, \code{\link{plot.EnvironmentalOutlier}}
#' \code{\link{outlier.tree}} in package `outliertree`
#'
#' @details
#' Please check more details in R documentation of function
#' \code{\link{outlier.tree}} in package `outliertree` and their GitHub.
#'
#' @references
#' \itemize{
#' \item{\href{https://arxiv.org/abs/2001.00636}{Cortes, David. "Explainable
#' outlier detection through decision tree conditioning."
#' \emph{arXiv preprint arXiv:2001.00636} (2020).}}
#' \item{\url{https://github.com/david-cortes/outliertree}}}
#'
#' @import checkmate
#' @importFrom methods is
#' @importFrom stars st_as_stars st_extract
#' @importFrom sf st_as_sf st_crs st_transform st_drop_geometry
#' @importFrom outliertree outlier.tree
#' @export
#' @examples
#' library(dplyr)
#' library(sf)
#' library(stars)
#' library(itsdm)
#'
#' data("occ_virtual_species")
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12, 16))
#'
#' occ_outliers <- suspicious_env_outliers(
#'   occ = occ_virtual_species, variables = env_vars,
#'   z_outlier = 3.5, outliers_print = 4L)
#'
#' occ_outliers
#' plot(occ_outliers)
#'
suspicious_env_outliers <- function(occ,
                                    occ_crs = 4326,
                                    variables,
                                    rm_outliers = FALSE,
                                    seed = 10L,
                                    ...) {
  # Check inputs
  checkmate::assert_multi_class(
    occ, c('data.frame', 'sf',
           'SpatialPointsDataFrame',
           null.ok = F))
  if (is(occ, 'data.frame') & (!is(occ, 'sf')) &
      (!is(occ, 'Spatial'))) {
    if(!all(c("x", "y") %in% colnames(occ))){
      stop("There must be x and y column in occ.")}}
  checkmate::assert_multi_class(
    variables, c('RasterStack', 'RasterLayer', 'stars'))

  # Reformat the inputs
  # Variables, use stars
  if (is(variables, 'RasterStack') | is(variables, 'RasterLayer')){
    variables <- st_as_stars(variables)
  }
  # Occurrence
  if (is(occ, 'data.frame') & (!is(occ, 'sf')) &
      (!is(occ, 'Spatial'))) {
    pts_occ <- occ %>%
      st_as_sf(coords = c('x', 'y'), crs = occ_crs)
  } else {
    pts_occ <- st_as_sf(occ)
  }
  if (st_crs(variables) != st_crs(pts_occ)){
    pts_occ <- st_transform(pts_occ, st_crs(variables))
  }

  # Prepare variable values
  var_occ <- st_extract(x = variables, at = pts_occ) %>%
    st_as_sf() %>% st_drop_geometry()

  # Detect suspicious environmental outliers
  set.seed(seed)
  outliers_model <- outlier.tree(
    var_occ,
    ...,
    save_outliers = T)

  outliers <- do.call(
    rbind, lapply(1:length(outliers_model$outliers_data), function(n) {
      lst <- outliers_model$outliers_data[[n]]
      if (!is.na(lst$outlier_score)){
        tibble(suspious_row = n,
               suspicous_column = lst$suspicous_value$column,
               suspicous_value = lst$suspicous_value$value,
               outlier_score = lst$outlier_score)}}))
  outliers_details <- outliers_model$outliers_data[unlist(outliers[, 'suspious_row'])]

  # Remove outliers
  if (!is.null(outliers)) {
    if (rm_outliers) {
      pts_occ <- pts_occ %>%
        slice(-unlist(outliers[, 'suspious_row']))
    }
    outliers <- pts_occ %>%
      slice(unlist(outliers[, 'suspious_row'])) %>%
      cbind(outliers)
  }

  # Return
  out <- list(outliers = outliers,
              outlier_details = outliers_details,
              pts_occ = pts_occ)
  class(out) <- append("EnvironmentalOutlier", class(out))
  out
}
