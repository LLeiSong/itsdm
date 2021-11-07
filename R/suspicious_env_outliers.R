#' @title Function to detect suspicious outliers based on environmental variables.
#' @description Call outlier.tree to detect suspicious outliers in observations.
#' @param occ (data.frame, sf, SpatialPointsDataFrame)
#' The observation dataset to analyze.
#' There must be column x and y for coordinates if it is a data.frame.
#' @param occ_crs (numeric or crs) The EPSG number or crs object of occurrence CRS.
#' The default value is 4326, which is the geographic coordinate system.
#' @param variables (RasterStack or stars) The stack of environmental variables.
#' @param rm_outliers (logical) The option to remove the suspicious outliers or not.
#' The default is FALSE.
#' @param ... Other arguments passed to function `outlier.tree`.
#' @return (EnvironmentalOutlier) a list that contains `sf` points of observations
#', outliers and outlier details returned from `outlier.tree`.
#'If `rm_outliers` is `TRUE`, outliers are deleted from points of
#' observations. If `FALSE`, the full observations are returned.
#' @import checkmate
#' @importFrom stars st_as_stars st_extract
#' @importFrom sf st_as_sf st_crs st_transform st_drop_geometry
#' @importFrom outliertree outlier.tree
#' @export
#' @examples
#' suspicious_env_outliers(occ = occ, variables = env_vars)
#'
suspicious_env_outliers <- function(occ,
                                    occ_crs = 4326,
                                    variables,
                                    rm_outliers = FALSE,
                                    ...) {
  # Check inputs
  checkmate::assert_multi_class(
    occ, c('data.frame', 'sf',
           'SpatialPointsDataFrame',
           null.ok = F))
  if (is(occ, 'data.frame')) {
    if(!all(c("x", "y") %in% colnames(occ))){
      stop("There must be x and y column in occ.")}}
  checkmate::assert_multi_class(
    variables, c('RasterStack', 'RasterLayer', 'stars'))

  # Reformat the inputs
  # Variables, use stars
  if (is(variables, 'RasterStack') | is(variables, 'RasterLayer')){
    variables <- st_as_stars(variabels)
  }
  # Occurrence
  if (is(occ, 'data.frame')) {
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
  if (rm_outliers) {
    pts_occ <- pts_occ %>%
      slice(-unlist(outliers[, 'suspious_row']))
  }
  outliers <- pts_occ %>%
    slice(unlist(outliers[, 'suspious_row'])) %>%
    cbind(outliers)

  # Return
  out <- list(outliers = outliers,
              outlier_details = outliers_details,
              pts_occ = pts_occ)
  class(out) <- append("EnvironmentalOutlier", class(out))
  out
}
