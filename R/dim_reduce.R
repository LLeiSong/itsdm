#' @title Remove environmental variables that have high correlation with others.
#' @description Select environmental variables that have pairwise Pearson
#' correlation lower than a user-defined threshold.
#' NOTE that it only works on numeric variables, does not work on
#' categorical variables.
#' @param img_stack (\code{stars} or \code{RasterStack})
#' The image stack to work on.
#' @param threshold (\code{numeric}) The threshold number of Pearson
#' correlation that indicates two variables are strongly correlated.
#' The default is 0.5.
#' @param preferred_vars (\code{vector} of \code{character}) The preferred
#' variables \bold{in order} in dimension reduction. The preferred variables
#' will move to the beginning before the reduction. So make sure they are
#' placed in order. Furthermore, setting preferred_vars does not guarantee
#' they can survive. For example, one preferred variable that is placed later
#' has strong correlation with former preferred variable.
#' @param samples (\code{\link{sf}} or \code{sp}) The samples to reduce
#' dimension.
#' If not \code{NULL}, it can take \code{\link{sf}}, \code{\link{sfc}},
#' \code{SpatialPointsDataFrame}, \code{SpatialPoints}, etc.
#' If \code{NULL}, the whole raster stack would be used.
#' The default is \code{NULL}.
#' @return (\code{ReducedImageStack}) A list of
#' \itemize{
#' \item{threshold (\code{numeric}) The threshold set in function inputs}
#' \item{img_reduced (\code{stars}) The image stack after dimension reduction}
#' \item{cors_original (\code{\link{data.frame}}) A table of Pearson
#' correlations between all variables.}
#' \item{cors_reduced (\code{\link{data.frame}}) A table of Pearson
#' correlations between variables after dimension reduction.}}
#' @import checkmate
#' @importFrom sf st_as_sf st_crop
#' @importFrom raster stack layerStats mask rasterize subset
#' @importFrom stars st_as_stars
#' @importFrom dplyr between select
#' @importFrom rlang is_empty
#' @importFrom methods is as
#' @export
#' @examples
#' library(sf)
#' library(itsdm)
#' library(stars)
#' library(dplyr)
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars()
#' img_reduced <- dim_reduce(env_vars, threshold = 0.7,
#'   preferred_vars = c('bio1', 'bio12'))
#'
dim_reduce <- function(img_stack = NULL,
                       threshold = 0.5,
                       preferred_vars = NULL,
                       samples = NULL) {
    # Check inputs
    stopifnot(is.numeric(threshold) & between(threshold, 0, 1))
    stopifnot(is(img_stack, 'stars') | is(img_stack, 'RasterStack'))
    if (is.null(samples)) {
        message("No samples set, use whole image.")
    } else{
        if (!(is(samples, "sf") | is(samples, 'sfc') |
              is(samples, 'SpatialPoints') |
              is(samples, "SpatialPointsDataFrame"))) {
            stop("Only support sf or sp.")
        }
    }
    checkmate::assert_vector(preferred_vars, null.ok = T)

    # Convert to raster to calculate correlations
    if_stars <- is(img_stack, 'stars')
    if (if_stars) {
        if (length(dim(img_stack)) == 2) {
            img <- stack(as(img_stack, 'Raster'))
        } else {
          img_stack <- split(img_stack)
          img <- stack(as(img_stack, 'Raster'))
        }
    } else img <- img_stack

    # Check preferred variables are all in image stack
    if (is.null(preferred_vars)) preferred_vars <- names(img)
    if (!all(preferred_vars %in% names(img))) {
        stop('Some of preferred_vars are not in image stack.')
    }

    # Extract samples if set any
    if (!is.null(samples)){
        samples <- st_as_sf(samples)
        if (st_crs(samples) != st_crs(img)){
          samples <- st_transform(samples, st_crs(img))
        }
        samples <- rasterize(samples, img[[1]], 1)
        img <- mask(img, samples)
    }

    # Calculate correlations
    stat <- "pearson" # Just use pearson because it is standardized.
    cors <- layerStats(img, stat, na.rm = T)
    ps_cor <- data.frame(cors[[1]])
    ids <- match(preferred_vars, names(ps_cor))
    ids <- c(ids, setdiff(1:nrow(ps_cor), ids))
    ps_cor <- ps_cor[ids, ids]
    i <- 1
    while (TRUE) {
        if(i > ncol(ps_cor)) break
        row_index <- which(abs(ps_cor[, i]) > threshold &
                               abs(ps_cor[, i]) < 0.9999)
        if(!is_empty(row_index)) ps_cor <- ps_cor[-row_index, -row_index]
        i <- i + 1
    }

    # Subset images and make object
    if (if_stars) {
        img_reduced <- img_stack %>% select(row.names(ps_cor)) %>%
          merge(name = 'band')
        names(img_reduced) <- 'reduced_image'
    } else {
      img_reduced <- raster::subset(img_stack, row.names(ps_cor))
        }
    img_reduced <- list(threshold = threshold,
                        img_reduced = img_reduced,
                        cors_original = cors,
                        cors_reduced = ps_cor)
    class(img_reduced) <- 'ReducedImageStack'

    # Print and return
    img_reduced
}

# dim_reduce end
