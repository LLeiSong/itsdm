#' A function to remove auto-correlation between raster layers.
#' @description This function allows you to reduce dimensions of raster layers
#' based on pearson correlation and a user-defined threshold.
#' @param img_stack (Stars or RasterStack) The image stack to work on.
#' @param threshold (numeric) The number of the threshold to reduce.
#' The default is 0.5.
#' @param samples (sf or sp) The samples to reduce dimension.
#' If NULL, the whole rasterstack would be used. The default is NULL.
#' @return ReducedImageStack.
#' @importFrom sf st_as_sf
#' @importFrom raster stack layerStats mask rasterize subset
#' @importFrom stars st_as_stars
#' @importFrom dplyr between
#' @importFrom purrr is_empty
#' @export
#' @examples
#' worldclim <- worldclim2(var = "bio")
#' dim_reduce(worldclim)
dim_reduce <- function(img_stack = NULL,
                       threshold = 0.5,
                       samples = NULL) {
    # Check inputs
    stopifnot(is.numeric(threshold) & between(threshold, 0, 1))
    stopifnot(is(img_stack, 'stars') | is(img_stack, 'raster'))
    if (is.null(samples)) {
        message("No samples set, use whole image.")
    } else{
        if (!(is(samples, "sf") | is(samples, 'sfc') |
              is(samples, 'SpatialPoints') |
              is(samples, "SpatialPointsDataFrame"))) {
            stop("Only support sf or sp.")
        }
    }

    # Convert to raster to calculate correlations
    if_stars <- is(img_stack, 'stars')
    if (if_stars) {
        if (length(dim(img_stack)) == 2) {
            img_stack <- stack(as(img_stack, 'Spatial'))
        } else {
            img_stack <- stack(as(split(img_stack), 'Spatial'))
        }
    }

    if (!is.null(samples)){
        samples <- st_as_sf(samples)
        samples <- rasterize(samples, img_stack[[1]], 1)
        img_stack <- mask(img_stack, samples)
    }

    # Calculate correlations
    stat <- "pearson" # Just use pearson because it is standardized.
    cors <- layerStats(img_stack, stat, na.rm = T)
    cors <- data.frame(cors[[1]])
    ps_cor <- cors
    for (i in 1:ncol(ps_cor)){
        if(i > ncol(ps_cor)){
            break}
        row_index <- which(abs(ps_cor[, i]) > threshold &
                               abs(ps_cor[, i]) < 1)
        if(!is_empty(row_index)) ps_cor <- ps_cor[-row_index, -row_index]
    }

    # Subset images and make object
    img_reduced <- raster::subset(img_stack, row.names(ps_cor))
    if (if_stars) img_reduced <- st_as_stars(img_reduced)
    img_reduced <- list(img_reduced = img_reduced,
                        cors_original = cors,
                        cors_reduced = ps_cor)
    class(img_reduced) <- 'ReducedImageStack'
    # Print and return
    cat('Reduced correlations:\n')
    print(img_reduced$cors_reduced)
    invisible(img_reduced)
}

# dim_reduce end
