#' A function to parse historic BIOs from CMCC-BioClimInd.
#' @description This function allows you to parse historic CMCC-BioClimInd
#' bioclimatic variables with a setting of boundary and a few other options.
#' @param bry (sf or sp) The boundary to mask the data.
#' if NULL, no clip to the global map. The default is NULL.
#' @param path (character) The path to save the downloaded imagery.
#' If NULL, then use the current working directory. The default is NULL.
#' @param nm_mark (character) the name mark of clipped images.
#' The default is "clip". It would be ignored if bry is NULL.
#' @param return_stack (logical) if TRUE, stack the imagery together and return.
#' If the area is large and resolution is high, it is better not to stack them.
#' The default is TRUE.
#' #' @return stars if return_stack is TRUE.
#' The images would be saved as a single file.
#' @references Noce, Sergio, Luca Caporaso, and Monia Santini.
#' "A new global dataset of bioclimatic indicators."
#' Scientific data 7.1 (2020): 1-12.
#' \url{https://doi.pangaea.de/10.1594/PANGAEA.904278?format=html#download}
#' @import ncdf4
#' @importFrom glue glue
#' @importFrom raster stack
#' @importFrom sf st_as_sf st_make_valid
#' @importFrom stars read_stars write_stars st_as_stars
#' @export
#' @examples
#' cmcc_bioclim()
#' cmcc_bioclim(return_stack = FALSE)
cmcc_bioclim <- function(bry = NULL,
                         path = NULL,
                         nm_mark = "clip",
                         return_stack = TRUE) {
    # Check the inputs
    ## bry
    if (is.null(bry)) {
        nm_mark <- 'global'
        message("No bry set, download global map.")
    } else {
        if (!(is(bry, "sf") | is(bry, 'sfc') |
              is(bry, "SpatialPolygonsDataFrame") |
              is(bry, 'SpatialPolygons'))) {
            stop("Only support sf or sp.")
        }
    }
    ## path
    if (is.null(path)) {
        path <- getwd()
    } else {
        if (!dir.exists(path)) {
            stop("Path does not exist!")
        }
    }

    # Set up
    path <- file.path(path, "cmcc_bioclim")
    dir.create(path, showWarnings = FALSE)

    # Download and extract historical variables
    url_base <- "https://hs.pangaea.de/model/NoceS-etal_2019"
    invisible(lapply(paste0('BIO', 1:35), function(var){
         if (!file.exists(sprintf("%s/%s_HIST_1960_99.nc", path, var))){
             zip_name <- sprintf("%s.zip", var)
             url <- file.path(url_base, zip_name)

             # Download to local
             temp <- tempfile()
             options(timeout = 1e5)
             download.file(url, temp)

             # Extract hist file from downloaded zip
             decompression <- system2("unzip",
                                      args = c("-j", "-o", temp,
                                               sprintf("%s_HIST_1960_99.nc", var),
                                               sprintf("-d %s", path)),
                                      stdout = TRUE)
             if (grepl("Warning message", tail(decompression, 1))) {
                 print(decompression)
             }
             unlink(temp)
         }
    }))

    ## Check unzipped files
    imgs_in <- list.files(path, pattern = "*.nc", full.names = T)
    imgs <- file.path(
            path, sprintf("BIO%s_HIST_1960_99.nc", 1:35))
    if (length(intersect(imgs, imgs_in)) != 35) {
        stop("Wrong file numbers unzipped.")}

    # Read files
    clip_imgs <- stack(imgs) %>% st_as_stars()
    clip_imgs <- st_set_dimensions(
        clip_imgs, 'band', values = paste0('bio', 1:35))
    names(clip_imgs) <- 'cmcc_bioclim_hist'

    if (!is.null(bry)) {
        # Read files
        bry <- st_as_sf(bry) %>% st_make_valid()
        clip_imgs <- st_crop(clip_imgs, bry)
    }

    ## Save out
    rst_name <- paste(nm_mark, 'cmcc_bioclim_hist.tif', sep = '_')
    rst_path <- file.path(path, rst_name)
    write_stars(clip_imgs, rst_path)

    # Clean temporary files
    unlink(imgs)

    # Return
    if (return_stack == TRUE) {
        clip_imgs
    } else {
        message(glue("Files are written to {path}."))
    }
}

# cmcc_bioclim end
