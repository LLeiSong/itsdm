#' @title Download historic Bioclimatic indicators (BIOs) named CMCC-BioClimInd.
#' @description Parse historic CMCC-BioClimInd bioclimatic indicators
#' optionally with a setting of boundary and a few other options.
#' @param bry (\code{\link[sf:sf]{sf}}) The boundary to mask the
#' downloaded original data. If \code{NULL}, it would get global map.
#' If not \code{NULL}, it can take \code{\link[sf:sf]{sf}},
#' \code{\link[sf:sfc]{sfc}}, etc.
#' The default is \code{NULL}.
#' @param path (\code{character}) The path to save the downloaded imagery.
#' If \code{NULL}, it would use the current working directory.
#' The default is \code{NULL}.
#' @param nm_mark (\code{character}) the name mark of clipped images.
#' The default is "clip". It would be ignored if \code{bry} is \code{NULL}.
#' @param return_stack (\code{logical}) if \code{TRUE}, stack the imagery
#' together and return.
#' If the area is large and resolution is high, it is better not to stack them.
#' The default is \code{TRUE}.
#' @return if \code{return_stack} is \code{TRUE}, the images would be
#' returned as a \code{stars}. Otherwise, nothing to return, but the user
#' would receive a message of where the images are.
#' @references
#' Noce, Sergio, Luca Caporaso, and Monia Santini."A new global dataset of
#' bioclimatic indicators. "\emph{Scientific data} 7.1 (2020): 1-12.
#' \doi{10.1038/s41597-020-00726-5}
#'
#' @details
#' \href{https://doi.pangaea.de/10.1594/PANGAEA.904278?format=html}{Web
#' page page for this dataset}
#'
#' @note The function is experimental at the moment, because the download server
#' of this dataset is not as stable as Worldclim yet. If it fails due to slow
#' internet, try to set a larger timeout option,
#' e.g., using `options(timeout = 1e3)`.
#'
#' @import ncdf4
#' @importFrom raster stack
#' @importFrom sf st_as_sf st_make_valid st_crop
#' @importFrom stars read_stars write_stars st_as_stars st_set_dimensions
#' @importFrom utils download.file tail
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(sf)
#' library(itsdm)
#' bry <- st_polygon(
#'   list(rbind(c(29.34, -11.72), c(29.34, -0.95),
#'              c(40.31, -0.95), c(40.31, -11.72),
#'              c(29.34, -11.72)))) %>%
#'   st_sfc(crs = 4326)
#'
#' cmcc_bios <- cmcc_bioclim(bry = bry,
#'   nm_mark = 'tza', path = tempdir())
#'}
#'
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
        if (!(is(bry, "sf") | is(bry, 'sfc'))) {
            stop("Only support sf.")
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
             dl <- try(download.file(url, temp))
             if (inherits(dl, "try-error")) {
               Sys.sleep(10)
               download.file(url, temp)
             }

             # Extract hist file from downloaded zip
             decompression <- system2(
               "unzip",
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
        message(sprintf("Files are written to %s.", path))
    }
}

# cmcc_bioclim end
