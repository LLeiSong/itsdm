#' A function to parse future BIOs from CMCC-BioClimInd.
#' @description This function allows you to parse future CMCC-BioClimInd
#' bioclimatic variables with a setting of boundary and a few other options.
#' @param bry (sf or sp) The boundary to mask the data.
#' if NULL, no clip to the global map. The default is NULL.
#' @param path (character) The path to save the downloaded imagery.
#' If NULL, then use the current working directory. The default is NULL.
#' @param esm (character) The option for  Earth System Models (ESMs).
#' Should be one of "CMCC-CESM", 'GFDL-ESM2M', 'HadGEM2-ES',
#' 'IPSL-CM5A-LR', 'MIROC-ESM-CHEM', 'NorESM1-M'.
#' The default is CMCC-CESM.
#' @param rcp (numeric) The option of Representative Concentration Pathways (RCPs).
#' Should be 45 or 85. Only 85 is available for CMCC-CESM.
#' The default is 85.
#' @param interval (character) The option for time interval.
#' Should be one of "2040-2079", "2060-2099".
#' The default is 2040-2079.
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
#' future_cmcc_bioclim()
#' future_cmcc_bioclim(return_stack = FALSE)
future_cmcc_bioclim <- function(bry = NULL,
                                path = NULL,
                                esm = 'CMCC-CESM',
                                rcp = 85,
                                interval = "2040-2079",
                                nm_mark = "clip",
                                return_stack = TRUE) {
    # Check the inputs
    stopifnot(esm %in% c("CMCC-CESM", 'GFDL-ESM2M', 'HadGEM2-ES',
                         'IPSL-CM5A-LR', 'MIROC-ESM-CHEM', 'NorESM1-M'))
    stopifnot(rcp %in% c(45, 85))
    if(esm == "CMCC-CESM" & rcp == 45){
        stop("There is no such var to download.")}
    stopifnot(interval %in% c("2040-2079", "2060-2099"))
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
    ## Convert inputs
    esm <- switch(esm,
                  `CMCC-CESM` = "CMCC",
                  `GFDL-ESM2M` = "GFDL",
                  `HadGEM2-ES` = "HADGEM",
                  `IPSL-CM5A-LR` = 'IPSL',
                  `MIROC-ESM-CHEM` = 'MIROC',
                  `NorESM1-M` = 'NORESM')
    interval <- switch(interval,
                       `2040-2079` = '2040_79',
                       `2060-2099` = '2060_99')

    ## Path
    path <- file.path(path, "cmcc_bioclim")
    dir.create(path, showWarnings = FALSE)

    # Download and extract historical variables
    url_base <- "https://hs.pangaea.de/model/NoceS-etal_2019"
    invisible(lapply(paste0('BIO', 1:35), function(var){
        if (!file.exists(
            sprintf("%s/%s_%s_%s_%s.nc",
                    path, var, esm, rcp, interval))){
            zip_name <- sprintf("%s.zip", var)
            url <- file.path(url_base, zip_name)

            # Download to local
            temp <- tempfile()
            options(timeout = 1e5)
            download.file(url, temp)

            # Extract hist file from downloaded zip
            decompression <- system2(
                "unzip",
                args = c("-j", "-o", temp,
                         sprintf("%s_%s_%s_%s.nc",
                                 var, esm, rcp, interval),
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
        path, sprintf("BIO%s_%s_%s_%s.nc", 1:35, esm, rcp, interval))
    if (length(intersect(imgs, imgs_in)) != 35) {
        stop("Wrong file numbers unzipped.")}

    # Read files
    clip_imgs <- stack(imgs) %>% st_as_stars()
    clip_imgs <- st_set_dimensions(
        clip_imgs, 'band', values = paste0('bio', 1:35))
    names(clip_imgs) <- sprintf('bioclim_%s_%s_%s', esm, rcp, interval)

    if (!is.null(bry)) {
        # Read files
        bry <- st_as_sf(bry) %>% st_make_valid()
        clip_imgs <- st_crop(clip_imgs, bry)
    }

    ## Save out
    rst_name <- paste(nm_mark,
                      sprintf('bioclim_%s_%s_%s.tif', esm, rcp, interval),
                      sep = '_')
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

# future_cmcc_bioclim end