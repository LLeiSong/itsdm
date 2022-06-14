#' A function to parse the future climate from worldclim version 2.1.
#' @description This function allows you to parse worldclim version 2.1
#' future climatic files with a setting of boundary and a few other options.
#' @param var (\code{character}) The option for the variable to download.
#' Should be one of tmin, tmax, prec, bioc.
#' The default is tmin.
#' @param res (\code{numeric}) The option for the resolution of image to
#' download. Should be one of 2.5, 5, 10. The default is 10.
#' @param gcm (\code{character}) The option for global climate models.
#' Should be one of "BCC-CSM2-MR", "CNRM-CM6-1","CNRM-ESM2-1", "CanESM5",
#' "GFDL-ESM4", "IPSL-CM6A-LR","MIROC-ES2L", "MIROC6", "MRI-ESM2-0".
#' The default is 'BCC-CSM2-MR'.
#' @param ssp (\code{character}) The option for Shared Socio-economic Pathways.
#' Should be one of "ssp126", "ssp245", "ssp370", "ssp585".
#' The default is "ssp585".
#' @param interval (\code{character}) The option for time interval.
#' Should be one of "2021-2040", "2041-2060", "2061-2080", "2081-2100".
#' The default is "2021-2040".
#' @param bry (\code{\link{sf}} or \code{sp}) The boundary to mask the
#' downloaded original data. If \code{NULL}, it would get global map.
#' If not \code{NULL}, it can take \code{\link{sf}},  \code{\link{sfc}},
#' \code{SpatialPolygonsDataFrame}, \code{SpatialPolygons}, etc.
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
#' Fick, Stephen E., and Robert J.
#' Hijmans. "WorldClim 2: new 1-km spatial resolution climate surfaces for
#' global land areas." \emph{International journal of climatology}
#' 37.12 (2017): 4302-4315.\doi{10.1002/joc.5086}
#'
#' @details
#' \href{https://worldclim.org/data/index.html}{Web page page for this dataset}
#'
#' @note
#' If it fails due to slow internet, try to set a larger timeout option,
#' e.g., using `options(timeout = 1e3)`.
#'
#' @importFrom sf st_as_sf st_make_valid st_crop
#' @importFrom stars read_stars write_stars st_set_dimensions
#' @importFrom utils tail download.file
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' future_worldclim2("tmin", 10, "BCC-CSM2-MR",
#'   "ssp585", "2021-2040",
#'   path = tempdir(), return_stack = FALSE)
#'}
#'
future_worldclim2 <- function(var = "tmin",
                              res = 10,
                              gcm = "BCC-CSM2-MR",
                              ssp = "ssp585",
                              interval = "2021-2040",
                              bry = NULL,
                              path = NULL,
                              nm_mark = "clip",
                              return_stack = TRUE) {
    ## Check the inputs
    stopifnot(res %in% c(2.5, 5, 10))
    stopifnot(var %in% c('tmin', 'tmax', 'prec', 'bioc'))
    stopifnot(gcm %in% c("BCC-CSM2-MR", "CNRM-CM6-1",
                         "CNRM-ESM2-1", "CanESM5",
                         "GFDL-ESM4", "IPSL-CM6A-LR",
                         "MIROC-ES2L", "MIROC6",
                         "MRI-ESM2-0"))
    stopifnot(ssp %in% c("ssp126", "ssp245",
                         "ssp370", "ssp585"))
    stopifnot(interval %in% c("2021-2040", "2041-2060",
                              "2061-2080", "2081-2100"))
    if(gcm == "GFDL-ESM4" & ssp == "ssp245"){
        stop("There is no any var to download.")}
    if (gcm == "GFDL-ESM4" & ssp == "ssp585" &
        var != "prec"){
        stop("There is no such var to download.")}

    if (is.null(bry)) {
      nm_mark <- "global"
        message("No bry set, download global map...")
    } else{
        if (!(is(bry, "sf") | is(bry, 'sfc') |
              is(bry, "SpatialPolygonsDataFrame") |
              is(bry, 'SpatialPolygons'))) {
            stop("Only support sf or sp.")
        }
    }
    if (is.null(path)) {
        path <- getwd()
    } else {
        if (!dir.exists(path)) {
            stop("Path does not exist!")
        }
    }

    ## Set up
    path <- file.path(path, "wc2.1")
    dir.create(path, showWarnings = FALSE)

    ## Prepare url and file name
    url_base <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/fut"
    zip_name <- sprintf("%sm/wc2.1_%sm_%s_%s_%s_%s.zip",
                        res, res, var, gcm, ssp, interval)
    url <- file.path(url_base, zip_name)

    ## Download to local
    temp <- tempfile()
    dl <- try(download.file(url, temp))
    if (inherits(dl, "try-error")) {
      Sys.sleep(10)
      download.file(url, temp)
    }

    # Define file number
    n <- ifelse(var == "bioc", 19, 12)

    if (is.null(bry)) {
        ## More stable way to unzip a huge file
        decompression <- system2(
            "unzip",
            args = c("-o", "-j", temp, sprintf("-d %s", path)),
            stdout = TRUE)
        if (grepl("Warning message", tail(decompression, 1))) {
            print(decompression)
        }; unlink(temp)

        fname <- sprintf("wc2.1_%sm_%s_%s_%s_%s.tif",
                         res, var, gcm, ssp, interval)
        fpath <- file.path(path, fname)
        if (!file.exists(fpath)) stop('Fail to extract file.')

        ## Read imgs as stars
        if (return_stack == TRUE) clip_imgs <- read_stars(fpath)
    } else {
        temp_path <- file.path(path, "global")
        dir.create(temp_path, showWarnings = FALSE)

        ## More stable way to unzip a huge file
        decompression <- system2(
            "unzip",
            args = c("-o", "-j", temp, sprintf("-d %s", temp_path)),
            stdout = TRUE)
        if (grepl("Warning message", tail(decompression, 1))) {
            print(decompression)
        }; unlink(temp)

        fname <- sprintf("wc2.1_%sm_%s_%s_%s_%s.tif",
                         res, var, gcm, ssp, interval)
        fpath <- file.path(temp_path, fname)
        if (!file.exists(fpath)) stop('Fail to extract file.')

        ## Read imgs as stars and clip
        clip_imgs <- read_stars(fpath, RasterIO = list(bands = c(1:n)))
        bry <- st_as_sf(bry) %>% st_make_valid()
        clip_imgs <- st_crop(clip_imgs, bry)

        ## Save out
        ### No 30s resolution, so stacking them together is fine.
        rst_name <- paste(nm_mark, names(clip_imgs), sep = "_")
        rst_path <- file.path(path, rst_name)
        write_stars(clip_imgs, rst_path)

        ## Clean the temp folder
        unlink(temp_path, recursive = TRUE)

        if (return_stack != TRUE) rm(clip_imgs)
    }

    if (return_stack == TRUE) {
        clip_imgs <- st_set_dimensions(
            clip_imgs, 'band', values = paste0(var, 1:n))
        clip_imgs
    } else {
      message(sprintf("Files are written to %s.", path))
    }
}

# future_worldclim2 end
