#' @title Download environmental variables made by worldclim version 2.1.
#' @description Parse historic worldclim version 2.1 variables
#' with a setting of boundary and a few other options.
#' @param var (\code{character}) The option for the variable to download,
#' should be one of tvag, tmin, tmax, prec, srad, wind, vapr and bio.
#' The default is 'tmin'.
#' @param res (\code{numeric}) The option for the resolution of image to download.
#' Should be one of 0.5, 2.5, 5, 10 in minute degree.
#' The default is 10.
#' @param bry (\code{\link{sf}} or \code{sp}) The boundary to mask the downloaded original data.
#' If \code{NULL}, it would get global map. If not \code{NULL}, it can take \code{\link{sf}},
#' \code{\link{sfc}}, \code{SpatialPolygonsDataFrame}, \code{SpatialPolygons}, etc.
#' The default is \code{NULL}.
#' @param path (\code{character}) The path to save the downloaded imagery.
#' If \code{NULL}, it would use the current working directory.
#' The default is \code{NULL}.
#' @param nm_mark (\code{character}) the name mark of clipped images.
#' The default is "clip". It would be ignored if \code{bry} is \code{NULL}.
#' @param return_stack (\code{logical}) if \code{TRUE}, stack the imagery together and return.
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
#' @importFrom methods is
#' @importFrom utils download.file tail
#' @importFrom sf st_as_sf st_make_valid st_is_valid st_crop
#' @importFrom stars read_stars write_stars
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' library(itsdm)
#'
#' bry <- sf::st_polygon(
#'   list(rbind(c(29.34, -11.72), c(29.34, -0.95),
#'              c(40.31, -0.95), c(40.31, -11.72),
#'              c(29.34, -11.72)))) %>%
#'   st_sfc(crs = 4326)
#'
#' bios <- worldclim2(var = "tmin", res = 10,
#'   bry = bry, nm_mark = 'exp', path = tempdir())
#'}
#'
worldclim2 <- function(var = "tmin",
                       res = 10,
                       bry = NULL,
                       path = NULL,
                       nm_mark = "clip",
                       return_stack = TRUE) {
    # Check the inputs
    ## Vars
    stopifnot(var %in% c(
        "tavg", "tmin", "tmax", "prec",
        "srad", "wind", "vapr", "bio"
    ))
    ## res
    stopifnot(res %in% c(0.5, 2.5, 5, 10))
    ## bry
    if (is.null(bry)) {
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
    if (res == 0.5) {
        res <- "30s"
    } else {
        res <- sprintf("%sm", res)
    }
    path <- file.path(path, "wc2.1")
    dir.create(path, showWarnings = FALSE)

    # Prepare url and file name
    url_base <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base"
    zip_name <- sprintf("wc2.1_%s_%s.zip", res, var)
    url <- file.path(url_base, zip_name)

    # Download to local
    temp <- tempfile()
    dl <- try(download.file(url, temp))
    if (inherits(dl, "try-error")) {
      Sys.sleep(10)
      download.file(url, temp)
    }

    # Define file number
    n <- ifelse(var == "bio", 19, 12)

    # Process
    if (is.null(bry)) {
        # More stable way to unzip a huge file
        decompression <- system2("unzip",
                                 args = c("-j", "-o", temp, sprintf("-d %s", path)),
                                 stdout = TRUE)
        if (grepl("Warning message", tail(decompression, 1))) {
            print(decompression)
        }
        unlink(temp)

        if (return_stack == TRUE) {
            ## Check unzip files
            imgs_in <- list.files(
                path,
                pattern = "*.tif", full.names = T
            )
            if(var == 'bio'){
                imgs <- file.path(
                    path, sprintf("wc2.1_%s_%s_%s.tif", res, var, 1:n)
                )
            } else {
                imgs <- file.path(
                    path, sprintf("wc2.1_%s_%s_%02d.tif", res, var, 1:n)
                )
            }
            if (length(intersect(imgs, imgs_in)) != n) {
                stop("Wrong file numbers unzipped.")
            }

            ## Read imgs as stars
            if (return_stack == TRUE) clip_imgs <- read_stars(imgs)
        }
    } else {
        # Unzip to a temporary path
        ## Define temp dir
        temp_path <- file.path(path, "global")
        dir.create(temp_path, showWarnings = FALSE)

        ## Unzip
        decompression <- system2(
            "unzip",
            args = c("-j", "-o", temp, sprintf("-d %s", temp_path)),
            stdout = TRUE)
        if (grepl("Warning message", tail(decompression, 1))) {
            print(decompression)
        }
        unlink(temp)

        # Clip the imagery to the boundary
        ## Check unzip files
        imgs_in <- list.files(temp_path, pattern = "*.tif", full.names = T)
        if(var == 'bio'){
            imgs <- file.path(
                temp_path, sprintf("wc2.1_%s_%s_%s.tif", res, var, 1:n)
            )
        } else {
            imgs <- file.path(
                temp_path, sprintf("wc2.1_%s_%s_%02d.tif", res, var, 1:n)
            )
        }
        if (length(intersect(imgs, imgs_in)) != n) {
            stop("Wrong file numbers unzipped.")}

        ## Read imgs as stars and clip
        clip_imgs <- read_stars(imgs)
        bry <- st_as_sf(bry)
        if(!st_is_valid(bry)) {
            stop('Not valid boundary.')
        }
        clip_imgs <- st_crop(clip_imgs, bry)
        if (inherits(clip_imgs, "stars_proxy")) {
          clip_imgs <- st_as_stars(clip_imgs)}

        ## Save out
        invisible(lapply(1:n, function(n) {
            rst_name <- paste(nm_mark, names(clip_imgs)[n], sep = "_")
            rst_path <- file.path(path, rst_name)
            write_stars(clip_imgs, rst_path, layer = n)
        }))

        ## Clean the temp folder
        unlink(temp_path, recursive = TRUE)

        if (return_stack != TRUE) rm(clip_imgs)
    }

    if (return_stack == TRUE) {
        names(clip_imgs) <- paste0(var, 1:n)
        clip_imgs <- merge(clip_imgs, name = 'band')
        names(clip_imgs) <- sprintf("wc2.1_%s_%s.tif", res, var)
        clip_imgs
    } else {
        message(sprintf("Files are written to %s.", path))
    }
}

# worldclim2 end
