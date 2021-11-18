## code to prepare `landcover` dataset goes here
## Copernicus global land cover map
# Download data from https://lcviewer.vito.be/download
tiles <- c('W020N40', 'E000N40', 'E020N40',
           'W020N20', 'E000N20', 'E020N20', 'E040N20',
                      'E000N00', 'E020N00', 'E040N00',
                      'E000S20', 'E020S20')
base_url <- file.path('https://s3-eu-west-1.amazonaws.com',
                      'vito.landcover.global/v3.0.1/2019')
base_name <- paste0('PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-',
                    'Classification-map_EPSG-4326.tif')
landcover <- do.call(st_mosaic, lapply(tiles, function(tile) {
  temp <- tempfile()
  urlpath <- file.path(base_url, tile,
                       sprintf('%s_%s', tile, base_name))
  download.file(urlpath, temp)
  read_stars(temp)
}))

# Warp to bios and boundary
fname <- 'extdata/bioclim_africa_10min.tif'
stars_template <- system.file(fname, package = 'itsdm') %>%
  read_stars() %>% slice('band', 1)
data("mainland_africa")

landcover <- st_warp(landcover, stars_template,
                     use_gdal = T, method = 'near') %>%
  st_crop(., mainland_africa) %>%
  setNames('landcover') %>%
  mutate(landcover = ifelse(landcover == 0, NA, landcover)) %>%
  mutate(landcover = as.factor(landcover))

write_stars(landcover, NA_value = 255, type = 'Byte',
            'inst/extdata/landcover_africa_10min.tif')

# Do not let it lazy loading, keep it as an external GeoTiff file
# usethis::use_data(landcover, overwrite = TRUE)
