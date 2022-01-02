## code to prepare `wdpa` dataset goes here
library(itsdm, quietly = T)
library(dplyr, quietly = T)
library(stars)

# Get a raster template
data("mainland_africa")
bios <- worldclim2(var = 'bio',
                   bry = mainland_africa,
                   path = tempdir(),
                   nm_mark = 'africa') %>%
  st_normalize()
stars_template <- bios %>% slice('band', 1) %>%
  mutate('wc2.1_10m_bio.tif' = 0)

# Get protected area
# Manually download from website. Note that the link will change everyday.
# https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_Nov2021_Public_AF_shp.zip
# UNEP-WCMC and IUCN (2021), Protected Planet: The World Database on Protected
# Areas (WDPA) and World Database on Other Effective Area-based Conservation
# Measures (WD-OECM) [Online], November 2021, Cambridge, UK: UNEP-WCMC and IUCN.
# Available at: www.protectedplanet.net.
# Download
urlpath <- file.path('https://d1gam3xoknrgr2.cloudfront.net/current',
                     'WDPA_WDOECM_Nov2021_Public_AF_shp.zip')
temp <- tempfile(); temp_dir <- tempdir()
options(timeout = 1e5)
download.file(urlpath, temp)
unzip(temp, exdir = temp_dir)

# Read polygons only
base_name <- 'WDPA_WDOECM_Nov2021_Public_AF_shp'
invisible(lapply(0:2, function(n) {
  unzip(file.path(temp_dir, sprintf('%s_%s.zip', base_name, n)),
        exdir = file.path(temp_dir, sprintf('%s_%s', base_name, n)))}))
wdpa <- do.call(rbind, lapply(0:2, function(n) {
  read_sf(file.path(temp_dir, sprintf('%s_%s', base_name, n),
                    sprintf('%s-polygons.shp', base_name))) %>%
  st_make_valid()
})); unlink(temp)

# Rasterize polygons to match with bios
wdpa <- st_rasterize(wdpa %>% mutate(wdpa = 1) %>% select(wdpa),
                     template = stars_template) %>%
  st_crop(., mainland_africa) %>%
  mutate(wdpa = factor(wdpa))
write_stars(wdpa, NA_value = 255, type = 'Byte',
            'inst/extdata/wdpa_africa_10min.tif')

# Do not let it lazy loading, keep it as an external GeoTiff file
# usethis::use_data(wdpa, overwrite = TRUE)
