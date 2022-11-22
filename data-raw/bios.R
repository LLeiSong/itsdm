## code to prepare `bios` dataset goes here
library(itsdm, quietly = T)
library(rnaturalearth, quietly = T)

# Get Africa continent
tza <- ne_countries(
  continent = 'africa', returnclass = 'sf') %>%
  filter(iso_a3 == 'TZA') # Tanzania

# Union countries to continent
tza_bry <- st_buffer(tza, 0.1) %>%
  st_union() %>%
  st_as_sf() %>%
  rename(geometry = x) %>%
  st_make_valid()

bios <- worldclim2(
  var = 'bio', bry = tza_bry,
  path = tempdir(),
  nm_mark = 'tanzania')
write_stars(bios, 'inst/extdata/bioclim_tanzania_10min.tif')

bios_future <- future_worldclim2(
  var = 'bioc', bry = tza_bry,
  path = tempdir(),
  nm_mark = 'tanzania') %>%
  slice("band", c(1, 5, 12))
write_stars(bios_future, 'inst/extdata/future_bioclim_tanzania_10min.tif')

# Do not let it lazy loading, keep it as an external GeoTiff file
# usethis::use_data(bios, overwrite = TRUE)
