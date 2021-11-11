## code to prepare `bios` dataset goes here
library(itsdm, quietly = T)
library(rnaturalearth, quietly = T)

# Get Africa continent
af <- ne_countries(
  continent = 'africa', returnclass = 'sf') %>%
  filter(admin != 'Madagascar') # remove Madagascar

# Union countries to continent
af_bry <- st_buffer(af, 0.1) %>%
  st_union() %>%
  st_as_sf() %>%
  rename(geometry = x) %>%
  st_make_valid()

bios <- worldclim2(var = 'bio', bry = af_bry, nm_mark = 'africa')
write_stars(bios, 'inst/extdata/bioclim_africa_10min.tif')

# Do not let it lazy loading, keep it as an external GeoTiff file
# usethis::use_data(bios, overwrite = TRUE)
