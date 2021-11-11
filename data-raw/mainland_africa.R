## code to prepare `mainland_africa` dataset goes here
library(rnaturalearth, quietly = T)

# Get Africa continent
mainland_africa <- ne_countries(
  continent = 'africa', returnclass = 'sf') %>%
  filter(admin != 'Madagascar') # remove Madagascar

# Union countries to continent
mainland_africa <- mainland_africa %>%
  st_buffer(0.1) %>%
  st_union() %>%
  st_as_sf() %>%
  rename(geometry = x) %>%
  st_make_valid()

# Save out
usethis::use_data(mainland_africa, overwrite = TRUE)
