## code to prepare `occ_virtual_species` dataset goes here
library(itsdm, quietly = T)
library(virtualspecies, quietly = T)

bios <- system.file('extdata/bioclim_africa_10min.tif', package = 'itsdm') %>%
  read_stars() %>% slice('band', c(1, 12))

# Convert to raster for virtualspecies
bios <- stack(as(split(bios), 'Spatial'))

# Formatting of the response functions
set.seed(10)
my.parameters <- formatFunctions(
  bio1 = c(fun = 'dnorm', mean = 25, sd = 10),
  bio12 = c(fun = 'dnorm', mean = 1000, sd = 500))

# Generation of the virtual species
set.seed(10)
my.species <- generateSpFromFun(
  raster.stack = bios,
  parameters = my.parameters,
  plot = F)

# Conversion to presence-absence
set.seed(10)
my.species <- convertToPA(
  my.species,
  beta = 0.7,
  plot = F)

# Pseudo presence-only observations
set.seed(10)
po.points <- sampleOccurrences(
  my.species,
  n = 2000,
  type = "presence only")
occ_virtual_species <- po.points$sample.points %>%
  select(x, y)

# Save out
usethis::use_data(occ_virtual_species, overwrite = TRUE)
