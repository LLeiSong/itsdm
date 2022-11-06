## code to prepare `occ_virtual_species` dataset goes here
library(itsdm, quietly = TRUE)
library(virtualspecies, quietly = TRUE)
library(dplyr, quietly = TRUE)

bios <- system.file('extdata/bioclim_tanzania_10min.tif', package = 'itsdm') %>%
  read_stars() %>% slice('band', c(1, 12))

# Convert to raster for virtualspecies
bios <- stack(as(split(bios), 'Spatial'))

# Formatting of the response functions
set.seed(10)
my.parameters <- formatFunctions(
  bio1 = c(fun = 'dnorm', mean = 22, sd = 5),
  bio12 = c(fun = 'dnorm', mean = 1000, sd = 200))

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

# Pseudo presence-absence observations
set.seed(10)
po.points <- sampleOccurrences(
  my.species,
  n = 500,
  type = "presence-absence")
occ_virtual_species <- po.points$sample.points %>%
  dplyr::select(x, y, Observed) %>%
  rename(observation = Observed)

# Split the table into train and evaluation
# so the sample data can fit to multiple types of
# input
occ_virtual_species <- occ_virtual_species %>%
  mutate(id = row_number())

set.seed(11)
occ <- occ_virtual_species %>% sample_frac(0.7) %>%
  mutate(usage = "train")
occ_test <- occ_virtual_species %>% filter(! id %in% occ$id) %>%
  mutate(usage = "eval")
occ_virtual_species <- rbind(occ, occ_test) %>%
  select(-id); rm(occ, occ_test)

# Save out
usethis::use_data(occ_virtual_species, overwrite = TRUE)
