# itsdm <img src='man/figures/hexagon_sticker.png' align="right" height="120"/>

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--11--16-yellowgreen.svg)](/commits/main)
<!-- badges: end -->

## Overview

`itsdm` calls isolation forest and variations such as SCiForest and EIF to do species distribution modeling. It provides features including:

- A few function to download environmental variables.
- Outlier tree based suspicious environmental outliers detection.
- Isolation forest based environmental suitability modeling.
- Response curves of environmental variable.
- Variable importance analysis.
- Presence-only model evaluation.
- Method to convert predicted suitability to presence-absence map.
- Variable contribution analysis for interested observations.

## Installation

You can install the development version of itsdm from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LLeiSong/itsdm")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(itsdm)

# Using a pseudo presence-only occurrence dataset of
# virtual species provided in this package
data("occ_virtual_species")

# Split to training and test
occ_virtual_species <- occ_virtual_species %>%
  mutate(id = row_number())
set.seed(11)
occ <- occ_virtual_species %>% sample_frac(0.7)
occ_test <- occ_virtual_species %>% filter(! id %in% occ$id)
occ <- occ %>% select(-id)
occ_test <- occ_test %>% select(-id)

# Get environmental variables
env_vars <- system.file(
  'extdata/bioclim_africa_10min.tif',
  package = 'itsdm') %>% read_stars() %>%
  %>% slice('band', c(1, 12))

# Train the model
mod <- isotree_po(
  occ = occ, occ_test = occ_test,
  variables = env_vars, ntrees = 200,
  sample_rate = 0.8, ndim = 0L,
  seed = 123L)

# Check results
## Suitability
ggplot() +
  geom_stars(data = mod$prediction) +
  scale_fill_viridis_c('Predicted suitability',
                       na.value = 'transparent') +
  coord_equal() +
  theme_linedraw()

## Plot independent response curves
plot(it_sdm$independent_responses, 
  target_var = c('bio1', 'bio12', 'bio15'))
```
