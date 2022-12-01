# itsdm <img src='man/figures/hexagon_sticker.png' align="right" height="120"/>

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/LLeiSong/itsdm/workflows/R-CMD-check/badge.svg)](https://github.com/LLeiSong/itsdm/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/itsdm)](https://CRAN.R-project.org/package=itsdm)
<!-- badges: end -->

## Overview

`itsdm` calls isolation forest and variations such as SCiForest and EIF to model species distribution. It provides features including:

- A few functions to download environmental variables.
- Outlier tree-based suspicious environmental outliers detection.
- Isolation forest-based environmental suitability modeling.
- Non-spatial response curves of environmental variables.
- Spatial response maps of environmental variables.
- Variable importance analysis.
- Presence-only model evaluation.
- Method to convert predicted suitability to presence-absence map.
- Variable contribution analysis for the target observations.
- Method to analyze the spatial impacts of changing environment.

## Installation

Install the CRAN release of `itsdm` with

```r
install.packages("itsdm")
```

You can install the development version of itsdm from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("LLeiSong/itsdm")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(itsdm)
library(dplyr)
library(stars)
library(ggplot2)

# Using a pseudo presence-only occurrence dataset of
# virtual species provided in this package
data("occ_virtual_species")
obs_df <- occ_virtual_species %>% filter(usage == "train")
eval_df <- occ_virtual_species %>% filter(usage == "eval")
x_col <- "x"
y_col <- "y"
obs_col <- "observation"
obs_type <- "presence_absence"

# Format the observations
obs_train_eval <- format_observation(
  obs_df = obs_df, eval_df = eval_df,
  x_col = x_col, y_col = y_col, obs_col = obs_col,
  obs_type = obs_type)

# Get environmental variables
env_vars <- system.file(
  'extdata/bioclim_tanzania_10min.tif',
  package = 'itsdm') %>% read_stars() %>%
  slice('band', c(1, 6, 12, 15))

# Train the model
mod <- isotree_po(
  obs_mode = "presence_absence",
  obs = obs_train_eval$obs,
  obs_ind_eval = obs_train_eval$eval,
  variables = env_vars, ntrees = 200,
  sample_size = 0.8, ndim = 2,
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
plot(mod$independent_responses, 
     target_var = c('bio1', 'bio12'))
```

The Shapley values-based analysis can apply to external models. Here is an example to analyze impacts of the bio12 decreasing 200 mm to species distribution based on Random Forest (RF) prediction:

```
# Prepare data
data("occ_virtual_species")
obs_df <- occ_virtual_species %>% 
  filter(usage == "train")

env_vars <- system.file(
  'extdata/bioclim_tanzania_10min.tif',
  package = 'itsdm') %>% read_stars() %>%
  slice('band', c(1, 5, 12)) %>% 
  split()

model_data <- stars::st_extract(
  env_vars, at = as.matrix(obs_df %>% select(x, y))) %>% 
  as.data.frame()
names(model_data) <- names(env_vars)
model_data <- model_data %>% 
  mutate(occ = obs_df[['observation']])
model_data$occ <- as.factor(model_data$occ)

mod_rf <- randomForest(
  occ ~ .,
  data = model_data,
  ntree = 200)

pfun <- function(X.model, newdata) {
  # for data.frame
  predict(X.model, newdata, type = "prob")[, "1"]
}

# Use a fixed value
climate_changes <- detect_envi_change(
  model = mod_rf,
  var_occ = model_data %>% select(-occ),
  variables = env_vars,
  target_var = "bio12",
  bins = 20,
  var_future = -200,
  pfun = pfun)
```

## Contributor

1. [David Cortes](https://github.com/david-cortes), helps to improve the flexibility of calling `isotree`.

We are welcome any helps! Please make a pull request or reach out to lsong@clarku.edu if you want to make any contribution.

## Funding
This package is part of project "Combining Spatially-explicit Simulation of Animal Movement and Earth Observation to Reconcile Agriculture and Wildlife Conservation". This project is funded by NASA FINESST program (award number: 80NSSC20K1640).
