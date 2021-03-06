% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.VariableContribution}
\alias{plot.VariableContribution}
\title{Function to plot variable contribution for target observations.}
\usage{
\method{plot}{VariableContribution}(x, plot_each_obs = FALSE, num_features = 5, ...)
}
\arguments{
\item{x}{(\code{VariableContribution}) The \code{VariableContribution} object to plot.
It could be the return of function \code{\link{variable_contrib}}.}

\item{plot_each_obs}{(\code{logical}) The option of plot type. If \code{TRUE}, it will
plot variable contribution for every observation. Otherwise, it will plot
variable contribution violin plot for all observations.}

\item{num_features}{(\code{integer}) A number of most important features to plot.
Just work if plot_each_obs is \code{TRUE}.}

\item{...}{Not used.}
}
\value{
\code{ggplot2} figure of Variable Contribution.
}
\description{
Plot variable contribution for target observation separately
or together using ggplot2.
}
\examples{
\donttest{
# Using a pseudo presence-only occurrence dataset of
# virtual species provided in this package
library(dplyr)
library(sf)
library(stars)
library(itsdm)

data("occ_virtual_species")
occ_virtual_species <- occ_virtual_species \%>\%
  mutate(id = row_number())

set.seed(11)
occ <- occ_virtual_species \%>\% sample_frac(0.7)
occ_test <- occ_virtual_species \%>\% filter(! id \%in\% occ$id)
occ <- occ \%>\% select(-id)
occ_test <- occ_test \%>\% select(-id)

env_vars <- system.file(
  'extdata/bioclim_tanzania_10min.tif',
  package = 'itsdm') \%>\% read_stars() \%>\%
  slice('band', c(1, 5, 12, 16))

mod <- isotree_po(
  occ = occ, occ_test = occ_test,
  variables = env_vars, ntrees = 50,
  sample_size = 0.8, ndim = 1L,
  seed = 123L, response = FALSE,
  spatial_response = FALSE,
  check_variable = FALSE)

var_contribution <- variable_contrib(
  model = mod$model,
  var_occ = mod$var_train \%>\% st_drop_geometry(),
  var_occ_analysis = mod$var_train \%>\%
    st_drop_geometry() \%>\% slice(1:10))
plot(var_contribution,
  plot_each_obs = TRUE,
  num_features = 3)
plot(var_contribution)
}

}
\seealso{
\code{\link{variable_contrib}}
}
