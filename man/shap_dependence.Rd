% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shap_dependence.R
\name{shap_dependence}
\alias{shap_dependence}
\title{Calculate Shapley value based variable dependence.}
\usage{
shap_dependence(model, var_occ, shap_nsim = 100, visualize = FALSE, seed = 10L)
}
\arguments{
\item{model}{(\code{isolation_forest}). The isolation forest SDM.
It could be the item \code{model} of \code{POIsotree} made by function \code{\link{isotree_po}}.}

\item{var_occ}{(\code{data.frame}, \code{tibble}) The \code{data.frame} style table that
include values of environmental variables at occurrence locations.}

\item{shap_nsim}{(\code{integer}) The number of Monte Carlo repetitions in SHAP
method to use for estimating each Shapley value. When the number of variables
is large, a smaller shap_nsim could be used. See details in documentation of
function \code{\link{explain}} in package \code{fastshap}.
The default is 100.}

\item{visualize}{(\code{logical}) if \code{TRUE}, plot the variable dependence plots.
The default is \code{FALSE}.}

\item{seed}{(\code{integer}) The seed for any random progress. The default is \code{10L}.}
}
\value{
(\code{ShapDependence}) A list of
\itemize{
\item{dependences_cont (\code{list}) A list of Shapley values of continuous variables}
\item{dependences_cat (\code{list}) A list of Shapley values of categorical variables}
\item{feature_values (\code{data.frame}) A table of feature values}
}
}
\description{
Calculate the variable dependence using Shapley values.
}
\details{
The values show how each environmental variable independently
affects the modeling prediction. They show how the Shapley value of each variable
changes as its value is varied.
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

var_dependence <- shap_dependence(
  model = mod$model,
  var_occ = mod$var_train \%>\% st_drop_geometry())
}

}
\references{
\itemize{
\item{Strumbelj, Erik,
and Igor Kononenko. "Explaining prediction models and individual predictions
with feature contributions." \emph{Knowledge and information systems}
41.3 (2014): 647-665.\doi{10.1007/s10115-013-0679-x}}
\item{\href{http://proceedings.mlr.press/v119/sundararajan20b.html}{Sundara
rajan, Mukund, and Amir Najmi. "The many Shapley values for model explanation
." \emph{International Conference on Machine Learning}. PMLR, 2020.}}
\item{\url{https://github.com/bgreenwell/fastshap}}
\item{\url{https://github.com/slundberg/shap}}
}
}
\seealso{
\code{\link{plot.ShapDependence}}
\code{\link{explain}} in \code{fastshap}
}
