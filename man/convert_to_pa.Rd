% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_to_pa.R
\name{convert_to_pa}
\alias{convert_to_pa}
\title{Convert predicted suitability to presence-absence map.}
\usage{
convert_to_pa(
  suitability,
  method = "logistic",
  beta = 0.5,
  alpha = -0.05,
  a = 1,
  b = 0,
  species_prevalence = NA,
  threshold = 0.5,
  seed = 10L,
  visualize = TRUE
)
}
\arguments{
\item{suitability}{(\code{stars} or \code{RasterLayer}) The suitability raster.}

\item{method}{(\code{character}) The conversion method, must be one of
'threshold', 'logistic', and 'linear'. The default is 'logistic'.}

\item{beta}{(\code{numeric}) Works for 'threshold' or 'logistic' method.
If \code{method} is threshold, then \code{beta} is the threshold value to cutoff.
If \code{method} is logistic, it is the sigmoid midpoint. The default is \code{0.5}.}

\item{alpha}{(\code{numeric}) Works for logistic method.
It is the logistic growth rate or steepness of the curve.
The default is \code{-.05}.}

\item{a}{(\code{numeric}) Works for linear method. It is the slope of the line.
The default is \code{1}.}

\item{b}{(\code{numeric}) Works for linear method.
It is the intercept of the line. The default is \code{0}.}

\item{species_prevalence}{(\code{numeric} or \code{NA}) Works for all three methods.
It is the species prevalence to classify suitability map.
It could be \code{NA}, when the will be calculated automatically
based on other arguments. The default is \code{NA}.}

\item{threshold}{(\code{numeric}) The threshold used to convert probability of
occurrence to presence-absence map. It ranges in \verb{[0, 1]}. The default is 0.5.}

\item{seed}{(\code{integer}) The seed for random progress. The default is \code{10L}}

\item{visualize}{(\code{logical}) If \code{TRUE}, plot map of suitability,
probability of occurrence, and presence-absence together.
The default is \code{TRUE}.}
}
\value{
(\code{PAConversion}) A list of
\itemize{
\item{suitability (\code{stars}) The input suitability map}
\item{probability_of_occurrence (\code{stars}) The map of occurrence probability}
\item{pa_conversion (\code{list}) A list of conversion arguments}
\item{pa_map (\code{stars}) The presence-absence map}}
}
\description{
Use threshold-based, logistic or linear conversion method to
convert predicted suitability map to presence-absence map.
}
\details{
Multiple methods and arguments could be used as a combination to do
the conversion.
}
\examples{
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
  sample_size = 0.8, ndim = 3L,
  seed = 123L, response = FALSE,
  spatial_response = FALSE,
  check_variable = FALSE)

# Threshold conversion
pa_thred <- convert_to_pa(mod$prediction,
  method = 'threshold', beta = 0.5)
pa_thred
plot(pa_thred)

# Logistic conversion
pa_log <- convert_to_pa(mod$prediction, method = 'logistic',
  beta = 0.5, alpha = -.05)

# Linear conversion
pa_lin <- convert_to_pa(mod$prediction, method = 'linear',
  a = 1, b = 0)

}
\references{
\href{https://github.com/Farewe/virtualspecies/blob/master/R/convertToPA.R}{c
onvertToPA in package \code{virtualspecies}}
}
\seealso{
\code{\link{plot.PAConversion}}
}
