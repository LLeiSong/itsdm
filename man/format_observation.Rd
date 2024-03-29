% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_observation.R
\name{format_observation}
\alias{format_observation}
\title{Format the occurrence dataset for usage in \pkg{itsdm}}
\usage{
format_observation(
  obs_df,
  eval_df = NULL,
  split_perc = 0.3,
  seed = 123,
  obs_crs = 4326,
  eval_crs = 4326,
  x_col = "x",
  y_col = "y",
  obs_col = "observation",
  obs_type = "presence_only"
)
}
\arguments{
\item{obs_df}{(\code{data.frame}). The \code{data.frame} style table that
include x and y coordinate and observation of training dataset.
This parameter is required as it is the training dataset.
Note: it only takes \code{data.frame} to reduce the risk of column name mismatch
between \code{data.frame} and other formats such as \code{tibble}.}

\item{eval_df}{(\code{data.frame} or \code{NULL}) The \code{data.frame} style table that
include x and y coordinate and observation of evaluation dataset.
Note: it only takes \code{data.frame} to reduce the risk of column name mismatch
between \code{data.frame} and other formats such as \code{tibble}.}

\item{split_perc}{(\code{numeric}) a \code{numeric} between 0 and 1 corresponding to the
percentage of data used to evaluate the models.
Only required if \code{eval_df} is \code{NULL}.}

\item{seed}{(\code{integer}) The seed to split train and evaluation set.
The default value is \code{123}. Only required if \code{eval_df} is \code{NULL}.}

\item{obs_crs}{(\code{integer}, \code{numeric}, \code{character}, or \code{crs})
The EPSG code, CRS string, or \code{sf::crs} object of the coordinate system of
the training dataset. It corresponds to \code{x_col} and \code{y_col} in \code{obs_df}.}

\item{eval_crs}{(\code{integer}, \code{numeric}, \code{character}, or \code{crs})
The EPSG code, CRS string, or \code{sf::crs} object of the coordinate system
of the evaluation dataset. Only required if  \code{eval_df} is not \code{NULL}.
It corresponds to \code{x_col} and \code{y_col} in \code{eval_df} if any.}

\item{x_col}{(\code{character}) The name of column that is x coordinate
in \code{obs_df} and \code{eval_df} if not \code{NULL}.}

\item{y_col}{(\code{character}) The name of column that is y coordinate
in \code{obs_df} and \code{eval_df} if not \code{NULL}.}

\item{obs_col}{(\code{character}) The name of column that represents observations
in \code{obs_df} and \code{eval_df} if not \code{NULL}.}

\item{obs_type}{(\code{character}) The type of observation to be formatted to.
Only can be one of \code{c("presence_only", "presence_absence")}.
Note that if "presence_only" is set, the absences in \code{obs_df} will be deleted.
This only affect \code{obs_df}, \code{eval_df} will keep the original type no matter it
is an independent one or is split from \code{eval_df}.}
}
\value{
(\code{FormatOccurrence}) A list of
\itemize{
\item{obs (\code{sf}) the formatted pts of observations.
The column of observation is "observation".}
\item{obs_type (\code{character}) the type of the observations,
presence_only or presence_absence.}
\item{has_eval (\code{logical}) whether evaluation dataset is set or generated.}
\item{eval (\code{sf}) the formatted pts of observations for evaluation if any.
The column of observation is "observation".}
\item{eval (\code{eval_type})  the type of the observations for evaluation,
presence_only or presence_absence.}
}
}
\description{
The focus of this function is to format the dataset
but to keep the dataset as original as possible.
Then the users can modify the data if they want before put it
into this function.
}
\examples{
library(dplyr)
library(itsdm)
data("occ_virtual_species")

# obs + eval, presence-absence
obs_df <- occ_virtual_species \%>\% filter(usage == "train")
eval_df <- occ_virtual_species \%>\% filter(usage == "eval")
x_col <- "x"
y_col <- "y"
obs_col <- "observation"
obs_type <- "presence_absence"

obs <- format_observation(
  obs_df = obs_df, eval_df = eval_df,
  x_col = x_col, y_col = y_col, obs_col = obs_col,
  obs_type = obs_type)

# obs + eval, presence-only
obs_df <- occ_virtual_species \%>\% filter(usage == "train")
eval_df <- occ_virtual_species \%>\% filter(usage == "eval")
x_col <- "x"
y_col <- "y"
obs_col <- "observation"
obs_type <- "presence_only"

obs <- format_observation(
  obs_df = obs_df, eval_df = eval_df,
  x_col = x_col, y_col = y_col, obs_col = obs_col,
  obs_type = obs_type)

# obs + eval, different crs, presence-only
obs_df <- occ_virtual_species \%>\% filter(usage == "train")
eval_df <- occ_virtual_species \%>\% filter(usage == "eval")
obs_crs <- 4326
# Fake one
eval_crs <- 20935
x_col <- "x"
y_col <- "y"
obs_col <- "observation"
obs_type <- "presence_only"

obs <- format_observation(
  obs_df = obs_df, eval_df = eval_df,
  obs_crs = obs_crs, eval_crs = eval_crs,
  x_col = x_col, y_col = y_col, obs_col = obs_col,
  obs_type = obs_type)

# obs + split, presence-absence
obs_df <- occ_virtual_species
split_perc <- 0.5
seed <- 123
obs_crs <- 4326
x_col <- "x"
y_col <- "y"
obs_col <- "observation"
obs_type <- "presence_absence"

obs <- format_observation(
  obs_df = obs_df, split_perc = split_perc,
  x_col = x_col, y_col = y_col,
  obs_col = obs_col, obs_type = obs_type)

# obs, presence-only, no eval
obs_df <- occ_virtual_species
eval_df <- NULL
split_perc <- 0
seed <- 123
obs_crs <- 4326
x_col <- "x"
y_col <- "y"
obs_col <- "observation"
obs_type <- "presence_only"

obs <- format_observation(
  obs_df = obs_df, eval_df = eval_df,
  split_perc = split_perc,
  x_col = x_col, y_col = y_col,
  obs_col = obs_col, obs_type = obs_type)

}
\seealso{
\code{\link{print.FormatOccurrence}}
}
