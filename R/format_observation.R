#' @title Format the occurrence dataset for usage in \pkg{itsdm}
#' @description The focus of this function is to format the dataset
#' but to keep the dataset as original as possible.
#' Then the users can modify the data if they want before put it
#' into this function.
#' @param obs_df (`data.frame`). The `data.frame` style table that
#' include x and y coordinate and observation of training dataset.
#' This parameter is required as it is the training dataset.
#' Note: it only takes `data.frame` to reduce the risk of column name mismatch
#' between `data.frame` and other formats such as `tibble`.
#' @param eval_df (`data.frame` or `NULL`) The `data.frame` style table that
#' include x and y coordinate and observation of evaluation dataset.
#' Note: it only takes `data.frame` to reduce the risk of column name mismatch
#' between `data.frame` and other formats such as `tibble`.
#' @param split_perc (`numeric`) a `numeric` between 0 and 1 corresponding to the
##' percentage of data used to evaluate the models.
##' Only required if `eval_df` is `NULL`.
#' @param seed (`integer`) The seed to split train and evaluation set.
#' The default value is `123`. Only required if `eval_df` is `NULL`.
#' @param obs_crs (`integer`, `numeric`, `character`, or `crs`)
#' The EPSG code, CRS string, or `sf::crs` object of the coordinate system of
#' the training dataset. It corresponds to `x_col` and `y_col` in `obs_df`.
#' @param eval_crs (`integer`, `numeric`, `character`, or `crs`)
#' The EPSG code, CRS string, or `sf::crs` object of the coordinate system
#' of the evaluation dataset. Only required if  `eval_df` is not `NULL`.
#' It corresponds to `x_col` and `y_col` in `eval_df` if any.
#' @param x_col (`character`) The name of column that is x coordinate
#' in `obs_df` and `eval_df` if not `NULL`.
#' @param y_col (`character`) The name of column that is y coordinate
#' in `obs_df` and `eval_df` if not `NULL`.
#' @param obs_col (`character`) The name of column that represents observations
#' in `obs_df` and `eval_df` if not `NULL`.
#' @param obs_type (`character`) The type of observation to be formatted to.
#' Only can be one of `c("presence_only", "presence_absence")`.
#' Note that if "presence_only" is set, the absences in `obs_df` will be deleted.
#' This only affect `obs_df`, `eval_df` will keep the original type no matter it
#' is an independent one or is split from `eval_df`.
#' @return (`FormatOccurrence`) A list of
#' \itemize{
#' \item{obs (`sf`) the formatted pts of observations.
#' The column of observation is "observation".}
#' \item{obs_type (`character`) the type of the observations,
#' presence_only or presence_absence.}
#' \item{has_eval (`logical`) whether evaluation dataset is set or generated.}
#' \item{eval (`sf`) the formatted pts of observations for evaluation if any.
#' The column of observation is "observation".}
#' \item{eval (`eval_type`)  the type of the observations for evaluation,
#' presence_only or presence_absence.}
#' }
#'
#' @seealso
#' \code{\link{print.FormatOccurrence}}
#'
#' @importFrom dplyr mutate row_number sample_frac select filter all_of rename
#' @importFrom sf st_as_sf st_crs st_transform
#' @export
#' @examples
#' library(dplyr)
#' library(itsdm)
#' data("occ_virtual_species")
#'
#' # obs + eval, presence-absence
#' obs_df <- occ_virtual_species %>% filter(usage == "train")
#' eval_df <- occ_virtual_species %>% filter(usage == "eval")
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#' obs_type <- "presence_absence"
#'
#' obs <- format_observation(
#'   obs_df = obs_df, eval_df = eval_df,
#'   x_col = x_col, y_col = y_col, obs_col = obs_col,
#'   obs_type = obs_type)
#'
#' # obs + eval, presence-only
#' obs_df <- occ_virtual_species %>% filter(usage == "train")
#' eval_df <- occ_virtual_species %>% filter(usage == "eval")
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#' obs_type <- "presence_only"
#'
#' obs <- format_observation(
#'   obs_df = obs_df, eval_df = eval_df,
#'   x_col = x_col, y_col = y_col, obs_col = obs_col,
#'   obs_type = obs_type)
#'
#' # obs + eval, different crs, presence-only
#' obs_df <- occ_virtual_species %>% filter(usage == "train")
#' eval_df <- occ_virtual_species %>% filter(usage == "eval")
#' obs_crs <- 4326
#' # Fake one
#' eval_crs <- 20935
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#' obs_type <- "presence_only"
#'
#' obs <- format_observation(
#'   obs_df = obs_df, eval_df = eval_df,
#'   obs_crs = obs_crs, eval_crs = eval_crs,
#'   x_col = x_col, y_col = y_col, obs_col = obs_col,
#'   obs_type = obs_type)
#'
#' # obs + split, presence-absence
#' obs_df <- occ_virtual_species
#' split_perc <- 0.5
#' seed <- 123
#' obs_crs <- 4326
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#' obs_type <- "presence_absence"
#'
#' obs <- format_observation(
#'   obs_df = obs_df, split_perc = split_perc,
#'   x_col = x_col, y_col = y_col,
#'   obs_col = obs_col, obs_type = obs_type)
#'
#' # obs, presence-only, no eval
#' obs_df <- occ_virtual_species
#' eval_df <- NULL
#' split_perc <- 0
#' seed <- 123
#' obs_crs <- 4326
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#' obs_type <- "presence_only"
#'
#' obs <- format_observation(
#'   obs_df = obs_df, eval_df = eval_df,
#'   split_perc = split_perc,
#'   x_col = x_col, y_col = y_col,
#'   obs_col = obs_col, obs_type = obs_type)
#'
# The focus of this function is to format the dataset.
# So it keeps the dataset as original as possible.
# The user can modify the dataset before put into this function
format_observation <- function(obs_df, # observation data frame
                               eval_df = NULL, # evaluation data frame
                               split_perc = 0.3, # percentage of obs for eval
                               seed = 123,
                               obs_crs = 4326, # crs of obs coordinates
                               eval_crs = 4326, # crs of eval coordinates
                               x_col = "x",
                               y_col = "y",
                               obs_col = "observation",
                               # Depends on users' need
                               obs_type = "presence_only") {

  # Check format of inputs
  checkmate::assert_data_frame(obs_df)
  checkmate::assert_data_frame(eval_df, null.ok = T)
  checkmate::assert_multi_class(
    obs_crs, c("integer", "numeric", 'character', 'crs'))
  checkmate::assert_character(x_col)
  checkmate::assert_character(y_col)
  checkmate::assert_character(obs_col)
  checkmate::assert_choice(
    obs_type, choices = c("presence_only", "presence_absence"))

  # Detailed check
  ## Check observation table
  if (!all(c(x_col, y_col, obs_col) %in% names(obs_df))){
    stop("Columns not found in observation data table.")}

  # Check the evaluation table
  if (is.null(eval_df)) {
    message(paste0('No independent evaluation. Will read split_perc and seed.'))
    checkmate::assert_number(split_perc, lower = 0, upper = 1)
    checkmate::assert_int(seed)
  } else {
    checkmate::assert_multi_class(
      eval_crs, c("integer", "numeric", 'character', 'crs'))
    if (!all(c(x_col, y_col, obs_col) %in% names(eval_df))){
      stop("Columns not found in evaluation data table.")}
  }

  # No problem in the inputs, start the process
  ## Step 1: specify train and test set
  ## Step 2: convert the train and test set into spatial format
  # Split the dataset
  has_eval <- TRUE
  if (is.null(eval_df)) {
    if (split_perc > 0) {
      eval_crs <- obs_crs
      obs_df <- obs_df %>%
        mutate(id = row_number())

      set.seed(seed)
      occ <- obs_df %>% sample_frac(split_perc)
      occ_left <- obs_df %>% filter(!.data$id %in% occ$id)
      obs_df <- occ_left %>% select(-c("id"))
      eval_df <- occ %>% select(-c("id"))
      rm(occ, occ_left)
    } else {
      warning("No evaluation set to use.")
      has_eval <- FALSE}
  }

  # train observation
  eval_type <- obs_type
  if (obs_type == "presence_only") {
    if (sum(obs_df[, obs_col] == 0) > 0) {
      warning(paste0("Set observation type as presence-only, ",
                     "extra absences may be deleted."))
      obs_df <- obs_df %>% filter(!!as.symbol(obs_col) == 1)}}
  obs_pts <- obs_df %>%
    select(all_of(x_col), all_of(y_col), all_of(obs_col)) %>%
    st_as_sf(coords = c(x_col, y_col), crs = obs_crs) %>%
    rename("observation" = all_of(obs_col))

  # eval
  ## Check the observation type first
  if (!is.null(eval_df)) {
    if (sum(eval_df[, obs_col] == 0) > 0) {
      eval_type <- "presence_absence"
    }
    eval_pts <- eval_df %>%
      select(all_of(x_col), all_of(y_col), all_of(obs_col)) %>%
      st_as_sf(coords = c(x_col, y_col), crs = eval_crs) %>%
      rename("observation" = all_of(obs_col))

    # Transform if need
    if (st_crs(eval_pts) != st_crs(obs_pts)) {
      eval_pts <- eval_pts %>% st_transform(crs = st_crs(obs_pts))
    }
  } else eval_pts <- NULL

  ## Gather and return the result
  ### So within the formatted sf, the column name of observation
  ### is observation
  out <- list(obs = obs_pts,
              obs_type = obs_type,
              has_eval = has_eval,
              eval = eval_pts,
              eval_type = if(is.null(eval_pts)) NULL else eval_type)

  class(out) <- append("FormatOccurrence", class(out))

  return(out)
}
