#' @title Calculate independent responses of each variables.
#' @description Calculate the independent responses of each variables
#' within the model.
#' @param model (Any predictive model). It is `isolation_forest` here. It could
#' be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' @param var_occ (`data.frame`, `tibble`) The `data.frame` style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (`stars`) The `stars` of environmental variables.
#' It should have multiple `attributes` instead of `dims`.
#' If you have `raster` object instead, you
#' could use \code{\link{st_as_stars}} to convert it to `stars` or use
#' \code{\link{read_stars}} directly read source data as a `stars`.
#' You also could use item `variables` of `POIsotree` made by function
#' \code{\link{isotree_po}}.
#' @param si (`integer`) The number of samples to generate response curves.
#' If it is too small, the response curves might be biased.
#' The default value is `1000`.
#' @param visualize (`logical`) if `TRUE`, plot the response curves.
#' The default is `FALSE`.
#'
#' @return (`IndependentResponse`) A list of
#' \itemize{
#' \item{responses_cont (`list`) A list of response values of continuous
#' variables}
#' \item{responses_cat (`list`) A list of response values of categorical
#' variables}
#' }
#'
#' @details
#' The values show how each environmental variable independently
#' affects the modeling prediction. They show how the predicted result
#' only using this variable changes as it is varied.
#'
#' @references
#' \itemize{
#' \item{Elith, Jane,
#' et al. "The evaluation strip: a new and robust method for plotting predicted
#' responses from species distribution models." \emph{Ecological modelling}
#' 186.3 (2005): 280-289.\doi{10.1016/j.ecolmodel.2004.12.007}}
#' }
#'
#' @seealso
#' \code{\link{plot.IndependentResponse}}
#'
#' @importFrom isotree isolation.forest
#' @importFrom dplyr select slice
#' @importFrom stars st_as_stars
#' @importFrom stats predict setNames
#' @importFrom rlang :=
#' @importFrom tidyselect all_of
#' @export
#' @examples
#' # Using a pseudo presence-only occurrence dataset of
#' # virtual species provided in this package
#' library(dplyr)
#' library(sf)
#' library(stars)
#' library(itsdm)
#'
#' data("occ_virtual_species")
#' obs_df <- occ_virtual_species %>% filter(usage == "train")
#' eval_df <- occ_virtual_species %>% filter(usage == "eval")
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#'
#' # Format the observations
#' obs_train_eval <- format_observation(
#'   obs_df = obs_df, eval_df = eval_df,
#'   x_col = x_col, y_col = y_col, obs_col = obs_col,
#'   obs_type = "presence_only")
#'
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 30,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' independent_responses <- independent_response(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables)
#' plot(independent_responses)
#'
independent_response <- function(model,
                                 var_occ,
                                 variables,
                                 si = 1000,
                                 visualize = FALSE) {

  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_int(si)
  checkmate::assert_class(variables, 'stars')
  stopifnot(length(dim(variables)) <= 2)
  checkmate::assert_logical(visualize)
  bands <- names(variables)
  stopifnot(all(bands %in% colnames(var_occ)))

  # Do full prediction
  ## Raster
  var_pred_full <- predict(variables, model)
  ## Stretch result to be comparable with other results
  var_pred_full <- 1 - var_pred_full

  # Reformat data
  ## Make models
  models <- lapply(bands, function(nm) {
    ind_var_occ <- var_occ %>% select(all_of(nm))
    # Remove and modify some arguments not works
    # for single-variable model, other arguments
    # inherit from model object.
    args_iforest <- c(
      list(data = ind_var_occ, seed = model$random_seed,
           nthreads = model$nthreads),
      model$params
    )
    args_iforest$ndim <- 1
    args_iforest$ncols_per_tree <- min(ncol(args_iforest$data),
                                       args_iforest$ncols_per_tree)
    if (args_iforest$new_categ_action == "impute") {
      args_iforest$new_categ_action <- "weighted"
    }
    do.call(isolation.forest, args_iforest)}) %>%
    setNames(bands)

  ## Split numeric and categorical
  isfacor <- as.vector(sapply(variables, is.factor))
  bands_cont <- bands[!isfacor]
  bands_cat <- bands[isfacor]

  # Numeric variables
  ## Not limited to data volume, could generate as many as possible
  ## pseudo observations in [min, max], so the function could make
  ## smoother curves.
  if (length(bands_cont) > 0) {
    mins <- sapply(bands_cont, function(nm) {
      min(variables %>% select(nm) %>% pull, na.rm = T)})
    maxs <- sapply(bands_cont, function(nm) {
      max(variables %>% select(nm) %>% pull, na.rm = T)})
    vals_cont <- do.call(cbind, lapply(1:length(mins), function(n) {
      seq(from = mins[n], to = maxs[n],
          length.out = si)})) %>%
      data.frame() %>%
      setNames(bands_cont)
    responses_cont <- lapply(names(vals_cont), function(nm) {
      vals_tmp <- vals_cont %>% select(nm)
      pred_tmp <- predict(models[[nm]], vals_tmp)
      pred_tmp <- 1 - pred_tmp
      pred_tmp <- .stars_stretch(var_pred_full, new_values = pred_tmp)
      data.frame(y = pred_tmp,
                 x = vals_cont %>% pull(nm)) %>%
        setNames(c("y", "x"))
    })
    names(responses_cont) <- bands_cont
    rm(mins, maxs, vals_cont)
  } else responses_cont <- NULL

  # Categorical variables
  ## The number of pseudo observations is limited to factor levels
  ## So categorical variables should generate one by one
  if (length(bands_cat) > 0) {
    responses_cat <- lapply(bands_cat, function(nm) {
      vals_this <- var_occ %>% pull(nm) %>%
        levels() %>% as.factor() %>%
        data.frame() %>% setNames(nm)

      pred_tmp <- predict(models[[nm]], vals_this)
      pred_tmp <- 1 - pred_tmp
      pred_tmp <- .stars_stretch(var_pred_full, new_values = pred_tmp)
      data.frame(y = pred_tmp,
                 x = vals_this %>% pull(nm)) %>%
        setNames(c("y", "x"))
    }) %>% setNames(bands_cat)
  } else responses_cat <- NULL

  # Put together
  responses <- list(responses_cont = responses_cont,
                    responses_cat = responses_cat)
  class(responses) <- append("IndependentResponse", class(responses))

  # Visualize
  if (visualize) {
    print(plot(responses))
  }

  # Return
  return(responses)
}

# independent_response end
