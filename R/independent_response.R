#' @title A function to calculate independent responses of each variables.
#' @description This function allows you to calculate the independent responses
#' of each variables within the model.
#' @param model (isolation_forest) The fitted SDM model.
#' @param var_occ (data.frame, tibble) the data.frame style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (stars) The stars of environmental variables. It should have
#' multiple attributes instead of dims. If you have `raster` object instead, you
#' could use `st_as_stars` to convert it to `stars` or use `read_stars` directly
#' read source data as a `stars`.
#' @param si (integer) the number of samples to generate response curves.
#' If it is too small, the response curves might be biased.
#' The default value is 1000.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @details The curves show how each environmental variable independently
#' affects the modeling prediction. The curves show how the predicted result
#' only using this variable changes as it is varied.
#' @return (IndependentResponse) a list of data.frame of response and variable values.
#' The response values correspond to suitability of this single variable.
#' @importFrom isotree isolation.forest
#' @importFrom dplyr select slice
#' @importFrom stars st_as_stars
#' @export
#' @examples
#' independent_response(model = mod$model, variables = env_vars)
independent_response <- function(model,
                                 var_occ,
                                 variables,
                                 si = 1000,
                                 visualize = FALSE) {

  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_int(si)
  checkmate::assert_class(variables, 'stars')
  stopifnot(length(dim(variables)) == 2)
  checkmate::assert_logical(visualize)
  bands <- names(variables)
  stopifnot(all(bands %in% colnames(var_occ)))

  # Reformat data
  ## Make models
  models <- lapply(bands, function(nm) {
    ind_var_occ <- var_occ %>% select(all_of(nm))
    # Remove and modify some arguments not works
    # for single-variable model, other arguments
    # inherit from model object.
    isolation.forest(
      ind_var_occ,
      ntrees = model$params$ntrees,
      sample_size = model$params$sample_size,
      ndim = 1,
      ntry = model$params$ntry,
      max_depth = model$params$max_depth,
      prob_pick_avg_gain = model$params$prob_pick_avg_gain,
      prob_pick_pooled_gain = model$params$prob_pick_pooled_gain,
      prob_split_avg_gain = model$params$prob_split_avg_gain,
      prob_split_pooled_gain = model$params$prob_split_pooled_gain,
      min_gain = model$params$min_gain,
      missing_action = model$params$missing_action,
      categ_split_type = model$params$categ_split_type,
      all_perm = model$params$all_perm,
      coef_by_prop = model$params$coef_by_prop,
      weights_as_sample_prob = model$params$weights_as_sample_prob,
      sample_with_replacement = model$params$sample_with_replacement,
      penalize_range = model$params$penalize_range,
      weigh_by_kurtosis = model$params$weigh_by_kurtosis,
      coefs = model$params$coefs,
      assume_full_distr = model$params$assume_full_distr,
      build_imputer = model$params$build_imputer,
      min_imp_obs = model$params$min_imp_obs,
      depth_imp = model$params$depth_imp,
      weigh_imp_rows = model$params$weigh_imp_rows)}) %>%
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
      pred_tmp <- .norm(pred_tmp)
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
      pred_tmp <- .norm(pred_tmp)
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
    plot(responses)
  }

  # Return
  return(responses)
}

# independent_response end
