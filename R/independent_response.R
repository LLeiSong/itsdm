#' A function to calculate independent responses of each variables.
#' @description This function allows you to calculate the independent responses
#' of each variables within the model.
#' @param var_occ (data.frame, tibble) the data.frame style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (RasterStack, RasterLayer or stars)
#' the stack of environmental variables. Could be raster or stars object.
#' @param si (integer) the number of samples to generate response curves.
#' If it is too small, the response curves might be biased.
#' The default value is 1000.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @details The curves show how each environmental variable independently
#' affects the modeling prediction. The curves show how the predicted result
#' only using this variable changes as it is varied.
#' @return a list of data.frame of response and variable values.
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
  checkmate::assert_multi_class(
    variables, c('RasterStack', 'RasterLayer', 'stars'))
  checkmate::assert_logical(visualize)
  bands <- st_get_dimension_values(variables, 'band')
  stopifnot(all(bands %in% colnames(var_occ)))

  # Reformat data
  # Variables, use stars
  if (is(variables, 'RasterStack') |
      is(variables, 'RasterLayer')){
    variables <- st_as_stars(variabels)}

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
      weigh_imp_rows = model$params$weigh_imp_rows)})

  mins <- sapply(1:length(bands), function(i) {
    min(variables %>% slice('band', i) %>% unlist(), na.rm = T)})
  maxs <- sapply(1:length(bands), function(i) {
    max(variables %>% slice('band', i) %>% unlist(), na.rm = T)})
  vals <- do.call(cbind, lapply(1:length(mins), function(n) {
    seq(from = mins[n], to = maxs[n],
        length.out = si)})) %>%
    data.frame() %>%
    setNames(bands)

  responses <- lapply(1:ncol(vals), function(x) {
    vals_tmp <- vals %>% select(bands[x])
    pred_tmp <- predict(models[[x]], vals_tmp)
    pred_tmp <- 1 - pred_tmp
    pred_tmp <- .norm(pred_tmp)
    data.frame(y = pred_tmp, x = vals_tmp) %>%
      setNames(c("y", "x"))
  })
  names(responses) <- bands
  class(responses) <- append("independent_response", class(responses))

  # Visualize
  if (visualize) {
    plot(responses)
  }

  # Return
  return(responses)
}

# independent_response end
