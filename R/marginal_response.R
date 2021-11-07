#' @title A function to calculate marginal responses of each variables.
#' @description This function allows you to calculate the marginal responses of
#' each variables within the model.
#' @param model Any predictive model. In this package, it could be isolation_forest.
#' In this package, it is isolation_forest.
#' @param var_occ (data.frame, tibble) the data.frame style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (RasterStack, RasterLayer or stars)
#' the stack of environmental variables. Could be raster or stars object.
#' @param si (integer) the number of samples to generate response curves.
#' If it is too small, the response curves might be biased.
#' The default value is 1000.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @details The curves show how each environmental variable affects the modeling
#' prediction. The curves show how the predicted result changes as each environmental
#' variable is varied while keeping all other environmental variables at average
#' sample value. The curves might be hard to interpret if there are strongly correlated
#' variables. The users could use `dim_reduce` function to remove the strong correlation
#' from original environmental variable stack.
#' @return (MarginalResponse) a list of data.frame of response and variable values.
#' The response values correspond to suitability of this single variable.
#' @importFrom dplyr select slice as_tibble pull
#' @importFrom stars st_as_stars
#' @export
#' @examples
#' marginal_response(model = mod$model, variables = env_vars)
marginal_response <- function(model,
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

  mins <- sapply(1:length(bands), function(i) {
    min(variables %>% slice('band', i) %>% unlist(), na.rm = T)})
  maxs <- sapply(1:length(bands), function(i) {
    max(variables %>% slice('band', i) %>% unlist(), na.rm = T)})
  vals <- do.call(cbind, lapply(1:length(mins), function(n) {
    seq(from = mins[n], to = maxs[n],
      length.out = si)})) %>%
    data.frame() %>%
    setNames(bands)

  # Get means and expand
  means <- lapply(bands, function(nm) {
    mean(var_occ %>% pull(nm), na.rm = T)}) %>%
    setNames(bands) %>%
    as_tibble() %>%
    slice(rep(1:n(), each = nrow(vals)))
  responses <- lapply(1:ncol(vals), function(x) {
    vals_tmp <- means
    vals_tmp[, x] <- vals[, x]
    pred_tmp <- predict(model, as.matrix(vals_tmp))
    pred_tmp <- 1 - pred_tmp
    pred_tmp <- .norm(pred_tmp)
    data.frame(y = pred_tmp, x = vals[, x]) %>%
      setNames(c("y", "x"))
  })
  names(responses) <- bands
  class(responses) <- append("MarginalResponse", class(responses))

  # Visualize
  if (visualize) {
    plot(responses)
  }

  # Return
  return(responses)
}

# marginal_response end
