#' A function to run extended isolation forest as SDM.
#' @description This function allows you to call extended isolation forest to do species distribution modeling.
#' @param occ (data.frame) The occurrence data frame for training.
#' There must be column x and y for coordinates. This argument is required
#' @param occ_test (data.frame) The occurrence data frame for testing.
#' There must be column x and y for coordinates if not NULL.
#' If NULL, no test will be used. The default is NULL.
#' @param eval_method (character) Method for evaluation.
#' Should be one of "po" (present-only) and "pa" (present-absence).
#' The default is po.
#' @param variables (RasterStack) the stack of environmental variables.
#' This argument is required.
#' @param ntrees (integer) the number of trees for the forest.
#' The default is 100.
#' @param sample_size (integer) a number for sampling size.
#' If a sample_rate is set, it will be ignored.
#' The default is NULL.
#' @param sample_rate (numeric) a rate for sampling size. Should be in [0, 1].
#' The default is NULL.
#' @param ExtensionLevel (integer) ExtensionLevel for isolation forest,
#' it must be smaller than the dimension of environmental variables.
#' The default is 0.
#' @param seed (integer) set the seed for random process.
#' The default is 100.
#' @param response (logical) if TRUE, generate response curves. The default is TRUE.
#' @param simple (logical) if TRUE, make the evaluation a simple data.frame;
#' if FALSE, make the evaluation a full list. The default is FALSE.
#' @param plots (logical) if TRUE, plot out the evaluation figures. The default is TRUE.
#' @param verbose (logical) if TRUE, print out the progress. The default is TRUE.
#' @return (eif_sdm) a list of model, occ var_train, pred_train, pred_test,
#' prediction, responses, eval_train, and eval_test.
#' @importFrom isotree isolation.forest
#' @export
#' @examples
#' isotree_po(occ = occ, variables = env_vars, ntrees = 500,
#' sample_size = 400, ndim = 0)
#'
isotree_po <- function(
  # SDM-related inputs
  occ, # sf, sp, or data.frame with x and y
  occ_test = NULL, # sf, sp, or data.frame with x and y
  occ_crs = 4326,
  variables, # rasterstack or stars. if stars, must have dimension band


  occ_type = 'po', # po or pa
  label_column = NULL, # label column, which must be in names(occ) if not NULL
  eval_method = "po", # po or pa or both.
  presence_convertion = 'logistic', # linear, logistic, cloglog


  # Isotree-related inputs
  ntrees = 100L, # number of trees
  sample_size = NA, # sample size
  sample_rate = NA, # sample rate
  ndim = 0L, # extension level
  # Other arguments of isolation.forest
  ...,
  # Other general inputs
  response = TRUE, # if generate response curves
  variable_importance = TRUE, # if generate var importance
  simple = FALSE, # if simplify evaluation results
  visualize = TRUE){ # if plot the curves

  # Check inputs - level 1
  checkmate::assert_multi_class(
    occ, c('data.frame', 'sf',
           'SpatialPointsDataFrame',
           null.ok = F))
  checkmate::assert_multi_class(
    occ_test, c('data.frame', 'sf',
                'SpatialPointsDataFrame',
                null.ok = T))
  checkmate::assert_multi_class(
    occ_crs, c('numeric', 'crs',
                null.ok = T))
  checkmate::assert_choice(occ_type, c('po', 'pa'))
  checkmate::assert_character(label_column, null.ok = T)
  checkmate::assert_multi_class(
    variables, c('RasterStack', 'RasterLayer', 'stars'))
  checkmate::assert_choice(eval_method, c('po', 'pa', 'both'))
  checkmate::assert_int(ntrees)
  checkmate::assert_integer(
    sample_size, lower = 0,
    upper = nrow(occ), na.ok = T)
  checkmate::assert_number(
    sample_size, lower = 0,
    upper = 1, null.ok = T)
  if (is(variables, 'RasterStack')) dim_max <- nlayers(variables)
  if (is(variables, 'RasterLayer')) dim_max <- 1
  if (is(variables, 'stars')) {
    dim_max <- length(st_get_dimension_values(var_stars, 'band'))}
  checkmate::assert_integer(
    ndim, lower = 0,
    upper = dim_max - 1, na.ok = T)
  checkmate::assert_int(seed)
  checkmate::assert_logical(response)
  checkmate::assert_logical(variable_importance)
  checkmate::assert_logical(simple)
  checkmate::assert_logical(visualize)

  # Check inputs - level 2
  ## Check related columns
  if (is(occ, 'data.frame')) {
    if(!all(c("x", "y") %in% colnames(occ))){
      stop("There must be x and y column in occ.")}}
  if (is(occ_test, 'data.frame')) {
    if(!all(c("x", "y") %in% colnames(occ_test))){
      stop("There must be x and y column in occ_test.")}}
  if (!is.null(label_column) &
      !(label_column %in% colnames(occ))){
    stop('No label_column found in occ.')}
  if (!is.null(label_column) &
      !(label_column %in% colnames(occ_test))){
    stop('No label_column found in occ_test.')}
  if (!is.na(sample_size) & !is.na(sample_rate)){
    stop('Only set sample_size or sample_rate.')
  } else if (is.na(sample_size) & is.na(sample_rate)) {
    message('No sample_size or sample_rate set. Use full sample.')
    sample_rate <- 1
  }

  # Reformat the inputs
  # Variables, use stars
  if (is(variables, 'RasterStack') | is(variables, 'RasterLayer')){
    variables <- st_as_stars(variabels)
  }
  # Occurrence
  if (is(occ, 'data.frame')) {
    pts_occ <- occ %>%
      st_as_sf(coords = c('x', 'y'), crs = occ_crs)
  } else {
    pts_occ <- st_as_sf(occ)
  }
  if (st_crs(variables) != st_crs(pts_occ)){
    pts_occ <- st_transform(pts_occ, st_crs(variables))
  }
  if(!is.null(occ_test)) {
    if (is(occ_test, 'data.frame')) {
      pts_occ_test <- occ_test %>%
        st_as_sf(coords = c('x', 'y'), crs = occ_crs)
    } else{
      pts_occ_test <- st_as_sf(occ_test)
    }
    if (st_crs(variables) != st_crs(pts_occ_test)){
      pts_occ_test <- st_transform(pts_occ_test, st_crs(variables))
    }
  }

  # Sample size
  if (is.na(sample_size)) sample_size <- nrow(occ) * sample_rate

  # Extract values
  occ_mat <- st_extract(x = variables, at = pts_occ) %>%
    st_as_sf()
  if (!is.null(occ_test)) {
    occ_test_mat <- st_extract(x = variables, at = pts_occ_test) %>%
      st_as_sf()
  } else {
    occ_test_mat <- NULL
  }

  # Train the model
  isotree_mod <- isolation.forest(
    occ_mat %>% st_drop_geometry(),
    ntrees = ntrees,
    sample_size = sample_size,
    ndim = ndim,
    ...)

  # Do prediction
  ## Raster
  var_pred <- predict(split(variables, 'band'), isotree_mod)
  # Stretch result to be comparable with other results
  var_pred <- 1 - var_pred
  var_pred <- .stars_stretch(var_pred)

  ## Training
  occ_pred <- st_extract(x = var_pred, at = pts_occ)

  ## Test
  if (!is.null(occ_test)){
    occ_test_pred <- st_extract(x = var_pred, at = pts_occ_test)
  } else {
    occ_test_pred <- NULL
  }

  # Generate response curves
  if (response) {
    marginal_responses <- marginal_response(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry(),
      variables = variables)
    independent_responses <- independent_response(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry(),
      variables = variables)
  } else {responses <- NULL}

  if (variable_importance) {
    vimp <- variable_importance(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry(),
      var_occ_test = occ_test_mat %>% st_drop_geometry())
  }

  # Make evaluation using presence-only
  ## Sample points from background
  set.seed(123)
  stars_mask <- var_pred
  stars_mask[[1]][!is.na(stars_mask[[1]])] <- 1
  pts_bg_occ <- pts_occ %>% mutate(value = 0) %>%
    st_rasterize(stars_mask)
  pts_bg_occ[[1]][pts_bg_occ[[1]] == 0] <- NA
  pts_bg_occ <- pts_bg_occ %>%
    st_xy2sfc(as_points = T) %>% st_as_sf() %>%
    sample_n(nrow(pts_occ)) %>% select(geometry)
  occ_bg_pred <- st_extract(x = var_pred, at = pts_bg_occ)
  rm(stars_mask, pts_bg_occ)

  ### Do the same to test if has any
  if(!is.null(occ_test)) {
    set.seed(123)
    stars_mask <- var_pred
    stars_mask[[1]][!is.na(stars_mask[[1]])] <- 1
    pts_bg_occ_test <- pts_occ_test %>% mutate(value = 0) %>%
      st_rasterize(stars_mask)
    pts_bg_occ_test[[1]][pts_bg_occ_test[[1]] == 0] <- NA
    pts_bg_occ_test <- pts_bg_occ_test %>%
      st_xy2sfc(as_points = T) %>% st_as_sf() %>%
      sample_n(nrow(pts_occ_test)) %>% select(geometry)
    occ_test_bg_pred <- st_extract(x = var_pred, at = pts_bg_occ_test)
    rm(stars_mask, pts_bg_occ_test)
  }

  ## Calculate
  eval_train <- evaluate_po(isotree_mod,
                            occ_pred$prediction,
                            occ_bg_pred$prediction,
                            na.omit(as.vector(var_pred[[1]])))
  if(!is.null(occ_test)){
    eval_test <- evaluate_po(isotree_mod,
                             occ_test_pred$prediction,
                             occ_test_bg_pred$prediction,
                             na.omit(as.vector(var_pred[[1]])))
  } else eval_test <- NULL

  # Return
  out <- list(model = isotree_mod,
              pts_occ = pts_occ,
              pts_occ_test = pts_occ_test,
              var_train = occ_mat,
              pred_train = occ_pred,
              eval_train = occ_eval,
              pred_test = occ_test_pred,
              eval_test = occ_test_eval,
              prediction = var_pred,
              marginal_responses = marginal_responses,
              independent_responses = independent_responses)
  class(out) <- append("isotree_presence_only", class(out))

  # Return
  out
}

# isotree_po end
