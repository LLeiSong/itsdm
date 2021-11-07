#' @title Function to run extended isolation forest as SDM.
#' @description Call isolation forest and its variations to do species distribution modeling.
#' @param occ (data.frame, sf, SpatialPointsDataFrame)
#' The occurrence dataset for training.
#' There must be column x and y for coordinates if it is a data.frame.
#' @param occ_test (data.frame, sf, SpatialPointsDataFrame)
#' The occurrence dataset for independent test.
#' There must be column x and y for coordinates if not NULL and it is a data.frame.
#' If NULL, no independent test will be used. The default is NULL.
#' @param occ_crs (numeric or crs) The EPSG number or crs object of occurrence CRS.
#' The default value is 4326, which is the geographic coordinate system.
#' @param variables (RasterStack or stars) The stack of environmental variables.
#' @param ntrees (integer) The number of trees for the isolation forest.
#' The default is 100.
#' @param sample_size (integer) A number for sampling size in `[2, nrow(occ)]`.
#' The default is NULL.
#' @param sample_rate (numeric) A rate for sampling size in [0, 1].
#' The default is NULL. Only set either `sample_size` or `sample_rate`.
#' @param ndim (integer) ExtensionLevel for isolation forest,
#' it must be no smaller than the dimension of environmental variables.
#' When it is 0, the model is a traditional isolation forest, otherwise the model
#' is an extended isolation forest. The default is 0.
#' @param ... Other arguments that `isotree::isolation.forest` needs.
#' @param response (logical) if TRUE, generate response curves.
#' The default is TRUE.
#' @param check_variables (logical) if TRUE, check the variable importance.
  #' The default is TRUE.
#' @param visualize (logical) if TRUE, generate the essential figures related to
#' the model. The default is FALSE.
#' @return (POIsotree) a list of model, occ var_train, pred_train,
#' pred_test,  prediction, responses, eval_train, and eval_test.
#' @import checkmate
#' @importFrom dplyr select slice mutate sample_n
#' @importFrom stars st_as_stars st_extract st_xy2sfc st_rasterize
#' @importFrom sf st_as_sf st_crs st_transform st_drop_geometry
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
  # Isotree-related inputs
  ntrees = 100L, # number of trees
  sample_size = NA, # sample size
  sample_rate = NA, # sample rate
  ndim = 0L, # extension level
  # Other arguments of isolation.forest
  ...,
  # Other general inputs
  response = TRUE, # if generate response curves
  check_variable = TRUE, # if generate var importance
  visualize = FALSE){ # if plot the curves

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
  checkmate::assert_multi_class(
    variables, c('RasterStack', 'RasterLayer', 'stars'))
  checkmate::assert_int(ntrees)
  checkmate::assert_int(
    sample_size, lower = 2,
    upper = nrow(occ), na.ok = T)
  checkmate::assert_number(
    sample_rate, lower = 0,
    upper = 1, null.ok = T)
  if (is(variables, 'RasterStack')) dim_max <- nlayers(variables)
  if (is(variables, 'RasterLayer')) dim_max <- 1
  if (is(variables, 'stars')) {
    dim_max <- length(st_get_dimension_values(variables, 'band'))}
  checkmate::assert_int(
    ndim, lower = 0,
    upper = dim_max - 1, na.ok = T)
  checkmate::assert_logical(response)
  checkmate::assert_logical(check_variable)
  checkmate::assert_logical(visualize)

  # Check inputs - level 2
  ## Check related columns
  if (is(occ, 'data.frame')) {
    if(!all(c("x", "y") %in% colnames(occ))){
      stop("There must be x and y column in occ.")}}
  if (is(occ_test, 'data.frame')) {
    if(!all(c("x", "y") %in% colnames(occ_test))){
      stop("There must be x and y column in occ_test.")}}
  if (!is.na(sample_size) & !is.na(sample_rate)){
    stop('Only set sample_size or sample_rate.')
  } else if (is.na(sample_size) & is.na(sample_rate)) {
    warning('No sample_size or sample_rate set. Use full sample.')
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
    variable_dependences <- variable_dependence(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry())
  } else {
    marginal_responses <- NULL
    independent_responses <- NULL
    variable_dependences <- NULL}

  if (check_variable) {
    vimp <- variable_analysis(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry(),
      var_occ_test = occ_test_mat %>% st_drop_geometry(),
      variables = variables)
  } else {
    vimp <- NULL
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
  var_occ_bg <- st_extract(x = split(variables, 'band'),
                           at = pts_bg_occ)
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
    var_occ_test_bg <- st_extract(x = split(variables, 'band'),
                                  at = pts_bg_occ_test)
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
              eval_train = eval_train,
              var_test = occ_test_mat,
              pred_test = occ_test_pred,
              eval_test = eval_test,
              prediction = var_pred,
              marginal_responses = marginal_responses,
              independent_responses = independent_responses,
              variable_dependence = variable_dependences,
              variable_importance = vimp)
  class(out) <- append("POIsotree", class(out))

  # Return
  out
}

# isotree_po end
