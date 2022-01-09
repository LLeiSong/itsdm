#' @title Function to run extended isolation forest as SDM.
#' @description Call isolation forest and its variations to do
#' species distribution modeling and optionally do model explanation.
#' @param occ (`data.frame`, `sf`, `SpatialPointsDataFrame`)
#' The occurrence dataset for training.
#' There must be column `x` and `y` for coordinates if
#' it is a regular `data.frame`.
#' @param occ_test (`data.frame`, `sf`, `SpatialPointsDataFrame`, or `NULL`)
#' The occurrence dataset for independent test. The same structure as `occ`.
#' If not `NULL`, there must be column `x` and `y` for coordinates when it is a
#' regular `data.frame`. If `NULL`, no independent test will be used.
#' The default is `NULL`.
#' @param occ_crs (`numeric` or \code{\link{crs}}) The EPSG number or
#' \code{\link{crs}} object of occurrence CRS.
#' The default value is `4326`, which is the geographic coordinate system.
#' @param variables (`RasterStack` or `stars`) The stack of environmental variables.
#' @param categ_vars (`vector` of `character` or `NULL`) The names of categorical
#' variables. Must be the same as the names in `variables`.
#' @param ntrees (`integer`) The number of trees for the isolation forest. It must
#' be integer, which you could use function \code{\link{as.integer}} to convert to.
#' The default is `100L`.
#' @param sample_size (`integer` or `NULL`) Alternative argument for `sample_rate`.
#' If not `NULL`, it should be a number for sampling size in `[2, nrow(occ)]`.
#' It must be integer, which you could use function \code{\link{as.integer}} to
#' convert to. The default is `NULL`.
#' Only set either `sample_size` or `sample_rate`.
#' @param sample_rate (`numeric` or `NULL`) Alternative argument for `sample_size`.
#' If not `NULL`, it should be a rate for sampling size in `[0, 1]`.
#' The default is `NULL`. Only set either `sample_size` or `sample_rate`.
#' @param ndim (`integer`) ExtensionLevel for isolation forest. It must
#' be integer, which you could use function \code{\link{as.integer}} to convert
#' to. Also, it must be no smaller than the dimension of environmental variables.
#' When it is 1, the model is a traditional isolation forest, otherwise the model
#' is an extended isolation forest. The default is 1.
#' @param seed (`integer`) The random seed used in the modeling. It should be an
#' integer. The default is `10L`.
#' @param ... Other arguments that \code{\link{isolation.forest}} needs.
#' @param offset (`numeric`) The offset to adjust fitted suitability. The default
#' is zero. Highly recommend to leave it as default.
#' @param response (`logical`) If `TRUE`, generate response curves.
#' The default is `TRUE`.
#' @param spatial_response (`logical`) If `TRUE`, generate spatial response maps.
#' The default is `TRUE` because it might be slow. NOTE that here SHAP-based map
#' is not generated because it is slow. If you want it be mapped, you could call
#' function \code{\link{spatial_response}} to make it.
#' @param check_variable (`logical`) If `TRUE`, check the variable importance.
#' The default is `TRUE`.
#' @param visualize (`logical`) If `TRUE`, generate the essential figures
#' related to the model. The default is `FALSE`.
#'
#' @return (`POIsotree`) A list of
#' \itemize{
#' \item{model (\code{\link{isolation.forest}}) The threshold set in
#' function inputs}
#' \item{variables (`stars`) The formatted image stack of
#' environmental variables}
#' \item{pts_occ (\code{\link{sf}}) A \code{\link{sf}} of training occurrence
#' dataset}
#' \item{pts_bg_occ (\code{\link{sf}}) A \code{\link{sf}} of background points
#' for training dataset evaluation or SHAP dependence plot}
#' \item{pts_occ_test (\code{\link{sf}} or `NULL`) A \code{\link{sf}} of test
#' occurrence dataset}
#' \item{pts_bg_occ_test (\code{\link{sf}} or `NULL`) A \code{\link{sf}} of
#' background points for test dataset evaluation or SHAP dependence plot}
#' \item{var_train (\code{\link{sf}}) A \code{\link{sf}} with values of each
#' environmental variables for training occurrence}
#' \item{pred_train (\code{\link{sf}}) A \code{\link{sf}} with values of
#' prediction for training occurrence}
#' \item{eval_train (`POEvaluation`) A list of presence-only evaluation metrics
#' based on training dataset. See details of `POEvaluation` in
#' \code{\link{evaluate_po}}}
#' \item{var_test (\code{\link{sf}} or `NULL`) A \code{\link{sf}} with values of each
#' environmental variables for test occurrence}
#' \item{pred_test (\code{\link{sf}} or `NULL`) A \code{\link{sf}} with values of
#' prediction for test occurrence}
#' \item{eval_test (`POEvaluation` or `NULL`) A list of presence-only evaluation metrics
#' based on test dataset.
#' See details of `POEvaluation` in \code{\link{evaluate_po}}}
#' \item{prediction (`stars`) The predicted environmental suitability}
#' \item{marginal_responses (`MarginalResponse` or `NULL`) A list of marginal response
#' values of each environmental variables.
#' See details in \code{\link{marginal_response}}}
#' \item{offset (`numeric`) The offset value set as inputs.}
#' \item{independent_responses (`IndependentResponse` or `NULL`) A list of independent
#' response values of each environmental variables.
#' See details in \code{\link{independent_response}}}
#' \item{shap_dependences (`ShapDependence` or `NULL`) A list of variable
#' dependence values of each environmental variables.
#' See details in \code{\link{shap_dependence}}}
#' \item{spatial_responses (`SpatialResponse` or `NULL`) A list of spatial variable
#' dependence values of each environmental variables.
#' See details in \code{\link{shap_dependence}}}
#' \item{variable_analysis (`VariableAnalysis` or `NULL`) A list of variable importance
#' analysis based on multiple metrics.
#' See details in \code{\link{variable_analysis}}}}
#'
#' @seealso
#' \code{\link{evaluate_po}}, \code{\link{marginal_response}},
#' \code{\link{independent_response}}, \code{\link{shap_dependence}},
#' \code{\link{spatial_response}}, \code{\link{variable_analysis}},
#' \code{\link{isolation.forest}}
#'
#' @references
#' \itemize{
#' \item{Liu, Fei
#' Tony, Kai Ming Ting, and Zhi-Hua Zhou. "Isolation forest."
#' \emph{2008 eighth ieee international conference on data mining}.IEEE, 2008.
#' \doi{10.1109/ICDM.2008.17}}
#' \item{Liu, Fei Tony, Kai Ming
#' Ting, and Zhi-Hua Zhou. "Isolation-based anomaly detection."
#' \emph{ACM Transactions on Knowledge Discovery from Data (TKDD)} 6.1 (2012): 1-39.
#' \doi{10.1145/2133360.2133363}}
#' \item{Liu, Fei Tony,
#' Kai Ming Ting, and Zhi-Hua Zhou. "On detecting clustered anomalies using
#' SCiForest." \emph{Joint European Conference on Machine Learning and
#' Knowledge Discovery in Databases}. Springer, Berlin, Heidelberg, 2010.
#' \doi{10.1007/978-3-642-15883-4_18}}
#' \item{Ha
#' riri, Sahand, Matias Carrasco Kind, and Robert J. Brunner. "Extended
#' isolation forest." \emph{IEEE Transactions on Knowledge and Data Engineering (2019)}.
#' \doi{10.1109/TKDE.2019.2947676}}
#' \item{\url{https://github.com/david-cortes/isotree}}
#' \item{References of related feature such as response curves and variable importance
#' will be listed under their own functions}
#' }
#'
#' @details
#' Please read details of algorithm \code{\link{isolation.forest}} on
#' \url{https://github.com/david-cortes/isotree}, and
#' the R documentation of function \code{\link{isolation.forest}}.
#'
#' @import checkmate
#' @importFrom raster nlayers
#' @importFrom dplyr select slice mutate sample_n
#' @importFrom stars st_as_stars st_extract st_xy2sfc st_rasterize
#' @importFrom sf st_as_sf st_crs st_transform st_drop_geometry
#' @importFrom isotree isolation.forest
#' @importFrom rlang := .data
#' @importFrom stats na.omit predict
#' @importFrom methods is
#' @export
#' @examples
#' \donttest{
#' # Using a pseudo presence-only occurrence dataset of
#' # virtual species provided in this package
#' library(dplyr)
#' library(sf)
#' library(stars)
#' library(itsdm)
#'
#' # Prepare data
#' data("occ_virtual_species")
#' occ_virtual_species <- occ_virtual_species %>%
#'   mutate(id = row_number())
#'
#' set.seed(11)
#' occ <- occ_virtual_species %>% sample_frac(0.7)
#' occ_test <- occ_virtual_species %>% filter(! id %in% occ$id)
#' occ <- occ %>% select(-id)
#' occ_test <- occ_test %>% select(-id)
#'
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12))
#'
#' # Modeling
#' mod_virtual_species <- isotree_po(
#'   occ = occ, occ_test = occ_test,
#'   variables = env_vars, ntrees = 10,
#'   sample_rate = 0.6, ndim = 1L,
#'   seed = 123L)
#'
#' # Check results
#' ## Evaluation based on training dataset
#' print(mod_virtual_species$eval_train)
#' plot(mod_virtual_species$eval_train)
#'
#' ## Response curves
#' plot(mod_virtual_species$marginal_responses)
#' plot(mod_virtual_species$independent_responses,
#'   target_var = c('bio1', 'bio5'))
#' plot(mod_virtual_species$shap_dependence)
#'
#' ## Relationships between target var and related var
#' plot(mod_virtual_species$shap_dependence,
#'   target_var = c('bio1', 'bio5'),
#'   related_var = 'bio12', smooth_span = 0)
#'
#' # Variable importance
#' mod_virtual_species$variable_analysis
#' plot(mod_virtual_species$variable_analysis)
#'}
#'
isotree_po <- function(
  # SDM-related inputs
  occ, # sf, sp, or data.frame with x and y
  occ_test = NULL, # sf, sp, or data.frame with x and y
  occ_crs = 4326,
  variables, # rasterstack or stars. if stars, must have dimension band
  categ_vars = NULL, # categorical variables
  # Isotree-related inputs
  ntrees = 100L, # number of trees
  sample_size = NA, # sample size
  sample_rate = 1.0, # sample rate
  ndim = 1L, # extension level
  seed = 10L, # random seed, must be integer
  # Other arguments of isolation.forest
  ...,
  # Probability conversion
  offset = 0.0,
  # Other general inputs
  response = TRUE, # if generate response curves
  spatial_response = TRUE, # if generate spatial response maps
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
    variables, c('RasterStack', 'stars'))
  checkmate::assert_vector(categ_vars, null.ok = T)
  checkmate::assert_int(ntrees)
  checkmate::assert_int(
    sample_size, lower = 2,
    upper = nrow(occ), na.ok = T)
  checkmate::assert_number(
    sample_rate, lower = 0,
    upper = 1, na.ok = T)
  if (is(variables, 'RasterStack')) dim_max <- nlayers(variables)
  if (is(variables, 'stars')) {
    if (length(dim(variables)) == 2){
      dim_max <- length(variables)
    } else if (length(dim(variables)) == 3) {
      if (!is.null(categ_vars)){
        warning(paste0('Categorical layers are merged to an dimension.',
                       ' Be careful if they are original values.'))
      }
      dim_max <- length(st_get_dimension_values(variables, 3))
    } else {
      stop('variables has more than 3 dimensions, do not know which one to use.')
    }
  }
  checkmate::assert_int(
    ndim, lower = 1,
    upper = dim_max, na.ok = T)
  checkmate::assert_logical(response)
  checkmate::assert_logical(spatial_response)
  checkmate::assert_logical(check_variable)
  checkmate::assert_logical(visualize)

  # Check categ_cols
  if (exists('categ_cols')) {
    stop('Set categ_vars instead.')
  }

  # Check inputs - level 2
  ## Check related columns
  if (is(occ, 'data.frame') & (!is(occ, 'sf')) & (!is(occ, 'Spatial'))) {
    if(!all(c("x", "y") %in% colnames(occ))){
      stop("There must be x and y column in occ.")}}
  if (is(occ_test, 'data.frame') & (!is(occ_test, 'sf')) &
      (!is(occ_test, 'Spatial'))) {
    if(!all(c("x", "y") %in% colnames(occ_test))){
      stop("There must be x and y column in occ_test.")}}
  if (!is.na(sample_size) & !is.na(sample_rate)){
    stop('Only set sample_size or sample_rate.')
  } else if (is.na(sample_size) & is.na(sample_rate)) {
    warning('No sample_size or sample_rate set. Use full sample.')
    sample_rate <- 1.0
  }

  # Reformat the inputs
  # Variables
  ## -- RasterStack - resume cat original values, convert to stars with multiple
  ##                  attributes, check if categ_vars are existing cat vars.
  ## -- Stars with 3 dims - split it, no need to check cat vars
  ## -- Stars with multiple attributes - check categ_vars are existing cat vars
  if (is(variables, 'RasterStack')){
    # Check categ_vars
    ## step 1
    stopifnot(all(categ_vars %in% names(variables)))
    ## step 2
    if (!identical(names(variables) %in% categ_vars,
                   is.factor(variables))) {
      warning('Categorical layers detected in RasterStack do not match with categ_vars.')
    }
    variables <- .remove_cats(variables)
    variables <- st_as_stars(variables) %>% split('band')
  } else {
    if (length(dim(variables)) == 3) {
      variables <- variables %>% split(3)
      # Check categ_vars
      stopifnot(all(categ_vars %in% names(variables)))
    } else {
      # Check categ_vars
      ## step 1
      stopifnot(all(categ_vars %in% names(variables)))
      # step 2
      isfacor <- as.vector(sapply(variables, is.factor))
      if (!identical(names(variables) %in% categ_vars,
                     isfacor)) {
        warning('Categorical layers detected in RasterStack do not match with categ_vars.')
      }
    }
  }

  # Convert to factors
  for (nm in categ_vars) {
    if (!is.factor(variables[[nm]])) {
      variables <- variables %>%
        mutate(!!nm := as.factor(variables[[nm]]))
    }
  }

  # Occurrence
  if (is(occ, 'data.frame') & (!is(occ, 'sf')) &
      (!is(occ, 'Spatial'))) {
    pts_occ <- occ %>%
      st_as_sf(coords = c('x', 'y'), crs = occ_crs)
  } else {
    pts_occ <- st_as_sf(occ)
  }
  if (st_crs(variables) != st_crs(pts_occ)){
    pts_occ <- st_transform(pts_occ, st_crs(variables))
  }
  if(!is.null(occ_test)) {
    if (is(occ_test, 'data.frame') & (!is(occ_test, 'sf')) &
        (!is(occ_test, 'Spatial'))) {
      pts_occ_test <- occ_test %>%
        st_as_sf(coords = c('x', 'y'), crs = occ_crs)
    } else{
      pts_occ_test <- st_as_sf(occ_test)
    }
    if (st_crs(variables) != st_crs(pts_occ_test)){
      pts_occ_test <- st_transform(pts_occ_test, st_crs(variables))
    }
  } else pts_occ_test <- NULL

  # Extract values
  occ_mat <- st_extract(x = variables, at = pts_occ) %>%
    st_as_sf() %>% na.omit()
  # Remove NAs from pts as well
  pts_occ <- occ_mat %>% select(.data$geometry)

  # Do the same thing to test dataset if there is one
  if (!is.null(occ_test)) {
    occ_test_mat <- st_extract(x = variables, at = pts_occ_test) %>%
      st_as_sf() %>% na.omit()
    pts_occ_test <- pts_occ_test %>% select(.data$geometry)
  } else {
    occ_test_mat <- NULL}

  # Sample size
  if (is.na(sample_size)) sample_size <- round(nrow(occ_mat) * sample_rate)

  # Train the model
  isotree_mod <- isolation.forest(
    occ_mat %>% st_drop_geometry(),
    ntrees = ntrees,
    sample_size = sample_size,
    ndim = ndim,
    seed = seed,
    ...)

  # Do prediction
  ## Raster
  var_pred <- probability(isotree_mod, variables, offset = offset)

  ## Training
  ### There should not are any NAs coming from variables
  occ_pred <- st_extract(x = var_pred, at = pts_occ)

  ## Test
  if (!is.null(occ_test)){
    occ_test_pred <- st_extract(x = var_pred, at = pts_occ_test)
  } else {
    occ_test_pred <- NULL}

  # Sample points from background
  set.seed(seed)
  stars_mask <- var_pred
  stars_mask[[1]][!is.na(stars_mask[[1]])] <- 1
  pts_bg_occ <- pts_occ %>% mutate(value = 0) %>%
    st_rasterize(stars_mask)
  pts_bg_occ[[1]][pts_bg_occ[[1]] == 0] <- NA
  pts_bg_occ <- pts_bg_occ %>%
    st_xy2sfc(as_points = T) %>% st_as_sf() %>%
    sample_n(nrow(pts_occ)) %>% select(.data$geometry)
  occ_bg_pred <- st_extract(x = var_pred, at = pts_bg_occ)
  var_occ_bg <- st_extract(x = variables,
                           at = pts_bg_occ)
  rm(stars_mask)

  ## Do the same to test if has any
  if(!is.null(occ_test)) {
    set.seed(seed)
    stars_mask <- var_pred
    stars_mask[[1]][!is.na(stars_mask[[1]])] <- 1
    pts_bg_occ_test <- pts_occ_test %>% mutate(value = 0) %>%
      st_rasterize(stars_mask)
    pts_bg_occ_test[[1]][pts_bg_occ_test[[1]] == 0] <- NA
    pts_bg_occ_test <- pts_bg_occ_test %>%
      st_xy2sfc(as_points = T) %>% st_as_sf() %>%
      sample_n(nrow(pts_occ_test)) %>% select(.data$geometry)
    occ_test_bg_pred <- st_extract(x = var_pred, at = pts_bg_occ_test)
    rm(stars_mask)
  } else {
    pts_bg_occ_test <- NULL
    occ_test_bg_pred <- NULL}

  # Generate response curves
  if (response) {
    # marginal variable dependence
    marginal_responses <- marginal_response(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry(),
      variables = variables,
      visualize = visualize)

    # independent variable dependence
    independent_responses <- independent_response(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry(),
      variables = variables,
      visualize = visualize)

    # Shapley value based variable dependence
    shap_var_occ <- rbind(
      occ_mat %>% st_drop_geometry(),
      var_occ_bg %>% st_drop_geometry())
    shap_dependences <- shap_dependence(
      model = isotree_mod,
      var_occ = shap_var_occ,
      visualize = visualize,
      seed = seed)
    rm(var_occ_bg, shap_var_occ)
  } else {
    marginal_responses <- NULL
    independent_responses <- NULL
    shap_dependences <- NULL}

  # Spatial dependence
  if (spatial_response) {
    spatial_responses <- spatial_response(
      model = isotree_mod,
      var_occ = occ_mat %>% st_drop_geometry(),
      variables = variables,
      shap_nsim = 0,
      seed = seed,
      visualize = visualize)
  } else {
    spatial_responses <- NULL}

  # Check variable importance
  if (check_variable) {
    # If there is just training
    if (is.null(pts_occ_test)) {
      vimp <- variable_analysis(
        model = isotree_mod,
        pts_occ = pts_occ,
        pts_occ_test = NULL,
        variables = variables,
        visualize = visualize,
        seed = seed)
    } else {
      # If there are both training and test
      vimp <- variable_analysis(
        model = isotree_mod,
        pts_occ = pts_occ,
        pts_occ_test = pts_occ_test,
        variables = variables,
        visualize = visualize,
        seed = seed)}

  } else {vimp <- NULL}

  # Make evaluation using presence-only
  ## Calculate
  eval_train <- suppressMessages(evaluate_po(isotree_mod,
                            occ_pred$prediction,
                            occ_bg_pred$prediction,
                            na.omit(as.vector(var_pred[[1]])),
                            visualize = visualize))
  rm(occ_bg_pred)
  if(!is.null(occ_test)){
    eval_test <- suppressMessages(evaluate_po(isotree_mod,
                             occ_test_pred$prediction,
                             occ_test_bg_pred$prediction,
                             na.omit(as.vector(var_pred[[1]])),
                             visualize = visualize))
    rm(occ_test_bg_pred)
  } else eval_test <- NULL

  # Return
  out <- list(
    model = isotree_mod, # model
    variables = variables, # the formatted stars
    pts_occ = pts_occ, # clean occurrence points
    pts_bg_occ = pts_bg_occ, # used background samples paired with training
    pts_occ_test = pts_occ_test, # clean test points
    pts_bg_occ_test = pts_bg_occ_test, # used background samples paired with test
    var_train = occ_mat, # table of training features
    pred_train = occ_pred, # prediction of training dataset
    eval_train = eval_train, # evaluation on training dataset
    var_test = occ_test_mat, # table of test dataset
    pred_test = occ_test_pred, # prediction of test dataset
    eval_test = eval_test, # evaluation on test dataset
    prediction = var_pred, # prediction map
    offset = offset, # the offset to adjust result
    marginal_responses = marginal_responses, # marginal response plot
    independent_responses = independent_responses, # independent response plot
    shap_dependences = shap_dependences, # Shapley value based response plot
    spatial_responses = spatial_responses, # spatial response maps
    variable_analysis = vimp) # variable evaluation
  class(out) <- append("POIsotree", class(out))

  # Return
  out
}

# isotree_po end
