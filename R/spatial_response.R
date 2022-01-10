#' @title Calculate spatial response or dependence plot.
#' @description Calculate spatially marginal, independence, and SHAP-based
#' response plots.
#' @param model (`isolation_forest`). It could
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
#' @param shap_nsim (`integer`) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value. See details in documentation
#' of function \code{\link{explain}} in package `fastshap`. Set it to 0 if you
#' don't want to make SHAP-based spatial dependence. When the number of variables
#' is large, a smaller shap_nsim could be used. Be cautious that making
#' SHAP-based spatial dependence will be slow because of Monte-Carlo
#' computation for all pixels. But it is worth the time because it is much more
#' informative. See details in documentation of function \code{\link{explain}}
#' in package `fastshap`. The default is 0. Usually a value 10 - 20 is enough.
#' @param seed (`integer`) The seed for any random progress. The default is `10L`.
#' @param visualize (`logical`) if `TRUE`, plot the response curves.
#' The default is `FALSE`.
#'
#' @return (`SpatialResponse`) A list of
#' \itemize{
#' \item{spatial_marginal_response (`list`) A list of `stars` object of spatially
#' marginal response of all variables}
#' \item{spatial_independent_response (`list`) A list of `stars` object of spatially
#' independent response of all variables}
#' \item{spatial_shap_dependence (`list`) A list of `stars` object of spatially
#' SHAP-based response of all variables}
#' }
#'
#' @details
#' The values show how each environmental variable affects the modeling
#' prediction in space. These maps could help to answer questions of where in
#' terms of environmental response. Compared to marginal dependence or
#' independent dependence maps, SHAP-based maps are way more informative because
#' SHAP-based dependence explain the contribution of each variable to final result.
#'
#' @seealso
#' \code{\link{plot.SpatialResponse}}
#'
#' @importFrom isotree isolation.forest
#' @importFrom dplyr select slice mutate as_tibble pull n %>% group_by
#' @importFrom stars st_as_stars
#' @importFrom stats predict setNames
#' @importFrom tidyselect all_of
#' @importFrom fastshap explain
#' @importFrom rlang :=
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
#' mod <- isotree_po(
#'   occ = occ, occ_test = occ_test,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' spatial_responses <- spatial_response(
#'   model = mod$model,
#'   var_occ = mod$var_train %>% st_drop_geometry(),
#'   variables = mod$variables,
#'   shap_nsim = 1)
#' plot(spatial_responses)
#'
spatial_response <- function(model,
                             var_occ,
                             variables,
                             shap_nsim = 0,
                             seed = 10L,
                             visualize = FALSE) {
  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_class(variables, 'stars')
  stopifnot(length(dim(variables)) <= 2)
  checkmate::assert_int(shap_nsim)
  checkmate::assert_int(seed)
  checkmate::assert_logical(visualize)
  bands <- names(variables)
  stopifnot(all(bands %in% colnames(var_occ)))

  # Reformat data
  ## Split numeric and categorical
  isfacor <- as.vector(sapply(variables, is.factor))
  bands_cont <- bands[!isfacor]
  bands_cat <- bands[isfacor]

  # Do full prediction
  ## Raster
  var_pred_full <- predict(variables, model)
  ## Stretch result to be comparable with other results
  var_pred_full <- 1 - var_pred_full

  ############ marginal response
  # Get means and expand
  means <- do.call(c, lapply(bands, function(nm) {
    if (nm %in% bands_cont) {
      val <- mean(var_occ %>% pull(nm), na.rm = T)
      variables[nm] %>% mutate('{nm}' := val)
    } else {
      val <- .mode(var_occ %>% pull(nm))
      variables[nm] %>% mutate('{nm}' := val)
    }}))

  mar_spatial <- lapply(bands, function(nm) {
    vals_tmp <- means
    vals_tmp <- vals_tmp %>% select(-nm)
    vals_tmp <- c(vals_tmp, variables %>% select(nm))
    pred_tmp <- predict(vals_tmp, model)
    pred_tmp <- 1 - pred_tmp
    .stars_stretch(var_pred_full, new_values = pred_tmp)
  })
  names(mar_spatial) <- bands

  ############ independent response
  models <- lapply(bands, function(nm) {
    ind_var_occ <- var_occ %>% select(nm)
    # Remove and modify some arguments not works
    # for single-variable model, other arguments
    # inherit from model object.
    args_iforest <- c(
      list(data=ind_var_occ, seed=model$random_seed, nthreads=model$nthreads),
      model$params
    )
    args_iforest$ndim <- 1
    args_iforest$ncols_per_tree <- min(ncol(args_iforest$data), args_iforest$ncols_per_tree)
    if (args_iforest$new_categ_action == "impute") {
      args_iforest$new_categ_action <- "weighted"
    }
    do.call(isolation.forest, args_iforest)}) %>%
    setNames(bands)

  ind_spatial <- lapply(bands, function(nm) {
    vals_tmp <- variables %>% select(nm)
    pred_tmp <- predict(vals_tmp, models[[nm]])
    pred_tmp <- 1 - pred_tmp
    .stars_stretch(var_pred_full, new_values = pred_tmp)
  })
  names(ind_spatial) <- bands

  ############ SHAP variable dependence
  if (shap_nsim > 0){
    # Do SHAP
    set.seed(seed)
    x_shap <- variables %>% as.data.frame()
    val_ids <- apply(x_shap[-c(1,2)], 1, function(x) all(!is.na(x)))

    # Calculate Shapley values for full records
    shap_explain <- explain(
      model,
      X = x_shap %>% filter(val_ids) %>% select(bands),
      nsim = shap_nsim,
      pred_wrapper = .pfun_shap)

    # Mosaic back the background pixels
    x_shap[] <- sapply(x_shap, as.numeric)
    x_shap[!val_ids, bands] <- NA
    x_shap[val_ids, ] <- cbind(
      x_shap %>% filter(val_ids) %>% select(.data$x, .data$y),
      shap_explain)
    shap_spatial <- lapply(bands, function(nm) {
      # Continuous variables
      if (nm %in% bands_cont) {
        # Burn values
        variables[nm] %>%
          mutate('{nm}' := x_shap %>% pull(nm)) %>%
          select(all_of(nm))
        # Categorical variables
        ## Each class should have the same value, so aggregate to
        ## each class with mean function
        } else if (nm %in% bands_cat) {
        # Get convert table
        convert_tb <- data.frame(
          shap = x_shap %>% pull(nm),
          class = variables[[nm]] %>% as.vector()) %>%
          group_by(class) %>%
          summarise(shap = mean(.data$shap, na.rm = T)) %>%
          mutate(class = as.integer(.data$class))

        # Generate stars template
        rst_raw <- variables[nm]
        rst_raw[[nm]] <- as.integer(levels(rst_raw[[nm]]))[rst_raw[[nm]]]

        # Replace class value with shap value
        for (i in 1:nrow(convert_tb)) {
          rst_raw[rst_raw == convert_tb[[i, 'class']]] <- convert_tb[i, 'shap']}
        rst_raw
        }
      })
    names(shap_spatial) <- bands
  } else {
    shap_spatial <- NULL
  }

  # Reunion
  spatial_deps <- list(spatial_marginal_response = mar_spatial,
                      spatial_independent_response = ind_spatial,
                      spatial_shap_dependence = shap_spatial)
  class(spatial_deps) <- append("SpatialResponse", class(spatial_deps))

  # Visualize
  if (visualize) {
    plot(spatial_deps)
  }

  # Return
  return(spatial_deps)
}

# spatial_response end
