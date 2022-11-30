#' @title Calculate shapley values-based spatial response.
#' @description Calculate spatially SHAP-based response figures.
#' They can help to diagnose both how and where the species
#' responses to environmental variables.
#' @param model (`isolation_forest` or other model). It could
#' be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' It also could be other user-fitted models as long as the `pfun` can work on it.
#' @param var_occ (`data.frame`, `tibble`) The `data.frame` style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (`stars`) The `stars` of environmental variables.
#' It should have multiple `attributes` instead of `dims`.
#' If you have `raster` object instead, you
#' could use \code{\link{st_as_stars}} to convert it to `stars` or use
#' \code{\link{read_stars}} directly read source data as a `stars`.
#' You also could use item `variables` of `POIsotree` made by function
#' \code{\link{isotree_po}}.
#' @param target_vars (a `vector` of `character`) The selected variables to
#' process. If it is `NULL`, all variables will be used.
#' @param shap_nsim (`integer`) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value. See details in documentation
#' of function \code{\link{explain}} in package `fastshap`.
#' When the number of variables is large, a smaller shap_nsim could be used.
#' Be cautious that making SHAP-based spatial dependence will be slow
#' because of Monte-Carlo computation for all pixels.
#' But it is worth the time because it is much more
#' informative. See details in documentation of function \code{\link{explain}}
#' in package `fastshap`. The default is 10. Usually a value 10 - 20 is enough.
#' @param seed (`integer`) The seed for any random progress. The default is `10L`.
#' @param pfun (`function`) The predict function that requires two arguments,
#' `object` and `newdata`.
#' It is only required when `model` is not \code{isolation_forest}.
#' The default is the wrapper function designed for iForest model in `itsdm`.
#'
#' @return (`SHAPSpatial`) A list of
#' \itemize{
#' A list of `stars` object of spatially SHAP-based response of all variables
#' }
#'
#' @details
#' The values show how each environmental variable affects the modeling
#' prediction in space. These maps could help to answer questions of where in
#' terms of environmental response.
#'
#' @seealso
#' \code{\link{spatial_response}}
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
#' # Prepare data
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
#'   slice('band', c(1, 5, 12))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' shap_spatial <- shap_spatial_response(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables,
#'   shap_nsim = 1)
#'
#' shap_spatial <- shap_spatial_response(
#'  model = mod$model,
#'  target_vars = c("bio1", "bio12"),
#'  var_occ = mod$vars_train,
#'  variables = mod$variables,
#'  shap_nsim = 1)
#'
#' \dontrun{
#' ##### Use Random Forest model as an external model ########
#' library(randomForest)
#'
#' # Prepare data
#' data("occ_virtual_species")
#' obs_df <- occ_virtual_species %>%
#'   filter(usage == "train")
#'
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12)) %>%
#'   split()
#'
#' model_data <- stars::st_extract(
#'   env_vars, at = as.matrix(obs_df %>% select(x, y))) %>%
#'   as.data.frame()
#' names(model_data) <- names(env_vars)
#' model_data <- model_data %>%
#'   mutate(occ = obs_df[['observation']])
#' model_data$occ <- as.factor(model_data$occ)
#'
#' mod_rf <- randomForest(
#'   occ ~ .,
#'   data = model_data,
#'   ntree = 200)
#'
#' pfun <- function(X.model, newdata) {
#'   # for data.frame
#'   predict(X.model, newdata, type = "prob")[, "1"]
#' }
#'
#' shap_spatial <- shap_spatial_response(
#'   model = mod_rf,
#'   target_vars = c("bio1", "bio12"),
#'   var_occ = model_data %>% select(-occ),
#'   variables = env_vars,
#'   shap_nsim = 10,
#'   pfun = pfun)
#' }
#'
shap_spatial_response <- function(model,
                                  var_occ,
                                  variables,
                                  target_vars = NULL,
                                  shap_nsim = 10,
                                  seed = 10,
                                  pfun = .pfun_shap){
  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_class(variables, 'stars')
  stopifnot(length(dim(variables)) <= 2)
  checkmate::assert_int(shap_nsim)
  checkmate::assert_int(seed)
  checkmate::assert_character(target_vars, null.ok = TRUE)
  bands <- names(variables)
  stopifnot(all(bands %in% colnames(var_occ)))
  stopifnot(all(target_vars %in% bands))

  # Check and reformat inputs
  if (is.null(target_vars)) target_vars <- bands

  ## Split numeric and categorical
  bands <- names(variables)
  isfacor <- as.vector(sapply(variables, is.factor))
  bands_cont <- bands[!isfacor]
  bands_cat <- bands[isfacor]

  # Do SHAP
  set.seed(seed)
  x_shap <- variables %>% as.data.frame()
  val_ids <- apply(x_shap[-c(1,2)], 1, function(x) all(!is.na(x)))

  # Calculate Shapley values for full records
  shap_explain <- explain(
    model,
    X = var_occ,
    newdata = x_shap %>% filter(val_ids) %>% select(all_of(bands)),
    nsim = shap_nsim,
    pred_wrapper = pfun)

  # Mosaic back the background pixels
  x_shap[] <- sapply(x_shap, as.numeric)
  x_shap[!val_ids, bands] <- NA
  x_shap[val_ids, ] <- cbind(
    x_shap %>% filter(val_ids) %>% select(.data$x, .data$y),
    shap_explain)

  shap_spatial <- lapply(target_vars, function(nm) {
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
  names(shap_spatial) <- target_vars
  class(shap_spatial) <- append("SHAPSpatial", class(shap_spatial))

  # return
  shap_spatial
}
# shap_spatial_response end
