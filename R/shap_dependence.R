#' @title Calculate Shapley value-based variable dependence.
#' @description Calculate how a species responses to environmental variables
#' using Shapley values.
#' @param model (\code{isolation_forest} or other model). The SDM.
#' It could be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' It also could be other user-fitted models as long as the `pfun` can work on it.
#' @param var_occ (`data.frame`, `tibble`) The `data.frame` style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (`stars`) The `stars` of environmental variables.
#' It should have multiple `attributes` instead of `dims`.
#' If you have `raster` object instead, you
#' could use \code{\link[stars:st_as_stars]{st_as_stars}} to convert it to `stars` or use
#' \code{\link[stars:read_stars]{read_stars}} directly read source data as a `stars`.
#' You also could use item `variables` of `POIsotree` made by function
#' \code{\link{isotree_po}}.
#' @param si (`integer`) The number of samples to generate response curves.
#' If it is too small, the response curves might be biased.
#' The default value is `1000`.
#' @param shap_nsim (`integer`) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value. When the number of variables
#' is large, a smaller shap_nsim could be used. See details in documentation of
#' function \code{\link[fastshap:explain]{explain}} in package `fastshap`.
#' The default is 100.
#' @param visualize (`logical`) if `TRUE`, plot the variable dependence plots.
#' The default is `FALSE`.
#' @param seed (`integer`) The seed for any random progress. The default is `10`.
#' @param pfun (`function`) The predict function that requires two arguments,
#' `object` and `newdata`.
#' It is only required when `model` is not \code{isolation_forest}.
#' The default is the wrapper function designed for iForest model in `itsdm`.
#'
#' @return (`ShapDependence`) A list of
#' \itemize{
#' \item{dependences_cont (`list`) A list of Shapley values of continuous variables}
#' \item{dependences_cat (`list`) A list of Shapley values of categorical variables}
#' \item{feature_values (`data.frame`) A table of feature values}
#' }
#'
#' @seealso
#' \code{\link{plot.ShapDependence}}
#' \code{\link[fastshap:explain]{explain}} in `fastshap`
#'
#' @details
#' The values show how each environmental variable independently
#' affects the modeling prediction. They show how the Shapley value of each variable
#' changes as its value is varied.
#'
#' @references
#' \itemize{
#' \item{Strumbelj, Erik,
#' and Igor Kononenko. "Explaining prediction models and individual predictions
#' with feature contributions." \emph{Knowledge and information systems}
#' 41.3 (2014): 647-665.\doi{10.1007/s10115-013-0679-x}}
#' \item{\href{http://proceedings.mlr.press/v119/sundararajan20b.html}{Sundara
#' rajan, Mukund, and Amir Najmi. "The many Shapley values for model explanation
#' ." \emph{International Conference on Machine Learning}. PMLR, 2020.}}
#' \item{\url{https://github.com/bgreenwell/fastshap}}
#' \item{\url{https://github.com/shap/shap}}
#' }
#'
#' @importFrom dplyr select
#' @importFrom fastshap explain
#' @importFrom stars st_apply
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
#'   variables = env_vars, ntrees = 10,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, nthreads = 1,
#'   response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_dependence <- shap_dependence(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables)
#' plot(var_dependence, target_var = "bio1", related_var = "bio16")
#' }
#'
#' \dontrun{
#' ##### Use Random Forest model as an external model ########
#' library(randomForest)
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
#' shap_dependences <- shap_dependence(
#'   model = mod_rf,
#'   var_occ = model_data %>% select(-occ),
#'   variables = env_vars,
#'   visualize = FALSE,
#'   seed = 10,
#'   pfun = pfun)
#' }
#'
shap_dependence <- function(model,
                            var_occ,
                            variables,
                            si = 1000,
                            shap_nsim = 100,
                            visualize = FALSE,
                            seed = 10,
                            pfun = .pfun_shap) {

  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_class(variables, 'stars')
  checkmate::assert_int(si)
  checkmate::assert_int(shap_nsim)
  checkmate::assert_logical(visualize)
  checkmate::assert_int(seed)

  # Get background samples
  # Because this is response, it is necessary to
  # cover the whole range of values. This is why
  # we need the background samples
  # Make a template
  rst_template <- st_apply(merge(variables), c("x", "y"),
                           "mean", na.rm = FALSE)
  rst_template[[1]][!is.na(rst_template[[1]])] <- 1

  # Observations
  obs_bg <- .bg_sampling(rst_template, NULL, seed, si)
  obs_bg_vars_mat <- st_extract(x = variables, at = obs_bg) %>%
    st_drop_geometry()
  rm(rst_template)

  # Do SHAP
  set.seed(seed)
  shap_explain <- explain(model, X = var_occ, nsim = shap_nsim,
                          newdata = obs_bg_vars_mat,
                          pred_wrapper = pfun)
  shap_explain <- as.data.frame(shap_explain) # For fastshap >= 0.1.0
  dependences <- lapply(names(var_occ), function(var) {
    data.frame(x = obs_bg_vars_mat %>% pull(var),
               y = shap_explain %>% pull(var))
  })
  names(dependences) <- names(var_occ)

  # Split numeric and categorical
  isfacor <- as.vector(sapply(var_occ, is.factor))
  bands_cont <- names(var_occ)[!isfacor]
  bands_cat <- names(var_occ)[isfacor]

  if (length(bands_cont) > 0) {
    dependences_cont <- dependences[bands_cont]
  } else dependences_cont <- NULL
  if (length(bands_cat) > 0) {
    dependences_cat <- dependences[bands_cat]
  } else dependences_cat <- NULL

  # Reunion
  dependences <- list(dependences_cont = dependences_cont,
                      dependences_cat = dependences_cat)
  dependences$feature_values <- obs_bg_vars_mat
  class(dependences) <- append("ShapDependence", class(dependences))

  # Visualize
  if (visualize) {
    print(plot(dependences))
  }

  # Return
  return(dependences)
}

# shap_dependence end
