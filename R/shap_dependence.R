#' @title Calculate Shapley value based variable dependence.
#' @description Calculate the variable dependence using Shapley values.
#' @param model (\code{isolation_forest}). The isolation forest SDM.
#' It could be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' @param var_occ (`data.frame`, `tibble`) The `data.frame` style table that
#' include values of environmental variables at occurrence locations.
#' @param shap_nsim (`integer`) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value. When the number of variables
#' is large, a smaller shap_nsim could be used. See details in documentation of
#' function \code{\link{explain}} in package `fastshap`.
#' The default is 100.
#' @param visualize (`logical`) if `TRUE`, plot the variable dependence plots.
#' The default is `FALSE`.
#' @param seed (`integer`) The seed for any random progress. The default is `10L`.
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
#' \code{\link{explain}} in `fastshap`
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
#' \item{\url{https://github.com/slundberg/shap}}
#' }
#'
#' @importFrom dplyr select
#' @importFrom fastshap explain
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' mod <- isotree_po(
#'   occ = occ, occ_test = occ_test,
#'   variables = env_vars, ntrees = 50,
#'   sample_rate = 0.8, ndim = 1L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_dependence <- shap_dependence(
#'   model = mod$model,
#'   var_occ = mod$var_train %>% st_drop_geometry())
#'}
#'
shap_dependence <- function(model,
                            var_occ,
                            shap_nsim = 100,
                            visualize = FALSE,
                            seed = 10L) {

  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_int(shap_nsim)
  checkmate::assert_logical(visualize)
  checkmate::assert_int(seed)

  # Do SHAP
  set.seed(seed)
  shap_explain <- explain(model, X = var_occ, nsim = shap_nsim,
                          pred_wrapper = .pfun_shap)
  dependences <- lapply(names(var_occ), function(var) {
    data.frame(x = var_occ %>% pull(var),
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
  dependences$feature_values <- var_occ
  class(dependences) <- append("ShapDependence", class(dependences))

  # Visualize
  if (visualize) {
    plot(dependences)
  }

  # Return
  return(dependences)
}

# shap_dependence end
