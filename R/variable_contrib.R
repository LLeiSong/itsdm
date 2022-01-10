#' @title Evaluate variable contributions for targeted observations.
#' @description Evaluate variable contribution for targeted observations according
#' to SHapley Additive exPlanations (SHAP).
#' @param model (\code{isolation_forest}) The isolation forest SDM.
#' It could be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' @param var_occ (`data.frame`, `tibble`) The `data.frame` style table that
#' include values of environmental variables at occurrence locations.
#' @param var_occ_analysis (`data.frame`, `tibble`) The `data.frame` style table that
#' include values of environmental variables at occurrence locations for analysis. It
#' could be either `var_occ` or its subset, or any new dataset.
#' @param shap_nsim (`integer`) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value. See details in documentation of
#' function \code{\link{explain}} in package `fastshap`.
#' @param visualize (`logical`) if `TRUE`, plot the response curves.
#' The default is `FALSE`.
#' @param seed (`integer`) The seed for any random progress. The default is `10L`.
#' @return (`VariableContribution`) A list of
#' \itemize{
#' \item{shapley_values (`data.frame`) A table of Shapley values of each variables for
#' all observations}
#' \item{feature_values (`tibble`) A table of values of each variables for all
#' observations}}
#'
#' @seealso
#' \code{\link{plot.VariableContribution}}
#' \code{\link{explain}} in `fastshap`
#'
#' @references
#' \itemize{
#' \item{\href{https://github.com/bgreenwell/fastshap}{https://github.com/
#' bgreenwell/fastshap}}
#' \item{\href{https://github.com/slundberg/shap}{https://github.com/slundberg/shap}}
#' }
#'
#' @importFrom fastshap explain
#' @export
#' @examples
#'
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
#'   sample_size = 0.8, ndim = 3L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_contribution <- variable_contrib(
#'   model = mod$model,
#'   var_occ = mod$var_train %>% st_drop_geometry(),
#'   var_occ_analysis = mod$var_train %>%
#'     st_drop_geometry() %>% slice(1:10))
#'\dontrun{
#' plot(var_contribution,
#'   num_features = 3,
#'   plot_each_obs = TRUE)
#'
#' # Plot together
#' plot(var_contribution)
#'}
#'
variable_contrib <- function(model,
                             var_occ, # Training, must set for model
                             var_occ_analysis,
                             shap_nsim = 100,
                             visualize = FALSE,
                             seed = 10) {
  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_data_frame(var_occ_analysis)
  checkmate::assert_int(shap_nsim)
  checkmate::assert_logical(visualize)
  checkmate::assert_int(seed)
  stopifnot(identical(colnames(var_occ_analysis),
                      setdiff(colnames(var_occ), 'var_occ')))

  # Use SHAP
  set.seed(seed)
  out <- explain(model, X = var_occ,
                 nsim = shap_nsim,
                 newdata = var_occ_analysis,
                 pred_wrapper = .pfun_shap)
  out <- list(shapley_values = out,
              feature_values = var_occ_analysis)
  class(out) <- append("VariableContribution", class(out))

  # Plot and return
  if (visualize) plot(out)
  return(out)
}
