#' @title Evaluate variable contributions for targeted observations.
#' @description Evaluate variable contribution for targeted observations according
#' to SHapley Additive exPlanations (SHAP).
#' @param model (\code{isolation_forest} or other model) The SDM.
#' It could be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' It also could be other user-fitted models as long as the `pfun` can work on it.
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
#' @param pfun (`function`) The predict function that requires two arguments,
#' `object` and `newdata`.
#' It is only required when `model` is not \code{isolation_forest}.
#' The default is the wrapper function designed for iForest model in `itsdm`.
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
#'   variables = env_vars, ntrees = 5,
#'   sample_size = 0.8, ndim = 1L,
#'   seed = 123L, nthreads = 1,
#'   response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_contribution <- variable_contrib(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   var_occ_analysis = mod$vars_train %>% slice(1:2))
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
                             seed = 10,
                             pfun = .pfun_shap) {
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
                 pred_wrapper = pfun)
  out <- as.data.frame(out) # For fastshap >= 0.1.0

  out <- list(shapley_values = out,
              feature_values = var_occ_analysis)
  class(out) <- append("VariableContribution", class(out))

  # Plot and return
  if (visualize) print(plot(out))
  return(out)
}
