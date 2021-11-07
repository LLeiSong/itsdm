#' A function to evaluate relative importance of each variable.
#' @description This function allows you to evaluate relative importance of
#' each variable within the model using methods: with this variable only,
#' with all variables except this variables, and permutation feature importance
#' method. Because the model is presence-only, we only use pearson correlation
#' with the result of full model as the metrics to rank the variables.
#' @param model (isolation_forest) The extended isolation forest SDM.
#' @param var_occ (data.frame, tibble) the data.frame style table that
#' include values of environmental variables at occurrence locations.
#' @param var_occ_analysis (data.frame, tibble) The data.frame style table that
#' include values of environmental variables at occurrence locations of analysis.
#' @param shap_nsim (integer) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @param seed (integer) The seed for any random progress. The default is 10.
#' @return (VariableContribution) a tibble of Shapley values of each variable for
#' each observation.
#' @importFrom fastshap explain
#' @export
#' @examples
#' variable_contrib(model = isotree_po,
#' var_occ = var_occ, variables = env_vars)
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
