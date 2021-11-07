#' @title A function to calculate variable dependence based on SHAP method.
#' @description This function allows you to calculate the variable dependence
#' based on Shapley values.
#' @param model Any predictive model. In this package, it could be isolation_forest.
#' @param var_occ (data.frame, tibble) the data.frame style table that
#' include values of environmental variables at occurrence locations.
#' @param shap_nsim (integer) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @param seed (integer) The seed for random progress.
#' @details The curves show how each environmental variable independently
#' affects the modeling prediction. The curves show how the predicted result
#' only using this variable changes as it is varied.
#' @return (VariableDependence) a list of data.frame of dependence and variable values.
#' The dependence values correspond to Shapley value.
#' @importFrom dplyr select
#' @importFrom fastshap explain
#' @export
#' @examples
#' variable_dependence(model = mod$model, var_occ = var_occ)
variable_dependence <- function(model,
                                var_occ,
                                shap_nsim = 100,
                                visualize = FALSE,
                                seed = 10) {

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
  dependences$feature_values <- var_occ
  class(dependences) <- append("VariableDependence", class(dependences))

  # Visualize
  if (visualize) {
    plot(dependences)
  }

  # Return
  return(dependences)
}

# independent_response end
