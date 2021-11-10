#' A function to evaluate relative importance of each variable.
#' @description This function allows you to evaluate relative importance of
#' each variable within the model using methods: with this variable only,
#' with all variables except this variables, and permutation feature importance
#' method. Because the model is presence-only, we only use pearson correlation
#' with the result of full model as the metrics to rank the variables.
#' @param model (isolation_forest) The extended isolation forest SDM.
#' @param var_occ (data.frame, tibble) the data.frame style table that
#' include values of environmental variables at occurrence locations.
#' @param var_occ_test (data.frame, tibble) The data.frame style table that
#' include values of environmental variables at occurrence locations of testing.
#' If NULL, no test will be used. The default is NULL.
#' @param variables (RasterStack or stars) The stack of environmental variables.
#' @param shap_nsim (integer) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @param seed (integer) The seed for any random progress. The default is 10.
#' @return (VariableAnalysis) a list of variable importance analysis according
#' to different metrics including jackknife of Pearson correlation and AUC_ratio,
#' and SHAP test.
#' @importFrom dplyr select tibble filter sample_n
#' @importFrom sf st_as_sf st_drop_geometry
#' @importFrom stars st_as_stars st_xy2sfc st_get_dimension_values
#' @importFrom fastshap explain
#' @export
#' @examples
#' variable_analysis(model = isotree_po,
#' var_occ = var_occ, variables = env_vars)
variable_analysis <- function(model,
                              var_occ,
                              var_occ_test, # Independent test
                              variables,
                              shap_nsim = 100,
                              visualize = FALSE,
                              seed = 10) {
  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_data_frame(var_occ_test)
  checkmate::assert_multi_class(
    variables, c('RasterStack', 'RasterLayer', 'stars'))
  checkmate::assert_int(shap_nsim)
  checkmate::assert_logical(visualize)
  checkmate::assert_int(seed)
  bands <- st_get_dimension_values(variables, 'band')
  stopifnot(all(bands %in% colnames(var_occ)))

  # Reformat data
  # Variables, use stars
  if (is(variables, 'RasterStack') |
      is(variables, 'RasterLayer')){
    variables <- st_as_stars(variables)}

  # Sampling from the whole image to speed things up.
  set.seed(seed)
  samples <- variables %>% slice('band', 1) %>%
    st_xy2sfc(as_points = T) %>% st_as_sf() %>%
    select(geometry) %>%
    sample_n(min(10000, nrow(.)))
  vars <- st_extract(x = split(variables, 'band'), at = samples) %>%
    st_drop_geometry()
  rm(samples, variables)

  # Original prediction
  ## Predict
  full_pred_occ <- 1 - predict(model, var_occ)
  full_pred_occ_test <- 1 - predict(model, var_occ_test)
  full_pred_var <- 1 - predict(model, vars)

  ## Stretch to 0 to 1.
  full_pred_occ <- .stretch(x = full_pred_var,
                            new_values = full_pred_occ)
  full_pred_occ_test <- .stretch(x = full_pred_var,
                                 new_values = full_pred_occ_test)
  full_pred_var <- .stretch(x = full_pred_var)

  ## AUC background and ratio
  full_auc_train <- .auc_ratio(full_pred_occ, full_pred_var)
  full_auc_test <- .auc_ratio(full_pred_occ_test, full_pred_var)

  # Process
  var_each <- do.call(rbind, lapply(bands, function(nm){
    # Model with only this variable
    ## Subset dataset
    this_var_occ <- var_occ %>% select(nm)
    this_var_occ_test <- var_occ_test %>% select(nm)
    this_vars <- vars %>% select(nm)

    ## Fit model
    this_model <- isolation.forest(
      this_var_occ,
      ntrees = model$params$ntrees,
      sample_size = model$params$sample_size,
      ndim = 1,
      ntry = model$params$ntry,
      max_depth = model$params$max_depth,
      prob_pick_avg_gain = model$params$prob_pick_avg_gain,
      prob_pick_pooled_gain = model$params$prob_pick_pooled_gain,
      prob_split_avg_gain = model$params$prob_split_avg_gain,
      prob_split_pooled_gain = model$params$prob_split_pooled_gain,
      min_gain = model$params$min_gain,
      missing_action = model$params$missing_action,
      categ_split_type = model$params$categ_split_type,
      all_perm = model$params$all_perm,
      coef_by_prop = model$params$coef_by_prop,
      weights_as_sample_prob = model$params$weights_as_sample_prob,
      sample_with_replacement = model$params$sample_with_replacement,
      penalize_range = model$params$penalize_range,
      weigh_by_kurtosis = model$params$weigh_by_kurtosis,
      coefs = model$params$coefs,
      assume_full_distr = model$params$assume_full_distr,
      build_imputer = model$params$build_imputer,
      min_imp_obs = model$params$min_imp_obs,
      depth_imp = model$params$depth_imp,
      weigh_imp_rows = model$params$weigh_imp_rows)

    ## Prediction
    this_occ_pred <- 1 - predict(this_model, this_var_occ)
    this_occ_test_pred <- 1 - predict(this_model, this_var_occ_test)
    this_vars_pred <- 1 - predict(this_model, this_vars)

    # Stretch
    this_occ_pred <- .stretch(x = this_vars_pred,
                              new_values = this_occ_pred)
    this_occ_test_pred <- .stretch(x = this_vars_pred,
                                   new_values = this_occ_test_pred)
    this_vars_pred <- .stretch(this_vars_pred)

    ## Calculate metrics
    r_only_train <- cor(full_pred_occ, this_occ_pred,
                        use = 'complete.obs')
    r_only_test <- cor(full_pred_occ_test, this_occ_test_pred,
                       use = 'complete.obs')
    auc_only_train <- .auc_ratio(this_occ_pred, this_vars_pred)
    auc_only_test <- .auc_ratio(this_occ_test_pred, this_vars_pred)

    # Model with variables except the chosen one
    ## Subset dataset
    nms <- setdiff(bands, nm)
    except_var_occ <- var_occ %>% select(all_of(nms))
    except_var_occ_test <- var_occ_test %>% select(all_of(nms))
    except_vars <- vars %>% select(all_of(nms))

    ## Fit model
    except_model <- isolation.forest(
      except_var_occ,
      ntrees = model$params$ntrees,
      sample_size = model$params$sample_size,
      ndim = ifelse(ncol(except_var_occ) >= model$params$ndim,
                    model$params$ndim, ncol(except_var_occ)),
      ntry = model$params$ntry,
      max_depth = model$params$max_depth,
      prob_pick_avg_gain = model$params$prob_pick_avg_gain,
      prob_pick_pooled_gain = model$params$prob_pick_pooled_gain,
      prob_split_avg_gain = model$params$prob_split_avg_gain,
      prob_split_pooled_gain = model$params$prob_split_pooled_gain,
      min_gain = model$params$min_gain,
      missing_action = model$params$missing_action,
      categ_split_type = model$params$categ_split_type,
      all_perm = model$params$all_perm,
      coef_by_prop = model$params$coef_by_prop,
      weights_as_sample_prob = model$params$weights_as_sample_prob,
      sample_with_replacement = model$params$sample_with_replacement,
      penalize_range = model$params$penalize_range,
      weigh_by_kurtosis = model$params$weigh_by_kurtosis,
      coefs = model$params$coefs,
      assume_full_distr = model$params$assume_full_distr,
      build_imputer = model$params$build_imputer,
      min_imp_obs = model$params$min_imp_obs,
      depth_imp = model$params$depth_imp,
      weigh_imp_rows = model$params$weigh_imp_rows)

    ## Prediction
    except_occ_pred <- 1 - predict(except_model, except_var_occ)
    except_occ_test_pred <- 1 - predict(except_model, except_var_occ_test)
    except_vars_pred <- 1 - predict(except_model, except_vars)

    ## Stretch
    except_occ_pred <- .stretch(x = except_vars_pred,
                                new_values = except_occ_pred)
    except_occ_test_pred <- .stretch(x = except_vars_pred,
                                     new_values = except_occ_test_pred)
    except_vars_pred <- .stretch(except_vars_pred)

    ## Calculate metrics
    r_except_train <- cor(full_pred_occ, except_occ_pred,
                          use = 'complete.obs')
    r_except_test <- cor(full_pred_occ_test, except_occ_test_pred,
                          use = 'complete.obs')
    auc_except_train <- .auc_ratio(except_occ_pred, except_vars_pred)
    auc_except_test <- .auc_ratio(except_occ_test_pred, except_vars_pred)

    # Clean up
    rm(this_var_occ, this_var_occ_test, this_vars, this_model, this_occ_pred,
       this_occ_test_pred, nms, except_var_occ, except_var_occ_test, except_vars,
       except_model, except_occ_pred, except_occ_test_pred)

    # Output
    tibble(variable = rep(nm, 8),
           metrics = c(rep('pearson_correlation', 4), rep('AUC_ratio', 4)),
           method = rep(c(rep('Only', 2), rep('Without', 2)), 2),
           usage = rep(c('Train', 'Test'), 4),
           value = c(r_only_train, r_only_test,
                     r_except_train, r_except_test,
                     auc_only_train, auc_only_test,
                     auc_except_train, auc_except_test))
  }))

  # NOTE: Permutation feature importance method does not work for isolation forest
  # because each split is independent to the rank of values. Its tree structure is
  # fundamentally different from tree structure such as random forest.
  # ALTERNATIVE: SHAP method that use Shapley values according to game theory.

  # SHAP feature importance
  set.seed(seed)
  shap_train <- explain(model, X = var_occ, nsim = shap_nsim,
                        pred_wrapper = .pfun_shap)
  set.seed(seed)
  shap_test <- explain(model, X = var_occ, nsim = shap_nsim,
                       newdata = var_occ_test,
                       pred_wrapper = .pfun_shap)

  # Output
  out <- list(variables = bands,
              pearson_correlation =
                var_each %>%
                filter(metrics == 'pearson_correlation') %>%
                select(-metrics),
              full_AUC_ratio =
                tibble(full_auc_train = full_auc_train,
                       full_auc_test = full_auc_test),
              AUC_ratio =
                var_each %>%
                filter(metrics == 'AUC_ratio') %>%
                select(-metrics),
              SHAP = list(train = shap_train,
                          test = shap_test))

  class(out) <- append("VariableAnalysis", class(out))

  # Visualize
  if (visualize) {
    plot(out)
  }

  # Return
  return(out)
}

# variable_analysis end
