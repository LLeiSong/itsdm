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
#' @param permution_cv (integer) The time to run permutation method.
#' The default is 5.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @return variable_importance.
#' @import dplyr
#' @import stars
#' @export
#' @examples
#' variable_importance(model = isotree_po,
#' var_occ = var_occ, variables = env_vars)
variable_importance <- function(model,
                                var_occ,
                                var_occ_bg, # Background points
                                var_occ_test = NULL, # independent test
                                var_occ_test_bg = NULL,
                                variables,
                                permution_cv = 5,
                                visualize = FALSE,
                                seed = 10) {
  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_data_frame(var_occ_test)
  checkmate::assert_multi_class(
    variables, c('RasterStack', 'RasterLayer', 'stars'))
  checkmate::assert_int(permution_cv)
  checkmate::assert_logical(visualize)
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
  full_pred_occ_bg <- 1 - predict(model, var_occ_bg)
  full_pred_occ_test <- 1 - predict(model, var_occ_test)
  full_pred_occ_test_bg <- 1 - predict(model, var_occ_test_bg)
  full_pred_var <- 1 - predict(model, vars)

  ## Stretch to 0 to 1.
  full_pred_occ <- .stretch(x = full_pred_var,
                            new_values = full_pred_occ)
  full_pred_occ_bg <- .stretch(x = full_pred_var,
                               new_values = full_pred_occ_bg)
  full_pred_occ_test <- .stretch(x = full_pred_var,
                                 new_values = full_pred_occ_test)
  full_pred_occ_test_bg <- .stretch(x = full_pred_var,
                                    new_values = full_pred_occ_test_bg)
  full_pred_var <- .stretch(x = full_pred_var)

  ## AUC background and ratio
  full_auc_train <- .auc_ratio(full_pred_occ, full_pred_var)
  full_auc_test <- .auc_ratio(full_pred_occ_test, full_pred_var)
  full_auc_pa_train <- rocit(score = c(full_pred_occ, full_pred_occ_bg),
                             class = c(rep(1, length(full_pred_occ)),
                                       rep(0, length(full_pred_occ_bg))))$AUC
  full_auc_pa_test <- rocit(score = c(full_pred_occ_test, full_pred_occ_test_bg),
                            class = c(rep(1, length(full_pred_occ_test)),
                                      rep(0, length(full_pred_occ_test_bg))))$AUC

  # Process
  out <- do.call(rbind, lapply(bands, function(nm){
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
    this_occ_bg_pred <- 1 - predict(this_model,
                                    var_occ_bg %>% select(nm))
    this_occ_test_pred <- 1 - predict(this_model, this_var_occ_test)
    this_occ_test_bg_pred <- 1 - predict(this_model,
                                         var_occ_test_bg %>% select(nm))
    this_vars_pred <- 1 - predict(this_model, this_vars)

    # Stretch
    this_occ_pred <- .stretch(x = this_vars_pred,
                              new_values = this_occ_pred)
    this_occ_bg_pred <- .stretch(x = this_vars_pred,
                              new_values = this_occ_bg_pred)
    this_occ_test_pred <- .stretch(x = this_vars_pred,
                                   new_values = this_occ_test_pred)
    this_occ_test_bg_pred <- .stretch(x = this_vars_pred,
                                   new_values = this_occ_test_bg_pred)
    this_vars_pred <- .stretch(this_vars_pred)

    ## Calculate metrics
    r_only_train <- cor(full_pred_occ, this_occ_pred,
                        use = 'complete.obs')
    r_only_test <- cor(full_pred_occ_test, this_occ_test_pred,
                       use = 'complete.obs')
    auc_only_train <- .auc_ratio(this_occ_pred, this_vars_pred)
    auc_only_test <- .auc_ratio(this_occ_test_pred, this_vars_pred)
    auc_only_pa_train <- rocit(score = c(this_occ_pred, this_occ_bg_pred),
                               class = c(rep(1, length(this_occ_pred)),
                                         rep(0, length(this_occ_bg_pred))))$AUC
    auc_only_pa_test <- rocit(score = c(this_occ_test_pred, this_occ_test_bg_pred),
                              class = c(rep(1, length(this_occ_test_pred)),
                                        rep(0, length(this_occ_test_bg_pred))))$AUC

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
    except_occ_bg_pred <- 1 - predict(except_model,
                                      var_occ_bg %>% select(all_of(nms)))
    except_occ_test_pred <- 1 - predict(except_model, except_var_occ_test)
    except_occ_test_bg_pred <- 1 - predict(except_model,
                                      var_occ_test_bg %>% select(all_of(nms)))
    except_vars_pred <- 1 - predict(except_model, except_vars)

    ## Stretch
    except_occ_pred <- .stretch(x = except_vars_pred,
                                new_values = except_occ_pred)
    except_occ_bg_pred <- .stretch(x = except_vars_pred,
                                   new_values = except_occ_bg_pred)
    except_occ_test_pred <- .stretch(x = except_vars_pred,
                                     new_values = except_occ_test_pred)
    except_occ_test_bg_pred <- .stretch(x = except_vars_pred,
                                        new_values = except_occ_test_bg_pred)
    except_vars_pred <- .stretch(except_vars_pred)

    ## Calculate metrics
    r_except_train <- cor(full_pred_occ, except_occ_pred,
                          use = 'complete.obs')
    r_except_test <- cor(full_pred_occ_test, except_occ_test_pred,
                          use = 'complete.obs')
    auc_except_train <- .auc_ratio(except_occ_pred, except_vars_pred)
    auc_except_test <- .auc_ratio(except_occ_test_pred, except_vars_pred)
    auc_except_pa_train <- rocit(score = c(except_occ_pred, except_occ_bg_pred),
                               class = c(rep(1, length(except_occ_pred)),
                                         rep(0, length(except_occ_bg_pred))))$AUC
    auc_except_pa_test <- rocit(score = c(except_occ_test_pred, except_occ_test_bg_pred),
                              class = c(rep(1, length(except_occ_test_pred)),
                                        rep(0, length(except_occ_test_bg_pred))))$AUC

    # Model with variables with permuted data
    r_permuted <- do.call(
      rbind, lapply(1:permution_cv, function(cv_n) {
      ## Permute variable
      set.seed(11 + cv_n)
      permuted_var_occ <- var_occ %>% select(-nm) %>%
        mutate(!!nm := sample(var_occ %>% pull(nm), nrow(var_occ)))
      set.seed(12 + cv_n)
      permuted_var_occ_bg <- var_occ_bg %>% select(-nm) %>%
        mutate(!!nm := sample(var_occ_bg %>% pull(nm), nrow(var_occ_bg)))
      set.seed(13 + cv_n)
      permuted_var_occ_test <- var_occ_test %>% select(-nm) %>%
        mutate(!!nm := sample(var_occ_test %>% pull(nm), nrow(var_occ_test)))
      set.seed(14 + cv_n)
      permuted_var_occ_test_bg <- var_occ_test_bg %>% select(-nm) %>%
        mutate(!!nm := sample(var_occ_test_bg %>% pull(nm), nrow(var_occ_test_bg)))
      set.seed(15 + cv_n)
      permuted_vars <- vars %>% select(-nm) %>%
        mutate(!!nm := sample(vars %>% pull(nm), nrow(vars)))

      ## Fit model
      permuted_model <- isolation.forest(
        permuted_var_occ,
        ntrees = model$params$ntrees,
        sample_size = model$params$sample_size,
        ndim = model$params$ndim,
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
      permuted_occ_pred <- 1 - predict(permuted_model, permuted_var_occ)
      permuted_occ_bg_pred <- 1 - predict(permuted_model, permuted_var_occ_bg)
      permuted_occ_test_pred <- 1 - predict(permuted_model, permuted_var_occ_test)
      permuted_occ_test_bg_pred <- 1 - predict(permuted_model, permuted_var_occ_test_bg)
      permuted_vars_pred <- 1 - predict(permuted_model, permuted_vars)

      ## Stretch
      permuted_occ_pred <- .stretch(x = permuted_vars_pred,
                                    new_values = permuted_occ_pred)
      permuted_occ_bg_pred <- .stretch(x = permuted_vars_pred,
                                       new_values = permuted_occ_bg_pred)
      permuted_occ_test_pred <- .stretch(x = permuted_vars_pred,
                                         new_values = permuted_occ_test_pred)
      permuted_occ_test_bg_pred <- .stretch(x = permuted_vars_pred,
                                            new_values = permuted_occ_test_bg_pred)
      permuted_vars_pred <- .stretch(permuted_vars_pred)

      ## Calculate metrics
      r_permuted_train <- cor(full_pred_occ, permuted_occ_pred,
                              use = 'complete.obs')
      r_permuted_test <- cor(full_pred_occ_test, permuted_occ_test_pred,
                             use = 'complete.obs')
      auc_permuted_train <- .auc_ratio(permuted_occ_pred, permuted_vars_pred)
      auc_permuted_test <- .auc_ratio(permuted_occ_test_pred, permuted_vars_pred)
      auc_permuted_pa_train <- rocit(score = c(permuted_occ_pred, permuted_occ_bg_pred),
                                   class = c(rep(1, length(permuted_occ_pred)),
                                             rep(0, length(permuted_occ_bg_pred))))$AUC
      auc_permuted_pa_test <- rocit(score = c(permuted_occ_test_pred, permuted_occ_test_bg_pred),
                                  class = c(rep(1, length(permuted_occ_test_pred)),
                                            rep(0, length(permuted_occ_test_bg_pred))))$AUC

      # Clean up
      rm(permuted_var_occ, permuted_model, permuted_occ_pred, permuted_occ_bg_pred,
         permuted_occ_test_pred, permuted_occ_test_bg_pred, permuted_vars_pred)
      tibble(cv = cv_n,
             cor_train = 1 - r_permuted_train,
             cor_test = 1 - r_permuted_test,
             auc_ratio_train = auc_permuted_train,
             auc_ratio_test = auc_permuted_test,
             auc_train = auc_permuted_pa_train,
             auc_test = auc_permuted_pa_test)
    }))

    # Get confident interval of permutation result
    ## Correlation
    r_mean_train <- mean(r_permuted$cor_train)
    r_mean_test <- mean(r_permuted$cor_test)

    ## AUC
    ar_mean_train <- mean(r_permuted$auc_ratio_train)
    ar_mean_test <- mean(r_permuted$auc_ratio_test)
    a_mean_train <- mean(r_permuted$auc_train)
    a_mean_test <- mean(r_permuted$auc_test)

    # Clean up
    rm(this_var_occ, this_var_occ_test, this_vars, this_model, this_occ_pred,
       this_occ_test_pred, nms, except_var_occ, except_var_occ_test, except_vars,
       except_model, except_occ_pred, except_occ_test_pred, r_permuted)

    # Output
    tibble(variable = rep(nm, 18),
           metrics = c(rep('Pearson correlation', 6), rep('AUC ratio', 6),
                       rep('AUC', 6)),
           method = rep(c('Only', 'Except', 'Permutation'), 6),
           usage = rep(c(rep('Train', 3), rep('Test', 3)), 3),
           value = c(r_only_train, r_except_train, r_mean_train,
                     r_only_test, r_except_test, r_mean_test,
                     full_auc_train - auc_only_train,
                     full_auc_train - auc_except_train,
                     full_auc_train - ar_mean_train,
                     full_auc_test - auc_only_test,
                     full_auc_test - auc_except_test,
                     full_auc_test - ar_mean_test,
                     full_auc_pa_train - auc_only_pa_train,
                     full_auc_pa_train - auc_except_pa_train,
                     full_auc_pa_train - a_mean_train,
                     full_auc_pa_test - auc_only_pa_test,
                     full_auc_pa_test - auc_except_pa_test,
                     full_auc_pa_test - a_mean_test))
  }))


  class(out) <- append("variable_importance", class(out))

  # Visualize
  if (visualize) {
    communicate(out)
    plot(out)
  }

  # Return
  return(out)
}

# variable_importance end
