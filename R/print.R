#' @title Print summary information from variable importance object.
#' @description Display the most general and informative characteristics of
#' a variable importance object.
#' @param x (VariableAnalysis) A variable importance object to be messaged.
#' It could be the return of function `variable_analysis`.
#' @param ... Not used.
#' @importFrom dplyr filter pull
#' @importFrom stringr str_pad
#' @return The same object that was passed as input.
#' @export
#' @examples
#' print(variable_analysis)
#'
print.VariableAnalysis <- function(x, ...){
  cat('Relative variable importance\n')
  cat(paste0(c(rep('=', 35), '\n'), collapse = ''))
  cat('Methods: Jackknife test and SHAP\n')
  cat(sprintf('Numer of variables: %d\n', length(x$variables)))
  cat(paste0(c(rep('=', 35), '\n'), collapse = ''))

  # Pearson correlation
  cat('Jackknife test\n')
  cat('Based on Pearson correlation (Max value is 1)\n')
  cor_x <- x$pearson_correlation
  cat('[Training dataset]:\n')
  cor_x_train <- cor_x %>%
    filter(usage == 'Train') %>%
    arrange(-value)

  n_max <- max(nchar(unique(cor_x_train$variable)))
  invisible(lapply(unique(cor_x_train %>%
                            filter(method == 'Only') %>%
                            pull(variable)), function(var){
    this_var <- cor_x_train %>% filter(variable == var)
    only <- this_var %>% filter(method == 'Only') %>% pull(value)
    without <- this_var %>% filter(method == 'Without') %>% pull(value)

    # print
    var <- ifelse(length(var) >= 20, var[1:20], var)
    n_space <- min(n_max, 20)
    cat(paste0(str_pad(var, n_space, side = 'right', pad = ' '),
               ' With only: ',
               paste(rep('#', round(only * 45)), collapse = ''),
               sprintf(' %s\n', round(only, 3))))
    cat(paste0(str_pad('', n_space, side = 'right', pad = ' '), ' Without  : ',
               paste(rep('#', round(without * 45)), collapse = ''),
               sprintf(' %s\n', round(without, 3))))
  }))

  cat('[Test dataset]:\n')
  cor_x_test <- cor_x %>%
    filter(usage == 'Test') %>%
    arrange(-value)
  invisible(lapply(unique(cor_x_test %>%
                  filter(method == 'Only') %>%
                  pull(variable)), function(var){
    this_var <- cor_x_test %>% filter(variable == var)
    only <- this_var %>% filter(method == 'Only') %>% pull(value)
    without <- this_var %>% filter(method == 'Without') %>% pull(value)

    # print
    var <- ifelse(length(var) >= 20, var[1:20], var)
    n_space <- min(n_max, 20)
    cat(paste0(str_pad(var, n_space, side = 'right', pad = ' '),
               ' With only: ',
               paste(rep('#', round(only * 45)), collapse = ''),
               sprintf(' %s\n', round(only, 3))))
    cat(paste0(str_pad('', n_space, side = 'right', pad = ' '), ' Without  : ',
               paste(rep('#', round(without * 45)), collapse = ''),
               sprintf(' %s\n', round(without, 3))))
  }))
  cat(paste0(c(rep('=', 70), '\n'), collapse = ''))

  # AUC ratio
  auc_x <- x$AUC_ratio
  full_auc_x <- x$full_AUC_ratio
  cat('Jackknife test\n')
  cat(sprintf('Based on AUC ratio (Max value of traing and test are %.3f and %.3f)\n',
      full_auc_x$full_auc_train, full_auc_x$full_auc_test))
  cat('[Training dataset]:\n')
  auc_x_train <- auc_x %>%
    filter(usage == 'Train') %>%
    arrange(-value)

  n_max <- max(nchar(unique(auc_x_train$variable)))
  invisible(lapply(unique(auc_x_train %>%
                            filter(method == 'Only') %>%
                            pull(variable)), function(var){
    this_var <- auc_x_train %>% filter(variable == var)
    only <- this_var %>% filter(method == 'Only') %>% pull(value)
    without <- this_var %>% filter(method == 'Without') %>% pull(value)

    # print
    var <- ifelse(length(var) >= 20, var[1:20], var)
    n_space <- min(n_max, 20)
    cat(paste0(str_pad(var, n_space, side = 'right', pad = ' '),
               ' With only: ',
               paste(rep('#', round(only * 45)), collapse = ''),
               sprintf(' %s\n', round(only, 3))))
    cat(paste0(str_pad('', n_space, side = 'right', pad = ' '), ' Without  : ',
               paste(rep('#', round(without * 45)), collapse = ''),
               sprintf(' %s\n', round(without, 3))))
  }))

  cat('[Test dataset]:\n')
  auc_x_test <- auc_x %>%
    filter(usage == 'Test') %>%
    arrange(-value)
  invisible(lapply(unique(auc_x_test %>%
                            filter(method == 'Only') %>%
                            pull(variable)), function(var){
    this_var <- auc_x_test %>% filter(variable == var)
    only <- this_var %>% filter(method == 'Only') %>% pull(value)
    without <- this_var %>% filter(method == 'Without') %>% pull(value)

    # print
    var <- ifelse(length(var) >= 20, var[1:20], var)
    n_space <- min(n_max, 20)
    cat(paste0(str_pad(var, n_space, side = 'right', pad = ' '),
               ' With only: ',
               paste(rep('#', round(only * 45)), collapse = ''),
               sprintf(' %s\n', round(only, 3))))
    cat(paste0(str_pad('', n_space, side = 'right', pad = ' '), ' Without  : ',
               paste(rep('#', round(without * 45)), collapse = ''),
               sprintf(' %s\n', round(without, 3))))
  }))
  cat(paste0(c(rep('=', 70), '\n'), collapse = ''))

  # SHAP
  ## Subset
  cat('SHAP (mean(|Shapley value|))\n')
  shap_x <- x$SHAP
  shap_x_train <- shap_x$train %>%
    summarise(across(all_of(names(.)), .abs_mean))
  shap_x_train <- t(
    apply(shap_x_train, 1,
          function(x) sort(x, decreasing = T))) %>%
    as.data.frame()
  shap_x_test <- shap_x$test %>%
    summarise(across(all_of(names(.)), .abs_mean))
  shap_x_test <- t(
    apply(shap_x_test, 1,
          function(x) sort(x, decreasing = T))) %>%
    as.data.frame()

  ## Training
  cat('[Training dataset]:\n')
  n_max <- max(nchar(names(shap_x_train)))
  val_max <- max(max(shap_x_train), max(shap_x_test))
  invisible(lapply(names(shap_x_train), function(var){
    val <- shap_x_train %>% pull(all_of(var))
    # print
    var <- ifelse(length(var) >= 20, var[1:20], var)
    n_space <- min(n_max, 20)
    cat(paste0(str_pad(var, n_space, side = 'right', pad = ' '),
               ' : ',
               paste(rep('#', round(val / val_max * 45)), collapse = ''),
               sprintf(' %s\n', round(val, 3))))
  }))

  ## Test
  cat('[Test dataset]:\n')
  invisible(lapply(names(shap_x_test), function(var){
    val <- shap_x_test %>% pull(all_of(var))
    # print
    var <- ifelse(length(var) >= 20, var[1:20], var)
    n_space <- min(n_max, 20)
    cat(paste0(str_pad(var, n_space, side = 'right', pad = ' '),
               ' : ',
               paste(rep('#', round(val / val_max * 45)), collapse = ''),
               sprintf(' %s\n', round(val, 3))))
  }))

  invisible(return(x))
}

#' @title Print summary information from presence-only evaluation object.
#' @description Display the most general and informative characteristics of
#' a presence-only evaluation object.
#' @param x (POEvaluation) A presence-only evaluation object to be messaged.
#' It could be the return of function `evaluate_po`.
#' @param ... Not used.
#' @importFrom stringr str_pad
#' @return The same object that was passed as input.
#' @export
#' @examples
#' print(evaluation_po)
#'
print.POEvaluation <- function(x, ...){
  # CVI
  cvi25 <- x$cvi$`cvi with 0.25`
  cvi05 <- x$cvi$`cvi with 0.5`
  cvi75 <- x$cvi$`cvi with 0.75`

  # CBI with 100 moving windows
  cbi <- x$boyce$Spearman.cor

  # AUC_ratio
  auc_r <- x$roc_ratio$auc_ratio

  # AUC_bg
  auc_bg <- x$roc_background$auc_background

  # print
  cat('Presence-only evaluation\n')
  cat(paste0(str_pad('CVI with 0.25 threshold:', 30, side = 'right', pad = ' '),
             sprintf('%.3f\n', cvi25)))
  cat(paste0(str_pad('CVI with 0.5 threshold:', 30, side = 'right', pad = ' '),
             sprintf('%.3f\n', cvi05)))
  cat(paste0(str_pad('CVI with 0.75 threshold:', 30, side = 'right', pad = ' '),
             sprintf('%.3f\n', cvi75)))
  cat(paste0(str_pad('CBI:', 30, side = 'right', pad = ' '),
             sprintf('%.3f\n', cbi)))
  cat(paste0(str_pad('AUC (ratio)',
                     30, side = 'right', pad = ' '),
             sprintf('%.3f\n', auc_r)))
  cat(paste0(str_pad('AUC (Presence-background)',
                     30, side = 'right', pad = ' '),
             sprintf('%.3f\n', auc_bg)))

  # return
  invisible(return(x))
}

#' @title Print summary information from ReducedImageStack object.
#' @description Display the most general and informative characteristics of
#' a ReducedImageStack object.
#' @param x (ReducedImageStack) A ReducedImageStack object to be messaged.
#' It could be the return of function `dim_reduce`.
#' @param ... Not used.
#' @importFrom stars st_get_dimension_values
#' @return The same object that was passed as input.
#' @export
#' @examples
#' print(img_reduced)
#'

print.ReducedImageStack <- function(x, ...) {
  cat('Dimension reduction\n')
  cat(sprintf('Correlation threshold: %s\n', x$threshold))
  cat(sprintf(
    'Original variables: %s\n',
    paste0(names(x$cors_original$mean), collapse = ', ')))
  cat(sprintf(
    'Variables after dimension reduction: %s\n',
    paste0(st_get_dimension_values(x$img_reduced, 'band'), collapse = ', ')))
  cat(paste0(paste(rep('=', 80), collapse = ''), "\n"))
  cat('Reduced correlations:\n')
  print(x$cors_reduced)

  # Return
  invisible(x)
}

#' @title Print summary information from PAConversion object.
#' @description Display the most general and informative characteristics of
#' a PAConversion object.
#' @param x (PAConversion) A PAConversion object to be messaged.
#' It could be the return of function `convert_to_pa`.
#' @param ... Not used.
#' @return The same object that was passed as input.
#' @export
#' @examples
#' print(pa_covert)
#'
print.PAConversion <- function(x, ...) {
  # threshold
  if (x$pa_conversion$conversion_method == 'threshold') {
    cat('Threshold conversion\n')
    cat(sprintf('cutoff = %s\n', x$pa_conversion$beta))
    cat(sprintf('species prevalence = %s\n', x$pa_conversion$species_prevalence))
  # logistic
  } else if (x$pa_conversion$conversion_method == 'logistic') {
    cat('Logistic conversion\n')
    cat(sprintf('beta = %s\n', x$pa_conversion$beta))
    cat(sprintf('alpha = %s\n', x$pa_conversion$alpha))
    cat(sprintf('species prevalence = %s\n', x$pa_conversion$species_prevalence))
  # linear
  } else if (x$pa_conversion$conversion_method == 'linear') {
    cat('Linear conversion\n')
    cat(sprintf('slope = %s\n', x$pa_conversion$a))
    cat(sprintf('intercept = %s\n', x$pa_conversion$b))
    cat(sprintf('species prevalence = %s\n', x$pa_conversion$species_prevalence))
  }

  # return
  invisible(x)
}

#' @title Print summary information from PAConversion object.
#' @description Display the most general and informative characteristics of
#' a PAConversion object.
#' @param x (EnvironmentalOutlier) A EnvironmentalOutlier object to be messaged.
#' It could be the return of function `suspicious_env_outliers`.
#' @param ... Not used.
#' @import outliertree
#' @return The same object that was passed as input.
#' @export
#' @examples
#' print(suspicious_outliers)
#'
print.EnvironmentalOutlier <- function(x, ...) {
  print(x$outlier_details)

  # Return
  invisible(x)
}
