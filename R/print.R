#' @title Print summary information from variable importance object.
#' @description Display the most general and informative characteristics of
#' a variable importance object.
#' @param x An object to be messaged.
#' @param ... Not used.
#' @importFrom dplyr filter pull
#' @importFrom stringr str_pad
#' @return The same object that was passed as input.
#' @export
#' @examples
#' print(variable_importance)
print.variable_importance <- function(x, ...){
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
  .abs_mean <- function(v) mean(abs(v))
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
}

#' @export
print.marginal_response <- function(x, ...){

}

#' @export
print.independent_response <- function(x, ...){

}

#' @export
print.evaluation_po <- function(x, ...){

}
