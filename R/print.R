#' @title Print summary information from variable importance object.
#' @description Display the most general and informative characteristics of
#' a variable importance object.
#' @param x (`VariableAnalysis`) A variable importance object to be messaged.
#' It could be the return of function \code{\link{variable_analysis}}.
#' @param ... Not used.
#' @importFrom tidyselect all_of
#' @importFrom dplyr filter pull across summarise
#' @importFrom stringr str_pad
#' @importFrom rlang .data
#' @return The same object that was passed as input.
#' @seealso
#' \code{\link{variable_analysis}}, \code{\link{plot.VariableAnalysis}}
#'
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
#'   sample_size = 0.8, ndim = 3L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_analysis <- variable_analysis(
#'   model = mod$model,
#'   pts_occ = mod$pts_occ,
#'   pts_occ_test = mod$pts_occ_test,
#'   variables = mod$variables)
#'
#' print(variable_analysis)
#'}
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
    filter(.data$usage == 'Train') %>%
    arrange(-.data$value)

  n_max <- max(nchar(unique(cor_x_train$variable)))
  invisible(lapply(unique(cor_x_train %>%
                            filter(.data$method == 'Only') %>%
                            pull(.data$variable)), function(var){
    this_var <- cor_x_train %>% filter(.data$variable == var)
    only <- this_var %>% filter(.data$method == 'Only') %>% pull(.data$value)
    without <- this_var %>% filter(.data$method == 'Without') %>% pull(.data$value)

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
    filter(.data$usage == 'Test') %>%
    arrange(-.data$value)
  invisible(lapply(unique(cor_x_test %>%
                  filter(.data$method == 'Only') %>%
                  pull(.data$variable)), function(var){
    this_var <- cor_x_test %>% filter(.data$variable == var)
    only <- this_var %>% filter(.data$method == 'Only') %>% pull(.data$value)
    without <- this_var %>% filter(.data$method == 'Without') %>% pull(.data$value)

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
    filter(.data$usage == 'Train') %>%
    arrange(-.data$value)

  n_max <- max(nchar(unique(auc_x_train$variable)))
  invisible(lapply(unique(auc_x_train %>%
                            filter(.data$method == 'Only') %>%
                            pull(.data$variable)), function(var){
    this_var <- auc_x_train %>% filter(.data$variable == var)
    only <- this_var %>% filter(.data$method == 'Only') %>% pull(.data$value)
    without <- this_var %>% filter(.data$method == 'Without') %>% pull(.data$value)

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
    filter(.data$usage == 'Test') %>%
    arrange(-.data$value)
  invisible(lapply(unique(auc_x_test %>%
                            filter(.data$method == 'Only') %>%
                            pull(.data$variable)), function(var){
    this_var <- auc_x_test %>% filter(.data$variable == var)
    only <- this_var %>% filter(.data$method == 'Only') %>% pull(.data$value)
    without <- this_var %>% filter(.data$method == 'Without') %>% pull(.data$value)

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

  invisible(x)
}

#' @title Print summary information from presence-only evaluation object.
#' @description Display the most general and informative characteristics of
#' a presence-only evaluation object.
#' @param x (`POEvaluation`) A presence-only evaluation object to be messaged.
#' It could be the return of function \code{\link{evaluate_po}}.
#' @param ... Not used.
#' @importFrom stringr str_pad
#' @return The same object that was passed as input.
#' @seealso
#' \code{\link{evaluate_po}}, \code{\link{plot.POEvaluation}}
#'
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
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' eval_train <- evaluate_po(mod$model,
#'   occ_pred = mod$pred_train$prediction,
#'   var_pred = na.omit(as.vector(mod$prediction[[1]])))
#'
#' print(eval_train)
#'}
#'
print.POEvaluation <- function(x, ...){
  # Presence-only
  po_eval <- x$po_evaluation
  # CVI
  cvi25 <- po_eval$cvi$`cvi with 0.25`
  cvi05 <- po_eval$cvi$`cvi with 0.5`
  cvi75 <- po_eval$cvi$`cvi with 0.75`

  # CBI with 100 moving windows
  cbi <- po_eval$boyce$cor

  # AUC_ratio
  auc_r <- po_eval$roc_ratio$auc_ratio

  # print
  cat(paste0(paste(rep('=', 35), collapse = ''), '\n'))
  cat('Presence-only evaluation:\n')
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
  if (!is.null(x$pb_evaluation)) {
    cat(paste0(paste(rep('=', 35), collapse = ''), '\n'))
    cat('Presence-background evaluation:\n')
    pb_eval <- x$pb_evaluation
    cat(paste0(str_pad('Sensitivity:', 30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$sensitivity)))
    cat(paste0(str_pad('Specificity:', 30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$specificity)))
    cat(paste0(str_pad('TSS:', 30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$TSS$`Optimal TSS`)))
    cat(paste0(str_pad('AUC:', 30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$roc$auc)))
    cat('Similarity indices:\n')
    cat(paste0(str_pad("Jaccard's similarity index:",
                       30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$`Jaccard's similarity index`)))
    cat(paste0(str_pad("S\u00F8rensen's similarity index:",
                       30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$`f-measure`)))
    cat(paste0(str_pad("Overprediction rate:",
                       30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$`Overprediction rate`)))
    cat(paste0(str_pad("Underprediction rate:",
                       30, side = 'right', pad = ' '),
               sprintf('%.3f\n', pb_eval$`Underprediction rate`)))
  }

  # return
  invisible(return(x))
}

#' @title Print summary information from ReducedImageStack object.
#' @description Display the most general and informative characteristics of
#' a ReducedImageStack object.
#' @param x (`ReducedImageStack`) A `ReducedImageStack` object to be messaged.
#' It could be the return of function \code{\link{dim_reduce}}.
#' @param ... Not used.
#' @importFrom stars st_get_dimension_values
#' @return The same object that was passed as input.
#' @seealso
#' \code{\link{dim_reduce}}
#'
#' @export
#' @examples
#' \donttest{
#' library(itsdm)
#' library(dplyr)
#' library(stars)
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars()
#'
#' img_reduced <- dim_reduce(env_vars, threshold = 0.7,
#'   preferred_vars = c('bio1', 'bio12'))
#'
#' print(img_reduced)
#'}
#'
print.ReducedImageStack <- function(x, ...) {
  cat('Dimension reduction\n')
  cat(sprintf('Correlation threshold: %s\n', x$threshold))
  msg <- sprintf(
    'Original variables: %s\n',
    paste0(names(x$cors_original$mean), collapse = ', '))
  msg <- paste(strwrap(msg, width = 80), collapse = '\n')
  cat(paste0(msg, "\n"))
  msg <- sprintf(
    'Variables after dimension reduction: %s\n',
    paste0(st_get_dimension_values(x$img_reduced, 'band'), collapse = ', '))
  msg <- paste(strwrap(msg, width = 80), collapse = '\n')
  cat(paste0(msg, "\n"))
  cat(paste0(paste(rep('=', 80), collapse = ''), "\n"))
  cat('Reduced correlations:\n')
  print(round(x$cors_reduced, 2))

  # Return
  invisible(x)
}

#' @title Print summary information from PAConversion object.
#' @description Display the most general and informative characteristics of
#' a PAConversion object.
#' @param x (`PAConversion`) A PAConversion object to be messaged.
#' It could be the return of function \code{\link{convert_to_pa}}.
#' @param ... Not used.
#' @return The same object that was passed as input.
#' @seealso
#' \code{\link{convert_to_pa}}, \code{\link{plot.PAConversion}}
#'
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
#'   sample_size = 0.8, ndim = 1L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' # Threshold conversion
#' pa_thred <- convert_to_pa(mod$prediction, method = 'threshold', beta = 0.5)
#' print(pa_thred)
#'}
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
#' @param x (`EnvironmentalOutlier`) A `EnvironmentalOutlier` object to be messaged.
#' It could be the return of function \code{\link{suspicious_env_outliers}}.
#' @param ... Not used.
#' @import outliertree
#' @return The same object that was passed as input.
#' @seealso
#' \code{\link{suspicious_env_outliers}}, \code{\link{plot.EnvironmentalOutlier}}
#'
#' @export
#' @examples
#' \donttest{
#' library(dplyr)
#' library(sf)
#' library(stars)
#' library(itsdm)
#'
#' data("occ_virtual_species")
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12, 16))
#'
#' occ_outliers <- suspicious_env_outliers(
#'   occ = occ_virtual_species, variables = env_vars,
#'   z_outlier = 5, outliers_print = 4L)
#'
#' print(occ_outliers)
#'}
#'
print.EnvironmentalOutlier <- function(x, ...) {
  print(x$outlier_details)

  # Return
  invisible(x)
}

#' @title Print summary information from POIsotree object.
#' @description Display the most general and informative characteristics of
#' a fitted POIsotree object.
#' @param x (`POIsotree`) The POIsotree object to be messaged.
#' It could be the return of function \code{\link{isotree_po}}.
#' @param ... Not used.
#' @return The same object that was passed as input.
#' @seealso
#' \code{\link{isotree_po}}
#'
#' @importFrom stringr str_pad
#' @importFrom tidyselect all_of
#' @importFrom dplyr filter pull across summarise
#' @importFrom rlang .data
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
#'   sample_size = 0.8, ndim = 1L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE)
#' print(mod)
#'}
#'
print.POIsotree <- function(x, ...){
  # Model itself
  cat(paste0(paste(rep('=', 60), collapse = ''), '\n'))
  cat('Species distribution model:\n')
  print(x$model)
  cat(sprintf('Variables are: %s.\n',
              paste(names(x$variables), collapse = ', ')))
  cat(sprintf('Number of training records: %s.\n',
              nrow(x$pts_occ)))
  cat(sprintf('Use independent test? %s.\n',
              !is.null(x$pts_occ_test)))
  if (!is.null(x$pts_occ_test)) {
    cat(sprintf('Number of test records: %s.\n',
                nrow(x$pts_occ_test)))}
  cat(sprintf('Has marginal responses? %s.\n',
              !is.null(x$marginal_responses)))
  cat(sprintf('Has independent responses? %s.\n',
              !is.null(x$independent_responses)))
  cat(sprintf('Has Shapley value based responses? %s.\n',
              !is.null(x$shap_dependences)))
  cat(sprintf('Has spatial responses? %s.\n',
              !is.null(x$spatial_responses)))
  cat(sprintf('Has variable analysis? %s.\n',
              !is.null(x$variable_analysis)))


  # Evaluation
  cat(paste0(paste(rep('=', 60), collapse = ''), '\n'))
  cat('Model evaluation:\n')
  cat('[Training dataset]:\n')
  print(x$eval_train)
  if (!is.null(x$eval_test)) {
    cat('[Test dataset]:\n')
    print(x$eval_test)
  }

  # Variable importance
  if (!is.null(x$variable_analysis)) {
    cat(paste0(paste(rep('=', 60), collapse = ''), '\n'))
    cat('Variable importance:\n')
    cat('Just show SHAP result, use print or plot of variable_analysis to see all.\n')

    # SHAP
    ## Subset
    cat('SHAP (mean(|Shapley value|))\n')
    shap_x <- x$variable_analysis$SHAP
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

  # Return
  invisible(x)
}
