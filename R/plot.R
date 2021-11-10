# An internal function to plot response curves
.plot_responses <- function(response_list,
                            smooth_span = 0.3){
  # Check inputs
  checkmate::assert_multi_class(
    response_list, c('MarginalResponse', 'IndependentResponse', 'List'))
  checkmate::assert_number(smooth_span)

  # Convert list to tibble
  response_df <- do.call(
    rbind,
    lapply(1:length(response_list),
           function(n) {
             response_list[[n]] %>%
               mutate(variable = names(response_list)[n])}))

  # Draw
  cex.axis <- 1
  cex.lab <- 1
  if (smooth_span == 0){
    invisibel(g <- ggplot(response_df, aes(x = x, y = y)) +
      geom_line(color = "black") +
      scale_y_continuous(limits = c(0, 1.0)) +
      facet_wrap(~variable, scales = 'free', ncol = 2) +
      xlab('Value') + ylab('Standardized suitability') +
      theme(axis.text = element_text(size = rel(cex.axis)),
            axis.title = element_text(size = rel(cex.lab)),
            plot.title = element_text(hjust = 0.5)) +
      theme_linedraw())
  } else {
    invisible(g <- ggplot(response_df, aes(x = x, y = y)) +
      geom_point(alpha = 0) +
      stat_smooth(method = 'loess', span = smooth_span, color = "black") +
      scale_y_continuous(limits = c(0, 1.0)) +
      facet_wrap(~variable, scales = 'free', ncol = 2) +
      xlab('Value') + ylab('Standardized suitability') +
      theme(axis.text = element_text(size = rel(cex.axis)),
            axis.title = element_text(size = rel(cex.lab)),
            plot.title = element_text(hjust = 0.5)) +
      theme_linedraw())
  }
  suppressMessages(suppressWarnings(print(g)))
}
# .plot_response end

#' @title Function to plot marginal response curves.
#' @description Plot marginal response curves using ggplot2.
#' @param x (MarginalResponse) The marginal response curve object to plot.
#' It could be the return of function `marginal_response`.
#' @param target_var (vector of character) The target variable to plot. It could be
#' NA. If it is NA, all variables will be plotted.
#' @param smooth_span (numeric) The span value for smooth fit in ggplot2.
#' When it is 0, no smooth applied. The default is 0.3.
#' @param ... Not used.
#' @return ggplot2 figure of response curves
#' @import ggplot2
#' @export
#' @examples
#' plot(marginal_responses)
#'
plot.MarginalResponse <- function(x,
                                  target_var = NA,
                                  smooth_span = 0.3,
                                  ...){
  # Checking
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    stopifnot(all(target_var %in% names(x)))}
  if (!all(is.na(target_var))) {
    cls <- class(x)
    x <- x[target_var]
    class(x) <- cls}

  # Plot
  .plot_responses(x, smooth_span)
}

#' @title Function to plot independent response curves.
#' @description Plot independent response curves using ggplot2.
#' @param x (IndependentResponse) The independent response curve object to plot.
#' It could be the return of function `independent_response`.
#' @param target_var (vector of character) The target variable to plot. It could be
#' NA. If it is NA, all variables will be plotted.
#' @param smooth_span (numeric) The span value for smooth fit in ggplot2.
#' When it is 0, no smooth applied. The default is 0.3.
#' @param ... Not used.
#' @return ggplot2 figure of response curves
#' @import ggplot2
#' @importFrom dplyr arrange slice
#' @export
#' @examples
#' plot(independent_responses)
#'
plot.IndependentResponse <- function(x,
                                     target_var = NA,
                                     smooth_span = 0.3,
                                     ...){
  # Checking
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    stopifnot(all(target_var %in% names(x)))}
  if (!all(is.na(target_var))) {
    cls <- class(x)
    x <- x[target_var]
    class(x) <- cls}

  # Plot
  .plot_responses(x, smooth_span)
}

#' @title Function to plot variable dependence obtained from SHAP test.
#' @description Plot variable dependence curves using ggplot2.
#' @param x (VariableDependence) The variable dependence object to plot.
#' It could be the return of function `variable_dependence`.
#' @param target_var (vector of character) The target variable to plot. It could be
#' NA. If it is NA, all variables will be plotted.
#' @param related_var (character) The dependent variable to plot together with
#' target variables. It could be NA. If it is NA, no related variable will be
#' plotted.
#' @param smooth_span (numeric) The span value for smooth fit in ggplot2.
#' When it is 0, no smooth applied. The default is 0.3.
#' @param ... Not used.
#' @return ggplot2 figure of dependent curves
#' @import ggplot2
#' @importFrom dplyr mutate
#' @export
#' @examples
#' plot(var_dependence)
#'
plot.VariableDependence <- function(x,
                                    target_var = NA,
                                    related_var = NA,
                                    smooth_span = 0.3,
                                    ...) {
  # Checking
  checkmate::assert_string(related_var, na.ok = T)
  if (!is.na(related_var)) stopifnot(related_var %in% names(x))
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    stopifnot(all(target_var %in% names(x)))}

  # Subset x if target_var is not NA
  # feature_values is the X in explain, so keep it.
  if (!all(is.na(target_var))) {
    x <- x[c(target_var, 'feature_values')] }

  if(is.na(related_var)) {
    # Transform data
    nms <- names(x)
    x_trans <- do.call(rbind, lapply(1:c(length(x) - 1), function(n) {
      x[[n]] %>% mutate(variable = nms[n])
    }))

    # Plot
    cex.axis <- 1
    cex.lab <- 1
    p <- ggplot(x_trans, aes(x = x, y = y)) +
      geom_point(size = 0.8) +
      xlab('Variable values') +
      ylab("Shapley value") +
      facet_wrap(~variable, scales = 'free', ncol = 2) +
      theme(axis.text = element_text(size = rel(cex.axis)),
            axis.title = element_text(size = rel(cex.lab)),
            plot.title = element_text(hjust = 0.5)) +
      theme_linedraw()
    if (smooth_span > 0) {
      p <- p +
        geom_smooth(color = 'red', span = smooth_span, alpha = 0)
    }
  } else {
    # Transform data
    nms <- names(x)
    vars <- x$feature_values
    x_trans <- do.call(rbind, lapply(1:c(length(x) - 1), function(n) {
      x[[n]] %>% mutate(related_var = vars %>% pull(related_var),
                        variable = rep(nms[n], nrow(.)))
    }))

    # Plot
    cex.axis <- 1
    cex.lab <- 1
    p <- ggplot(x_trans,
                aes(x = x, y = y, color = related_var)) +
      geom_point(size = 0.8) +
      xlab('Variable values') +
      ylab("Shapley value") +
      scale_color_viridis_c(related_var) +
      facet_wrap(~variable, scales = 'free', ncol = 2) +
      theme(axis.text = element_text(size = rel(cex.axis)),
            axis.title = element_text(size = rel(cex.lab)),
            plot.title = element_text(hjust = 0.5)) +
      theme_linedraw()
    if (smooth_span > 0) {
      p <- p +
        geom_smooth(color = 'red', span = smooth_span, alpha = 0)
    }
  }

  p
}

#' @title Function to plot variable contribution for target observations.
#' @description Plot variable contribution for target observation separately
#' or together using ggplot2.
#' @param x (VariableContribution) The VariableContribution object to plot.
#' It could be the return of function `variable_contrib`.
#' @param plot_each_obs (logical) The option of plot type. If `TRUE`, it will
#' plot variable contribution for every observation. Otherwise, it will plot
#' variable contribution violin plot for all observations.
#' @param num_features (integer) A number of most important features to plot.
#' Just work if plot_each_obs is `TRUE`.
#' @param ... Not used.
#' @return ggplot2 figure of Variable Contribution.
#' @import ggplot2
#' @importFrom dplyr arrange slice
#' @export
#' @examples
#' plot(independent_responses)
#'
plot.VariableContribution <- function(x,
                                      plot_each_obs = FALSE,
                                      # Just work for plot_each_obs is TRUE
                                      num_features = 5,
                                      ...) {
  # Checking
  checkmate::assert_int(num_features, lower = 1, upper = nrow(x$shapley_values))
  checkmate::assert_logical(plot_each_obs)
  if (isTRUE(plot_each_obs) & nrow(x$shapley_values) > 16) {
    stop(paste0('Too many observations in VariableContribution to plot separately. \n',
             'Consider to use less observations in VariableContribution or set ',
             'plot_each_obs to FALSE.'))}
  shapley_values <- x$shapley_values
  feature_values <- x$feature_values
  stopifnot(identical(names(shapley_values), names(feature_values)))

  if (plot_each_obs) {
    # Convert data
    values_cont <- do.call(rbind, lapply(1:nrow(shapley_values), function(n) {
      vals <- round(feature_values[n, ], 2)
      vals <- paste0(names(vals), ' = ', as.vector(vals))
      data.frame(num_obs  = paste0('Obs No.', n),
                 variable = vals,
                 shapley_value = unlist(shapley_values[n, ])) %>%
        arrange(-abs(shapley_value)) %>% slice(1:num_features)
    })); row.names(values_cont) <- NULL
    values_cont <- values_cont %>%
      mutate(num_obs = factor(
        num_obs, levels = paste0('Obs No.', 1:nrow(shapley_values))))

    # Plot
    ggplot(values_cont,
           aes(x = variable,
               y = shapley_value)) +
      geom_bar(aes(fill = abs(shapley_value)),
               stat = 'identity',
               position = position_dodge()) +
      ggtitle('Variable contribution') +
      ylab('Shapley value') +
      xlab('') +
      scale_fill_viridis_c('Abs value') +
      facet_wrap(~num_obs, scales = 'free', ncol = 2) +
      theme_minimal() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            plot.title.position = 'panel') +
      coord_flip()
  } else {
    # Convert data
    values_cont <- do.call(rbind, lapply(1:nrow(shapley_values), function(n) {
      data.frame(variable = names(shapley_values),
                 shapley_value = unlist(shapley_values[n, ]))
    })); row.names(values_cont) <- NULL

    # Plot
    ggplot(values_cont,
           aes(x = variable,
               y = shapley_value)) +
      geom_violin(position = position_dodge()) +
      geom_jitter(aes(color = abs(shapley_value)),
                  size = 1.5, position = position_jitter(0.2)) +
      scale_color_viridis_c('Abs value') +
      ggtitle('Variable contribution') +
      ylab('Shapley value') +
      xlab('Environmental variables') +
      scale_fill_viridis_c('Abs value') +
      theme_minimal() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            plot.title.position = 'panel') +
      coord_flip()
  }
}

#' @title Function to plot variable importance.
#' @description Display informative and detailed figures of variable importance.
#' @param x (VariableAnalysis) The variable importance object to plot.
#' It could be the return of function `variable_analysis`.
#' @param ... Not used.
#' @return a patchwork of ggplot figure of variable importance
#' according to multiple metrics.
#' @import ggplot2
#' @import patchwork
#' @export
#' @examples
#' plot(variable_analysis)
#'
plot.VariableAnalysis <- function(x, ...) {
  # Pearson correlation
  cor_x <- x$pearson_correlation

  # Training
  ## Rotate dataset
  var_order_train <- cor_x %>%
    filter(usage == 'Train') %>%
    filter(method == 'Only') %>%
    arrange(value) %>% pull(variable)
  cor_x_train <- cor_x %>%
    filter(usage == 'Train') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Train', 2),
                 value = c(1, 0))) %>%
    mutate(variable = factor(
      variable,
      levels = c('', var_order_train)),
      method = factor(
        method,
        levels = c('full', 'Only', 'Without')))

  ## Plot
  g_cor_train <- ggplot(cor_x_train,
         aes(x = variable,
             y = value, fill = method)) +
    geom_bar(stat = 'identity',
             position = position_dodge()) +
    ggtitle('Jackknife of correlation on training data') +
    ylab('Pearson correlation with full model') +
    xlab('Environmental variable') +
    scale_fill_manual(
      'Training',
      values = c('red', 'black', 'lightgray'),
      labels = c('With all variables (= 1)', 'With only variable', 'Without variable')) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(0.5,"line")) +
    coord_flip()

  # Test
  ## Rotate dataset
  var_order_test <- cor_x %>%
    filter(usage == 'Test') %>%
    filter(method == 'Only') %>%
    arrange(value) %>% pull(variable)
  cor_x_test <- cor_x %>%
    filter(usage == 'Test') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Test', 2),
                 value = c(1, 0))) %>%
    mutate(variable = factor(
      variable,
      levels = c('', var_order_test)),
      method = factor(
        method,
        levels = c('full', 'Only', 'Without')))
  ## Plot
  g_cor_test <- ggplot(cor_x_test) +
    geom_bar(aes(x = variable,
                 y = value, fill = method),
             stat = 'identity',
             position = position_dodge()) +
    ggtitle('Jackknife of correlation on test data') +
    ylab('Pearson correlation with full model') +
    xlab('Environmental variable') +
    scale_fill_manual(
      'Test',
      values = c('red', 'black', 'lightgray'),
      labels = c('With all variables (= 1)', 'With only variable', 'Without variable')) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(0.5,"line")) +
    coord_flip()

  # AUC ratio
  auc_r <- x$AUC_ratio
  full_auc_r <- x$full_AUC_ratio

  # Training
  ## Rotate dataset
  var_order_train <- auc_r %>%
    filter(usage == 'Train') %>%
    filter(method == 'Only') %>%
    arrange(value) %>% pull(variable)
  auc_r_train <- auc_r %>%
    filter(usage == 'Train') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Train', 2),
                 value = c(full_auc_r$full_auc_train, 0))) %>%
    mutate(variable = factor(
      variable,
      levels = c('', var_order_train)),
      method = factor(
        method,
        levels = c('full', 'Only', 'Without')))

  ## Plot
  g_auc_train <- ggplot(auc_r_train,
                    aes(x = variable,
                        y = value, fill = method)) +
    geom_bar(stat = 'identity',
             position = position_dodge()) +
    ggtitle('Jackknife of AUC on training data') +
    ylab(expression('AUC'['ratio'])) +
    xlab('Environmental variable') +
    scale_fill_manual(
      'Training',
      values = c('red', 'black', 'lightgray'),
      labels = c(sprintf('With all variables (= %s)',
                         round(full_auc_r$full_auc_train, 2)),
                 'With only variable', 'Without variable')) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(0.5,"line")) +
    coord_flip()

  # Test
  ## Rotate dataset
  var_order_test <- auc_r %>%
    filter(usage == 'Test') %>%
    filter(method == 'Only') %>%
    arrange(value) %>% pull(variable)
  auc_r_test <- auc_r %>%
    filter(usage == 'Test') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Test', 2),
                 value = c(full_auc_r$full_auc_test, 0))) %>%
    mutate(variable = factor(
      variable,
      levels = c('', var_order_test)),
      method = factor(
        method,
        levels = c('full', 'Only', 'Without')))

  ## Plot
  g_auc_test <- ggplot(auc_r_test,
                        aes(x = variable,
                            y = value, fill = method)) +
    geom_bar(stat = 'identity',
             position = position_dodge()) +
    ggtitle('Jackknife of AUC on test data') +
    ylab(expression('AUC'['ratio'])) +
    xlab('Environmental variable') +
    scale_fill_manual(
      'Test',
      values = c('red', 'black', 'lightgray'),
      labels = c(sprintf('With all variables (= %s)',
                         round(full_auc_r$full_auc_test, 2)),
                         'With only variable', 'Without variable')) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          legend.title = element_blank(),
          legend.position = "top",
          legend.key.size = unit(0.5,"line")) +
    coord_flip()

  # SHAP
  ## Training
  shap_x_train <- x$SHAP$train %>%
    summarise(across(all_of(names(.)), .abs_mean))
  shap_x_train <- apply(shap_x_train, 1,
          function(x) sort(x)) %>%
    data.frame(value = .) %>%
    mutate(variable = factor(row.names(.),
                             levels = row.names(.)))
  g_shap_train <- ggplot(shap_x_train,
                        aes(x = variable,
                            y = value)) +
    geom_bar(stat = 'identity',
             fill = 'black',
             position = position_dodge()) +
    ggtitle('SHAP on training data') +
    ylab('mean(|Shapley value|)') +
    xlab('Environmental variable') +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel') +
    coord_flip()

  ## Test
  shap_x_test <- x$SHAP$test %>%
    summarise(across(all_of(names(.)), .abs_mean))
  shap_x_test <- apply(shap_x_test, 1,
                        function(x) sort(x)) %>%
    data.frame(value = .) %>%
    mutate(variable = factor(row.names(.),
                             levels = row.names(.)))
  g_shap_test <- ggplot(shap_x_test,
                         aes(x = variable,
                             y = value)) +
    geom_bar(stat = 'identity',
             fill = 'black',
             position = position_dodge()) +
    ggtitle('SHAP on test data') +
    ylab('mean(|Shapley value|)') +
    xlab('Environmental variable') +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel') +
    coord_flip()

  # Ensemble figures
  (g_cor_train | g_cor_test) /
    (g_auc_train | g_auc_test) /
    (g_shap_train | g_shap_test)
}

#' @title Function to plot presence-only evaluation.
#' @description Display informative and detailed figures of continuous Boyce
#' index and AUC curves.
#' @param x (POEvaluation) The presence-only evaluation object to plot.
#' It could be the return of function `evaluate_po`.
#' @param ... Not used.
#' @return a patchwork of ggplot figure of AUC_ratio, AUC_background and CBI.
#' @import ggplot2
#' @import patchwork
#' @export
#' @examples
#' plot(evaluation_po)
#'
plot.POEvaluation <- function(x, ...) {
  cex.axis <- 1
  cex.lab <- 1
  # ROC ratio
  roc_r <- x$roc_ratio$roc_ratio
  p_roc_r <- ggplot(roc_r, aes(y = presence,x = cell)) +
    geom_line(aes(colour = "roc", linetype = 'roc'), size = 0.8) +
    geom_line(aes(y = cell, x = cell,
                  colour = "chance", linetype = "chance"),
              size = 0.8) +
    geom_text(x = 0.8, y = 0.1,
              label = sprintf("AUC: %s",
                              round(x$roc_ratio$auc_ratio, 3))) +
    ggtitle('Modified ROC curve') +
    labs(y = "1 - omission error",
         x = "Proportion of area predicted present") +
    scale_color_manual(
      '',
      values = c('roc' = 'black', 'chance' = 'grey'),
      labels = c(expression('Empirical ROC'['ratio']~'curve'),
                 'Chance line')) +
    scale_linetype_manual(
      '',
      values = c('roc' = 'solid', 'chance' = 'dashed'),
      labels = c(expression('Empirical ROC'['ratio']~'curve'),
                 'Chance line')) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          axis.text = element_text(size = rel(cex.axis)),
          axis.title = element_text(size = rel(cex.lab)),
          legend.position = 'top')

  # AUC background
  roc_bg <- x$roc_background$roc_background
  roc_bg <- data.frame(tpr = roc_bg$TPR,
                       fpr = roc_bg$FPR)
  p_roc_bg <- ggplot(roc_bg, aes(y = tpr,
                                 x = fpr)) +
    geom_line(aes(colour = "roc", linetype = 'roc'), size = 0.8) +
    geom_line(aes(y = fpr, x = fpr,
                  colour = "chance", linetype = "chance"),
              size = 0.8) +
    geom_text(x = 0.8, y = 0.1,
              label = sprintf("AUC: %s",
                              round(x$roc_background$auc_background, 3))) +
    ggtitle('ROC curve') +
    labs(y = "Sensitivity (TPR)",
         x = "1-Specificity (FPR)") +
    scale_color_manual(
      '',
      values = c('roc' = 'black', 'chance' = 'grey'),
      labels = c(expression('Empirical ROC'['ratio']~'curve'),
                 'Chance line')) +
    scale_linetype_manual(
      '',
      values = c('roc' = 'solid', 'chance' = 'dashed'),
      labels = c(expression('Empirical ROC'['ratio']~'curve'),
                 'Chance line')) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          axis.text = element_text(size = rel(cex.axis)),
          axis.title = element_text(size = rel(cex.lab)),
          legend.position = 'top')

  # CBI
  cbi_bins <- data.frame(f_ratio = x$boyce$F.ratio,
                         hs = x$boyce$HS)
  p_boy <- ggplot(cbi_bins, aes(y = f_ratio,
                                x = hs)) +
    geom_line(colour = "black", size = 0.8) +
    geom_hline(yintercept = 1, color = 'red', size = 0.8) +
    scale_x_continuous(n.breaks = 9) +
    geom_text(x = 0.2, y = max(cbi_bins$f_ratio),
              label = sprintf("CBI: %s", round(x$boyce$Spearman.cor, 3))) +
    labs(y = "P/E ratio",
         x = "Suitability") +
    ggtitle("Continuous boyces index") +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          axis.text = element_text(size = rel(cex.axis)),
          axis.title = element_text(size = rel(cex.lab)))

  # Ensemble
  (p_roc_r | p_roc_bg) / p_boy
}

#' @title Function to plot results of conversion to PA.
#' @description Display raster of suitability, probability of occurrence,
#' presence-absence binary map from PA conversion.
#' @param x (PAConversion) The PAConversion object to plot.
#' It could be the return of function `convert_to_pa`.
#' @param ... Not used.
#' @importFrom dplyr as_tibble
#' @import ggplot2
#' @import patchwork
#' @return a patchwork of ggplot figure of suitability, probability of occurrence,
#' presence-absence binary map.
#' @export
#' @examples
#' plot(pa_covert)
#'
plot.PAConversion <- function(x, ...) {
  g1 <- ggplot() +
    geom_raster(data = as_tibble(x$suitability),
                aes(x = x, y = y, fill = prediction)) +
    ggtitle('Suitability') +
    scale_fill_viridis_c('Value', na.value = "transparent") +
    coord_equal() +
    theme_classic() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5))
  g2 <- ggplot() +
    geom_raster(data = as_tibble(x$probability_of_occurrence),
                aes(x = x, y = y, fill = prediction)) +
    ggtitle('Probability of ocurrence') +
    scale_fill_viridis_c('Value', na.value = "transparent") +
    coord_equal() +
    theme_classic() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5))
  g3 <- ggplot() +
    geom_raster(data = as_tibble(x$pa_map),
                aes(x = x, y = y, fill = prediction)) +
    ggtitle('Presence-absence') +
    scale_fill_manual(values = c('yellow', 'red', 'none'),
                      labels = c('Absence', 'Presence', ''),
                      na.value = "transparent") +
    coord_equal() +
    theme_classic() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          legend.title = element_blank())
  g1 / g2 / g3
}

#' @title Function to plot suspicious outliers in an observation dataset.
#' @description Display observations and outliers in a dataset.
#' @param x (EnvironmentalOutlier) The PAConversion object to plot.
#' It could be the return of function `suspicious_env_outliers`.
#' @param overlay_raster (RasterLayer or stars) The environmental raster to plot
#' together with points.
#' @param pts_alpha (numeric) the alpha used by ggplot2 to show points.
#' @param ... Not used.
#' @import ggplot2
#' @importFrom stars st_as_stars
#' @export
#' @examples
#' plot(suspicious_outliers)
#'
plot.EnvironmentalOutlier <- function(x,
                                      overlay_raster = NULL,
                                      pts_alpha = 0.5,
                                      ...) {
  # Check inputs
  checkmate::assert_multi_class(
    overlay_raster, c('RasterLayer', 'stars'), null.ok = T)
  # Convert overlay_raster if it is a raster
  if (is(overlay_raster, 'RasterLayer')){
    overlay_raster <- st_as_stars(overlay_raster)}
  checkmate::assert_number(pts_alpha, lower = 0, upper = 1)

  if (is.null(overlay_raster)) {
    ggplot() +
      geom_sf(data = x$pts_occ, aes(color = 'Normal'), size = 0.8) +
      geom_sf(data = x$outliers, aes(color = 'Outlier')) +
      scale_color_manual(values = c('Normal' = 'blue', 'Outlier' = 'red')) +
      xlab('x') + ylab('y') +
      ggtitle('Environmental outliers') +
      theme_classic() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            legend.title = element_blank())
  } else {
    ggplot() +
      geom_stars(data = overlay_raster) +
      scale_fill_viridis_c(names(overlay_raster)[1],
                           na.value = 'transparent') +
      geom_sf(data = x$pts_occ, aes(color = 'Normal'),
              size = 0.8, alpha = pts_alpha) +
      geom_sf(data = x$outliers, aes(color = 'Outlier')) +
      scale_color_manual('',
                         values = c('Normal' = 'blue', 'Outlier' = 'red')) +
      ggtitle('Environmental outliers') +
      theme_classic() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5))
  }
}
