# An internal function to plot response curves
.plot_responses <- function(response_list,
                            smooth_span = 0.3){
  # Check inputs
  checkmate::assert_multi_class(
    response_list, c('marginal_response', 'independent_response'))
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
#' @param x (marginal_response) The marginal response curve object to plot.
#' It could be the return of function `marginal_response`.
#' @param smooth_span (numeric) The span value for smooth fit in ggplot2.
#' When it is 0, no smooth applied. The default is 0.3.
#' @param ... Not used.
#' @return ggplot2 figure of response curves
#' @import ggplot2
#' @export
#' @examples
#' plot(marginal_responses)
#'
plot.marginal_response <- function(x, smooth_span = 0.3, ...){
  .plot_responses(x, smooth_span)
}

#' @title Function to plot independent response curves.
#' @description Plot independent response curves using ggplot2.
#' @param x (independent_response) The independent response curve object to plot.
#' It could be the return of function `independent_response`.
#' @param smooth_span (numeric) The span value for smooth fit in ggplot2.
#' When it is 0, no smooth applied. The default is 0.3.
#' @param ... Not used.
#' @return ggplot2 figure of response curves
#' @import ggplot2
#' @export
#' @examples
#' plot(independent_responses)
#'
plot.independent_response <- function(x, smooth_span = 0.3, ...){
  .plot_responses(x, smooth_span)
}

#' @title Function to plot variable importance.
#' @description Display informative and detailed figures of variable importance.
#' @param x (variable_analysis) The variable importance object to plot.
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
plot.variable_analysis <- function(x, ...) {
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
    mutate( variable = factor(
      variable,
      levels = c('', var_order_train)))

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
      levels = c('', var_order_test)))
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
    mutate( variable = factor(
      variable,
      levels = c('', var_order_train)))

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
    mutate( variable = factor(
      variable,
      levels = c('', var_order_test)))

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
  shap_x_train <- shap_x$train %>%
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
  shap_x_test <- shap_x$test %>%
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
#' @param x (evaluation_po) The presence-only evaluation object to plot.
#' It could be the return of function `evaluate_po`.
#' @param ... Not used.
#' @return a patchwork of ggplot figure of AUC_ratio, AUC_background and CBI.
#' @import ggplot2
#' @import patchwork
#' @export
#' @examples
#' plot(evaluation_po)
#'
plot.evaluation_po <- function(x, ...) {
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
