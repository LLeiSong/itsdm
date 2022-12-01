# An internal function to plot response curves
#' @importFrom rlang .data
.plot_responses <- function(response_list,
                            smooth_span = 0.3){
  # Check inputs
  checkmate::assert_multi_class(
    response_list, c('MarginalResponse', 'IndependentResponse', 'List'))
  checkmate::assert_number(smooth_span)

  # Check categorical and continuous vars
  if (!is.null(response_list$responses_cont)) {
    responses_cont <- response_list$responses_cont

    # Convert list to tibble
    response_df <- do.call(
      rbind,
      lapply(1:length(responses_cont),
             function(n) {
               responses_cont[[n]] %>%
                 mutate(variable = names(responses_cont)[n])}))

    # Draw
    cex.axis <- 1
    cex.lab <- 1
    if (smooth_span == 0){
      invisible(g_cont <- ggplot(response_df, aes(x = .data$x, y = .data$y)) +
                  geom_line(color = "black") +
                  #scale_y_continuous(limits = c(0, 1.0)) +
                  facet_wrap(~variable, scales = 'free', ncol = 2) +
                  xlab('Value') + ylab('(Continuous variables)\nStandardized suitability') +
                  theme(axis.text = element_text(size = rel(cex.axis)),
                        axis.title = element_text(size = rel(cex.lab)),
                        plot.title = element_text(hjust = 0.5)) +
                  theme_linedraw())
    } else {
      invisible(g_cont <- ggplot(response_df, aes(x = .data$x, y = .data$y)) +
                  geom_point(alpha = 0) +
                  stat_smooth(method = 'loess', span = smooth_span, color = "black") +
                  #scale_y_continuous(limits = c(0, 1.0)) +
                  facet_wrap(~variable, scales = 'free', ncol = 2) +
                  xlab('Value') + ylab('(Continuous variables)\nStandardized suitability') +
                  theme(axis.text = element_text(size = rel(cex.axis)),
                        axis.title = element_text(size = rel(cex.lab)),
                        plot.title = element_text(hjust = 0.5)) +
                  theme_linedraw())
    }
  } else g_cont <- NULL

  if (!is.null(response_list$responses_cat)) {
    responses_cat <- response_list$responses_cat

    # Convert list to tibble
    response_df <- do.call(
      rbind,
      lapply(1:length(responses_cat),
             function(n) {
               responses_cat[[n]] %>%
                 mutate(variable = names(responses_cat)[n])}))
    var_order <- response_df %>%
      mutate(x = as.numeric(levels(.data$x))[.data$x]) %>%
      pull(.data$x) %>% unique() %>% sort()
    response_df <- response_df %>%
      mutate(x = factor(.data$x, levels = var_order))

    # Draw
    cex.axis <- 1
    cex.lab <- 1
    g_cat <- ggplot(response_df, aes(x = .data$x, y = .data$y)) +
      geom_bar(fill = "black",
               stat = 'identity',
               position = position_dodge()) +
      scale_y_continuous(limits = c(0, 1.0)) +
      facet_wrap(~variable, scales = 'free', ncol = 2) +
      xlab('Value') + ylab('(Categorical variables)\nStandardized suitability') +
      theme(axis.text = element_text(size = rel(cex.axis)),
            axis.title = element_text(size = rel(cex.lab)),
            plot.title = element_text(hjust = 0.5)) +
      theme_linedraw()
  } else g_cat <- NULL

  if (!is.null(g_cont) & !is.null(g_cat)){
    g_cont + g_cat +
      plot_layout(guides = "collect", nrow = 2, ncol = 1,
                  heights = c(ceiling(length(response_list$responses_cont) / 2),
                              ceiling(length(response_list$responses_cat) / 2)))
  } else if (!is.null(g_cont)) {
    g_cont
  } else if (!is.null(g_cat)) {
    g_cat
  }
}
# .plot_response end

#' @title Show marginal response curves.
#' @description Plot marginal response curves using ggplot2 by optionally set
#' target variable(s).
#' @param x (`MarginalResponse`) The marginal response curve object to plot.
#' It could be the return of function \code{\link{marginal_response}}.
#' @param target_var (`vector` of `character`) The target variable to plot. It could be
#' `NA`. If it is `NA`, all variables will be plotted.
#' @param smooth_span (`numeric`) The span value for smooth fit in `ggplot2`.
#' When it is `0`, no smooth applied. The default is `0.3`.
#' @param ... Not used.
#' @return `ggplot2` figure of response curves
#' @seealso
#' \code{\link{marginal_response}}
#'
#' @import ggplot2
#' @importFrom patchwork plot_layout
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' marginal_responses <- marginal_response(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables)
#' plot(marginal_responses, target_var = 'bio1')
#'}
#'
plot.MarginalResponse <- function(x,
                                  target_var = NA,
                                  smooth_span = 0.3,
                                  ...){
  # Checking
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    nms <- do.call(c, sapply(x, names))
    stopifnot(all(target_var %in% nms))}

  # Subset
  if (!all(is.na(target_var))) {
    cls <- class(x)
    x <- lapply(x, function(responses) {
      target_this <- intersect(target_var, names(responses))
      if (length(target_this) > 0) {
        responses[target_this]
      } else NULL})
    class(x) <- cls}

  # Plot
  .plot_responses(x, smooth_span)
}

#' @title Show independent response curves.
#' @description Plot independent response curves using ggplot2 by optionally
#' set target variable(s).
#' @param x (`IndependentResponse`) The independent response curve object to plot.
#' It could be the return of function \code{\link{independent_response}}.
#' @param target_var (`vector` of `character`) The target variable to plot. It could be
#' `NA`. If it is `NA`, all variables will be plotted.
#' @param smooth_span (`numeric`) The span value for smooth fit in `ggplot2`.
#' When it is `0`, no smooth applied. The default is `0.3`.
#' @param ... Not used.
#' @return `ggplot2` figure of response curves
#' @seealso
#' \code{\link{independent_response}}
#'
#' @import ggplot2
#' @importFrom patchwork plot_layout
#' @importFrom dplyr arrange slice
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' independent_responses <- independent_response(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables)
#' plot(independent_responses)
#'}
#'
plot.IndependentResponse <- function(x,
                                     target_var = NA,
                                     smooth_span = 0.3,
                                     ...){
  # Checking
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    nms <- do.call(c, sapply(x, names))
    stopifnot(all(target_var %in% nms))}

  # Subset
  if (!all(is.na(target_var))) {
    cls <- class(x)
    x <- lapply(x, function(responses) {
      target_this <- intersect(target_var, names(responses))
      if (length(target_this) > 0) {
        responses[target_this]
      } else NULL})
    class(x) <- cls}

  # Plot
  .plot_responses(x, smooth_span)
}

#' @title Show variable dependence plots and variable interaction plots
#' obtained from Shapley values.
#' @description Plot Shapley value-based variable dependence curves using
#' ggplot2 by optionally selecting target variable(s). It also can plot the
#' interaction between a related variable to the selected variable(s).
#' @param x (`ShapDependence`) The variable dependence object to plot.
#' It could be the return of function \code{\link{shap_dependence}}.
#' @param target_var (`vector` of `character`) The target variable to plot. It could be
#' `NA`. If it is `NA`, all variables will be plotted.
#' @param related_var (`character`) The dependent variable to plot together with
#' target variables. It could be `NA`. If it is `NA`, no related variable will be
#' plotted.
#' @param sample_prop (`numeric`) The proportion of points to sample for plotting.
#' It will be ignored if the number of points is less than 1000.
#' The default is `0.3`.
#' @param sample_bin (`integer`) The number of bins to use for stratified sampling.
#' @param smooth_line (`logical`) Whether to fit the smooth line or not.
#' It will be ignored if the number of points is less than 1000.
#' The default is 100.
#' @param seed (`integer`) The seed for sampling.
#' It will be ignored if the number of points is less than 1000.
#' The default is 123.
#' @param ... Other arguments passed on to \code{\link{geom_smooth}}. Mainly
#' `method` and `formula` to fit the smooth line. Note that the same arguments
#' will be used for all target variables. User could set variable one by one to
#' set the arguments separately.
#' @return `ggplot2` figure of dependent curves
#' @seealso
#' \code{\link{shap_dependence}}
#'
#' @details
#' If the number of samples is more than 1000, a stratified sampling is used to
#' thin the sample pool, and then plot its subset. The user could set a proportion
#' to sample and a number of bins for stratified sampling.
#'
#' @import ggplot2
#' @importFrom methods is
#' @importFrom patchwork plot_layout
#' @importFrom dplyr mutate mutate_if
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_dependence <- shap_dependence(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables)
#' plot(var_dependence, target_var = 'bio1', related_var = 'bio12')
#'}
#'
plot.ShapDependence <- function(x,
                                target_var = NA,
                                related_var = NA,
                                sample_prop = 0.3,
                                sample_bin = 100,
                                smooth_line = TRUE,
                                seed = 123,
                                ...) {
  # Checking
  nms <- names(x$feature_values)
  checkmate::assert_string(related_var, na.ok = T)
  if (!is.na(related_var)) {
    stopifnot(related_var %in% nms)}
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    stopifnot(all(target_var %in% nms))}
  if (nrow(x$feature_values) > 1000) {
    thin <- TRUE
    checkmate::assert_number(
      sample_prop, lower = 0, upper = 1, na.ok = TRUE, null.ok = TRUE)
    if (is.na(sample_prop) | is.null(sample_prop)) sample_prop <- 0.3
    checkmate::assert_int(
      sample_bin, lower = 0, upper = nrow(x$feature_values) %/% 2,
      na.ok = TRUE, null.ok = TRUE)
    if (is.na(sample_bin) | is.null(sample_bin)) sample_bin <- 100
    checkmate::assert_int(seed, na.ok = TRUE, null.ok = TRUE)
    if (is.na(seed) | is.null(seed)) seed <- 123
  } else thin <- FALSE
  checkmate::assert_flag(smooth_line)

  # Subset x if target_var is not NA
  # feature_values is the X in explain, so keep it.
  bands_cat <- names(x$dependences_cat)
  if (!all(is.na(target_var))) {
    cls <- class(x)
    x <- lapply(x, function(deps) {
      if (is(deps,'data.frame')) {
        deps
      } else {
        target_this <- intersect(target_var, names(deps))
        if (length(target_this) > 0) {
          deps[target_this]
        } else NULL}})
    class(x) <- cls}

  # Do not have related var
  if(is.na(related_var)) {
    # Continuous vars
    if (!is.null(x$dependences_cont)) {
      # Transform data
      x_cont <- x$dependences_cont
      nms <- names(x_cont)
      x_trans <- do.call(rbind, lapply(1:c(length(x_cont)), function(n) {
        x_cont[[n]] %>% mutate(variable = nms[n])
      }))

      if (thin){
        x_trans <- .shap_plot_thin(x_trans, sample_prop, sample_bin, seed)}

      # Plot
      cex.axis <- 1
      cex.lab <- 1
      g_cont <- ggplot(x_trans, aes(x = .data$x, y = .data$y)) +
        geom_point(size = 0.8) +
        xlab('Variable values') +
        ylab("(Continuous variables)\nShapley value") +
        facet_wrap(~variable, scales = 'free', ncol = 2) +
        theme(axis.text = element_text(size = rel(cex.axis)),
              axis.title = element_text(size = rel(cex.lab)),
              plot.title = element_text(hjust = 0.5)) +
        theme_linedraw()
      if (smooth_line) {
        g_cont <- g_cont +
          geom_smooth(color = 'red', alpha = 0, ...)
      }
    } else g_cont <- NULL

    # Categorical vars
    ## TODO: what if categorical var has many classes
    if (!is.null(x$dependences_cat)) {
      # Transform data
      x_cat <- x$dependences_cat
      nms <- names(x_cat)
      x_trans <- do.call(rbind, lapply(1:c(length(x_cat)), function(n) {
        x_cat[[n]] %>% mutate(variable = nms[n])
      }))

      if (thin){
        x_trans <- .shap_plot_thin(
          x_trans, sample_prop, sample_bin, seed, TRUE)}

      # Plot
      cex.axis <- 1
      cex.lab <- 1
      g_cat <- ggplot(x_trans, aes(x = .data$x, y = .data$y)) +
        geom_violin(position = position_dodge()) +
        geom_jitter(size = 0.8, position = position_jitter(0.2)) +
        xlab('Variable values') +
        ylab("(Categorical variables)\nShapley value") +
        facet_wrap(~variable, scales = 'free', ncol = 2) +
        theme(axis.text = element_text(size = rel(cex.axis)),
              axis.title = element_text(size = rel(cex.lab)),
              plot.title = element_text(hjust = 0.5)) +
        theme_linedraw()
    } else g_cat <- NULL

    # Merge together
    if (!is.null(g_cont) & !is.null(g_cat)){
      g_cont + g_cat +
        plot_layout(guides = "collect", nrow = 2, ncol = 1,
                    heights = c(ceiling(length(x$dependences_cont) / 2),
                                ceiling(length(x$dependences_cat) / 2)))
    } else if (!is.null(g_cont)) {
      g_cont
    } else if (!is.null(g_cat)) {
      g_cat
    }

    # Have related var
    ## related var is continuous
    ## related var is categorical
  } else {
    # Continuous vars
    if (!is.null(x$dependences_cont)) {
      # Transform data
      x_cont <- x$dependences_cont
      nms <- names(x_cont)
      vars <- x$feature_values
      x_trans <- do.call(rbind, lapply(1:length(x_cont), function(n) {
        x_cont[[n]] %>% mutate(related_var = vars %>% pull(related_var),
                               variable = rep(nms[n], nrow(.)))
      }))

      if (thin){
        x_trans <- .shap_plot_thin(x_trans, sample_prop, sample_bin, seed)}

      # Plot
      cex.axis <- 1
      cex.lab <- 1
      if (related_var %in% bands_cat) {
        g_cont <- ggplot(x_trans,
                         aes(x = .data$x, y = .data$y,
                             color = .data$related_var)) +
          geom_point(size = 0.8) +
          xlab('Variable values') +
          ylab("(Continuous variables)\nShapley value") +
          scale_color_viridis_d(related_var) +
          facet_wrap(~variable, scales = 'free', ncol = 2) +
          theme(axis.text = element_text(size = rel(cex.axis)),
                axis.title = element_text(size = rel(cex.lab)),
                plot.title = element_text(hjust = 0.5)) +
          theme_linedraw()
      } else {
        g_cont <- ggplot(x_trans,
                         aes(x = .data$x, y = .data$y,
                             color = .data$related_var)) +
          geom_point(size = 0.8) +
          xlab('Variable values') +
          ylab("(Continuous variables)\nShapley value") +
          scale_color_viridis_c(related_var) +
          facet_wrap(~variable, scales = 'free', ncol = 2) +
          theme(axis.text = element_text(size = rel(cex.axis)),
                axis.title = element_text(size = rel(cex.lab)),
                plot.title = element_text(hjust = 0.5)) +
          theme_linedraw()
      }
      if (smooth_line) {
        g_cont <- g_cont +
          geom_smooth(color = 'red', alpha = 0, ...)
      }
    } else g_cont <- NULL

    # Categorical vars
    ## TODO: what if categorical var has many classes
    if (!is.null(x$dependences_cat)) {
      # Transform data
      x_cat <- x$dependences_cat
      nms <- names(x_cat)
      vars <- x$feature_values
      x_trans <- do.call(rbind, lapply(1:c(length(x_cat)), function(n) {
        x_cat[[n]] %>% mutate(related_var = vars %>% pull(related_var),
                              variable = rep(nms[n], nrow(.)))
      }))

      if (thin){
        x_trans <- .shap_plot_thin(
          x_trans, sample_prop, sample_bin, seed, TRUE)}

      # Plot
      cex.axis <- 1
      cex.lab <- 1
      if (related_var %in% bands_cat) {
        g_cat <- ggplot(x_trans,
                        aes(x = .data$x, y = .data$y)) +
          geom_violin(position = position_dodge()) +
          geom_jitter(aes(color = related_var),
                      size = 0.8, position = position_jitter(0.2),
                      show.legend = FALSE) +
          xlab('Variable values') +
          ylab("(Categorical variables)\nShapley value") +
          scale_color_viridis_d(related_var) +
          facet_wrap(~variable, scales = 'free', ncol = 2) +
          theme(axis.text = element_text(size = rel(cex.axis)),
                axis.title = element_text(size = rel(cex.lab)),
                plot.title = element_text(hjust = 0.5)) +
          theme_linedraw()
      } else {
        g_cat <- ggplot(x_trans,
                        aes(x = .data$x, y = .data$y)) +
          geom_violin(position = position_dodge()) +
          geom_jitter(aes(color = related_var),
                      size = 0.8, position = position_jitter(0.2),
                      show.legend = FALSE) +
          xlab('Variable values') +
          ylab("(Categorical variables)\nShapley value") +
          scale_color_viridis_c(related_var) +
          facet_wrap(~variable, scales = 'free', ncol = 2) +
          theme(axis.text = element_text(size = rel(cex.axis)),
                axis.title = element_text(size = rel(cex.lab)),
                plot.title = element_text(hjust = 0.5)) +
          theme_linedraw()
      }

    } else g_cat <- NULL

    # Merge together
    if (!is.null(g_cont) & !is.null(g_cat)){
      g_cont + g_cat +
        plot_layout(guides = "collect", nrow = 2, ncol = 1,
                    heights = c(ceiling(length(x$dependences_cont) / 2),
                                ceiling(length(x$dependences_cat) / 2)))
    } else if (!is.null(g_cont)) {
      g_cont
    } else if (!is.null(g_cat)) {
      g_cat
    }
  }
}

#' @title Display spatial variable dependence maps.
#' @description Plot spatial variable dependence maps using ggplot2 by
#' optionally setting target variable(s).
#' @param x (`SpatialResponse`) The spatial variable dependence object to plot.
#' It could be the return of function \code{\link{spatial_response}}.
#' @param target_var (`vector` of `character`) The target variable to plot.
#' It could be `NA`. If it is `NA`, all variables will be plotted.
#' @param ... Not used.
#' @return `ggplot2` figure of dependent maps
#' @seealso
#' \code{\link{spatial_response}}
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom stars st_set_dimensions
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' spatial_responses <- spatial_response(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables,
#'   shap_nsim = 10)
#' plot(spatial_responses)
#' plot(spatial_responses, target_var = 'bio1')
#'}
#'
plot.SpatialResponse <- function(x,
                                 target_var = NA,
                                 ...) {
  # Checking
  nms <- names(x$spatial_marginal_response)
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    stopifnot(all(target_var %in% nms))}

  # Subset
  if (!all(is.na(target_var))) {
    x$spatial_marginal_response <- x$spatial_marginal_response[target_var]
    x$spatial_independent_response <- x$spatial_independent_response[target_var]
    x$spatial_shap_dependence <- x$spatial_shap_dependence[target_var]}

  if (length(names(x$spatial_marginal_response)) > 1) {
    # marginal plot
    m_plot <- ggplot() +
      geom_stars(
        data = do.call(c, x$spatial_marginal_response) %>%
          merge(name = 'band') %>%
          st_set_dimensions('band',
                            values = names(x$spatial_marginal_response))) +
      scale_fill_distiller('Value', palette = "YlOrRd",
                           direction = 1,
                           na.value = "transparent") +
      ggtitle(sprintf('Marginal effect of %s', 'variables')) +
      coord_equal() +
      facet_wrap(~band) +
      theme_void() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            strip.text.x = element_text(size = 12, face = "bold.italic"))

    # independent plot
    i_plot <- ggplot() +
      geom_stars(
        data = do.call(c, x$spatial_independent_response) %>%
          merge(name = 'band') %>%
          st_set_dimensions('band',
                            values = names(x$spatial_marginal_response))) +
      scale_fill_distiller('Value', palette = "YlOrRd",
                           direction = 1,
                           na.value = "transparent") +
      ggtitle(sprintf('Independent effect of %s', 'variables')) +
      coord_equal() +
      facet_wrap(~band) +
      theme_void() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            strip.text.x = element_text(size = 12, face = "bold.italic"))

    # SHAP plot
    if (!is.null(x$spatial_shap_dependence)) {
      s_plot <- ggplot() +
        geom_stars(
          data = do.call(c, x$spatial_shap_dependence) %>%
            merge(name = 'band') %>%
            st_set_dimensions('band',
                              values = names(x$spatial_marginal_response))) +
        scale_fill_distiller('Value', palette = "RdYlBu",
                             na.value = "transparent") +
        ggtitle(sprintf('SHAP-based effect of %s', 'variables')) +
        coord_equal() +
        facet_wrap(~band) +
        theme_void() +
        theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
              strip.text.x = element_text(size = 12, face = "bold.italic"))

      m_plot / i_plot / s_plot
    } else {
      m_plot / i_plot
    }
  } else {
    # marginal plot
    m_plot <- ggplot() +
      geom_stars(
        data = do.call(c, x$spatial_marginal_response)) +
      scale_fill_distiller('Value', palette = "YlOrRd",
                           direction = 1,
                           na.value = "transparent") +
      ggtitle(sprintf('Marginal effect of %s',
                      names(x$spatial_marginal_response))) +
      coord_equal() +
      theme_void() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            strip.text.x = element_text(size = 12, face = "bold.italic"))

    # independent plot
    i_plot <- ggplot() +
      geom_stars(
        data = do.call(c, x$spatial_independent_response)) +
      scale_fill_distiller('Value', palette = "YlOrRd",
                           direction = 1,
                           na.value = "transparent") +
      ggtitle(sprintf('Independent effect of %s',
                      names(x$spatial_marginal_response))) +
      coord_equal() +
      theme_void() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            strip.text.x = element_text(size = 12, face = "bold.italic"))

    # SHAP plot
    if (!is.null(x$spatial_shap_dependence)) {
      s_plot <- ggplot() +
        geom_stars(
          data = do.call(c, x$spatial_shap_dependence)) +
        scale_fill_distiller('Value', palette = "RdYlBu",
                             na.value = "transparent") +
        ggtitle(sprintf('SHAP-based effect of %s',
                        names(x$spatial_marginal_response))) +
        coord_equal() +
        theme_void() +
        theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
              strip.text.x = element_text(size = 12, face = "bold.italic"))

      m_plot / i_plot / s_plot
    } else {
      m_plot / i_plot
    }
  }
}

#' @title Exhibit variable contribution for target observations.
#' @description Use ggplot2 to plot variable contribution for each target
#' observation separately or summarize the overall variable contribution across
#' all selected observations.
#' @param x (`VariableContribution`) The `VariableContribution` object to plot.
#' It could be the return of function \code{\link{variable_contrib}}.
#' @param plot_each_obs (`logical`) The option of plot type. If `TRUE`, it will
#' plot variable contribution for every observation. Otherwise, it will plot
#' variable contribution violin plot for all observations.
#' @param num_features (`integer`) A number of most important features to plot.
#' Just work if plot_each_obs is `TRUE`.
#' @param ... Not used.
#' @return `ggplot2` figure of Variable Contribution.
#' @seealso
#' \code{\link{variable_contrib}}
#'
#' @import ggplot2
#' @importFrom dplyr arrange slice
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_contribution <- variable_contrib(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   var_occ_analysis = mod$vars_train %>% slice(1:10))
#'
#' # Plot variable contribution to each observation
#' plot(var_contribution,
#'      plot_each_obs = TRUE,
#'      num_features = 3)
#'
#' # Plot the summarized contribution
#' plot(var_contribution)
#'}
#'
plot.VariableContribution <- function(x,
                                      plot_each_obs = FALSE,
                                      # Just work for plot_each_obs is TRUE
                                      num_features = 5,
                                      ...) {
  # Checking
  checkmate::assert_logical(plot_each_obs)
  if (isTRUE(plot_each_obs)) {
    checkmate::assert_int(num_features,
                          lower = 1, upper = ncol(x$shapley_values))
  }
  if (isTRUE(plot_each_obs) & nrow(x$shapley_values) > 16) {
    stop(paste0('Too many observations in VariableContribution to plot separately. \n',
             'Consider to use less observations in VariableContribution or set ',
             'plot_each_obs to FALSE.'))}
  shapley_values <- x$shapley_values
  feature_values <- x$feature_values %>%
    mutate_if(is.numeric, function(x) round(x, 2)) %>%
    mutate_if(is.factor, function(x) as.character(x))
  stopifnot(identical(names(shapley_values), names(feature_values)))
  if (plot_each_obs) {
    checkmate::assert_int(num_features,
                          lower = 1, upper = ncol(shapley_values))}


  # Plot
  if (plot_each_obs) {
    # Convert data
    values_cont <- do.call(rbind, lapply(1:nrow(shapley_values), function(n) {
      vals <- feature_values[n, ]
      vals <- paste0(names(vals), ' = ', as.vector(vals))
      data.frame(num_obs  = paste0('Obs No.', n),
                 variable = vals,
                 shapley_value = unlist(shapley_values[n, ])) %>%
        arrange(-abs(.data$shapley_value)) %>% slice(1:num_features)
    })); row.names(values_cont) <- NULL
    values_cont <- values_cont %>%
      mutate(num_obs = factor(
        .data$num_obs, levels = paste0('Obs No.', 1:nrow(shapley_values))))

    # Plot
    ggplot(values_cont,
           aes(x = .data$variable,
               y = .data$shapley_value)) +
      geom_bar(aes(fill = abs(.data$shapley_value)),
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
           aes(x = .data$variable,
               y = .data$shapley_value)) +
      geom_violin(position = position_dodge()) +
      geom_jitter(aes(color = abs(.data$shapley_value)),
                  size = 0.8, position = position_jitter(0.2)) +
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

#' @title Display variable importance.
#' @description Display informative and detailed figures of variable importance.
#' @param x (`VariableAnalysis`) The variable importance object to plot.
#' It could be the return of function \code{\link{variable_analysis}}.
#' @param ... Not used.
#' @return A `patchwork` of `ggplot2` figure of variable importance
#' according to multiple metrics.
#' @seealso
#' \code{\link{variable_analysis}}, \code{\link{print.VariableAnalysis}}
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr summarise across
#' @importFrom tidyselect all_of
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' var_analysis <- variable_analysis(
#'   model = mod$model,
#'   pts_occ = mod$observation,
#'   pts_occ_test = mod$independent_test,
#'   variables = mod$variables)
#' plot(var_analysis)
#'}
#'
plot.VariableAnalysis <- function(x, ...) {
  # Pearson correlation
  cor_x <- x$pearson_correlation

  # Training
  ## Rotate dataset
  var_order_train <- cor_x %>%
    filter(.data$usage == 'Train') %>%
    filter(.data$method == 'Only') %>%
    arrange(.data$value) %>% pull(.data$variable)
  cor_x_train <- cor_x %>%
    filter(.data$usage == 'Train') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Train', 2),
                 value = c(1, 0))) %>%
    mutate(variable = factor(
      .data$variable,
      levels = c('', var_order_train)),
      method = factor(
        .data$method,
        levels = c('full', 'Only', 'Without')))

  ## Plot
  g_cor_train <- ggplot(cor_x_train,
         aes(x = .data$variable,
             y = .data$value, fill = .data$method)) +
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
    filter(.data$usage == 'Test') %>%
    filter(.data$method == 'Only') %>%
    arrange(.data$value) %>% pull(.data$variable)
  cor_x_test <- cor_x %>%
    filter(.data$usage == 'Test') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Test', 2),
                 value = c(1, 0))) %>%
    mutate(variable = factor(
      .data$variable,
      levels = c('', var_order_test)),
      method = factor(
        .data$method,
        levels = c('full', 'Only', 'Without')))
  ## Plot
  g_cor_test <- ggplot(cor_x_test) +
    geom_bar(aes(x = .data$variable,
                 y = .data$value, fill = .data$method),
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
    filter(.data$usage == 'Train') %>%
    filter(.data$method == 'Only') %>%
    arrange(.data$value) %>% pull(.data$variable)
  auc_r_train <- auc_r %>%
    filter(.data$usage == 'Train') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Train', 2),
                 value = c(full_auc_r$full_auc_train, 0))) %>%
    mutate(variable = factor(
      .data$variable,
      levels = c('', var_order_train)),
      method = factor(
        .data$method,
        levels = c('full', 'Only', 'Without')))

  ## Plot
  g_auc_train <- ggplot(auc_r_train,
                    aes(x = .data$variable,
                        y = .data$value, fill = .data$method)) +
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
    filter(.data$usage == 'Test') %>%
    filter(.data$method == 'Only') %>%
    arrange(.data$value) %>% pull(.data$variable)
  auc_r_test <- auc_r %>%
    filter(.data$usage == 'Test') %>%
    rbind(tibble(variable = rep('', 2),
                 method = c('full', 'Only'),
                 usage = rep('Test', 2),
                 value = c(full_auc_r$full_auc_test, 0))) %>%
    mutate(variable = factor(
      .data$variable,
      levels = c('', var_order_test)),
      method = factor(
        .data$method,
        levels = c('full', 'Only', 'Without')))

  ## Plot
  g_auc_test <- ggplot(auc_r_test,
                        aes(x = .data$variable,
                            y = .data$value,
                            fill = .data$method)) +
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
                        aes(x = .data$variable,
                            y = .data$value)) +
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
                         aes(x = .data$variable,
                             y = .data$value)) +
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

#' @title Show model evaluation.
#' @description Display informative and detailed figures of continuous Boyce
#' index, AUC curves, and TSS curve.
#' @param x (`POEvaluation`) The presence-only evaluation object to plot.
#' It could be the return of function \code{\link{evaluate_po}}.
#' @param ... Not used.
#' @return A `patchwork` of `ggplot2` figure of AUC_ratio, AUC_background and CBI.
#' @seealso
#' \code{\link{evaluate_po}}, \code{\link{print.POEvaluation}}
#'
#' @import ggplot2
#' @import patchwork
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' eval_train <- evaluate_po(
#'   mod$model,
#'   occ_pred = mod$pred_train$prediction,
#'   var_pred = na.omit(as.vector(mod$prediction[[1]])))
#'
#' plot(eval_train)
#'}
#'
plot.POEvaluation <- function(x, ...) {
  cex.axis <- 1
  cex.lab <- 1
  # Presence-only
  ## ROC ratio
  po_eval <- x$po_evaluation
  roc_r <- po_eval$roc_ratio$roc_ratio
  p_roc_r <- ggplot(roc_r, aes(y = .data$presence, x = .data$cell)) +
    geom_line(aes(colour = "roc", linetype = 'roc'), size = 0.8) +
    geom_line(aes(y = .data$cell, x = .data$cell,
                  colour = "chance", linetype = "chance"),
              size = 0.8) +
    geom_text(x = 0.8, y = 0.1,
              label = sprintf("AUC[ratio]: %s",
                              round(po_eval$roc_ratio$auc_ratio, 3)),
              parse = TRUE) +
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

  ## CBI
  cbi_bins <- data.frame(f_ratio = po_eval$boyce$F.ratio,
                         hs = po_eval$boyce$HS)
  p_boy <- ggplot(cbi_bins, aes(y = .data$f_ratio,
                                x = .data$hs)) +
    geom_line(colour = "black", size = 0.8) +
    geom_hline(yintercept = 1, color = 'red', size = 0.8) +
    scale_x_continuous(n.breaks = 9) +
    geom_text(x = 0.2, y = max(cbi_bins$f_ratio),
              label = sprintf("CBI: %s", round(po_eval$boyce$cor, 3))) +
    labs(y = "Predicted-to-expected (P/E) ratio",
         x = "Suitability") +
    ggtitle("Continuous P/E curve") +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
          plot.title.position = 'panel',
          axis.text = element_text(size = rel(cex.axis)),
          axis.title = element_text(size = rel(cex.lab)))

  # Presence-background
  ## AUC background
  if (!is.null(x$pb_evaluation)){
    pb_eval <- x$pb_evaluation
    roc_bg <- pb_eval$roc$roc
    roc_bg <- data.frame(tpr = roc_bg$TPR,
                         fpr = roc_bg$FPR)
    p_roc_bg <- ggplot(roc_bg, aes(y = .data$tpr,
                                   x = .data$fpr)) +
      geom_line(aes(colour = "roc", linetype = 'roc'), size = 0.8) +
      geom_line(aes(y = .data$fpr, x = .data$fpr,
                    colour = "chance", linetype = "chance"),
                size = 0.8) +
      geom_text(x = 0.8, y = 0.1,
                label = sprintf("AUC: %s",
                                round(pb_eval$roc$auc, 3))) +
      ggtitle('Receiver-operator curve') +
      labs(y = "Sensitivity (TPR)",
           x = "1-Specificity (FPR)") +
      scale_color_manual(
        '',
        values = c('roc' = 'black', 'chance' = 'grey'),
        labels = c('Empirical ROC curve',
                   'Chance line')) +
      scale_linetype_manual(
        '',
        values = c('roc' = 'solid', 'chance' = 'dashed'),
        labels = c('Empirical ROC curve',
                   'Chance line')) +
      theme_minimal() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            plot.title.position = 'panel',
            axis.text = element_text(size = rel(cex.axis)),
            axis.title = element_text(size = rel(cex.lab)),
            legend.position = 'top')

    ## Threshold-based
    pb_eval <- x$pb_evaluation
    tss_pb <- data.frame(thred = pb_eval$TSS$cutoff,
                         tss = pb_eval$TSS$tss)
    p_thred_bg <- ggplot(tss_pb, aes(y = .data$tss,
                                   x = .data$thred)) +
      geom_line(aes(colour = "tss", linetype = 'tss'), size = 0.8) +
      geom_line(aes(y = .data$tss,
                    x = rep(pb_eval$TSS$`Recommended threshold`,
                            length(.data$tss)),
                    colour = "best", linetype = "best"),
                size = 0.8) +
      geom_text(x = 0.5, y = 0.1,
                label = sprintf("Threshold: %s",
                                round(pb_eval$TSS$`Recommended threshold`, 3))) +
      ggtitle('Threshold-performance curve') +
      labs(y = "True skill statistic (TSS)",
           x = "Threshold") +
      scale_color_manual(
        '',
        values = c('tss' = 'black', 'best' = 'grey'),
        labels = c('TSS curve',
                   'Recommended thresold')) +
      scale_linetype_manual(
        '',
        values = c('tss' = 'solid', 'best' = 'dashed'),
        labels = c('TSS curve',
                   'Recommended thresold')) +
      theme_minimal() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            plot.title.position = 'panel',
            axis.text = element_text(size = rel(cex.axis)),
            axis.title = element_text(size = rel(cex.lab)),
            legend.position = 'top')

    # Ensemble
    (p_roc_r | p_boy) / (p_roc_bg | p_thred_bg)
  } else {
    (p_roc_r | p_boy)
  }
}

#' @title Display results of conversion to presence-absence (PA).
#' @description Display raster of suitability, probability of occurrence,
#' presence-absence binary map from presence-absence (PA) conversion.
#' @param x (`PAConversion`) The `PAConversion` object to plot.
#' It could be the return of function \code{\link{convert_to_pa}}.
#' @param ... Not used.
#' @return A `patchwork` of `ggplot2` figure of suitability, probability of occurrence,
#' presence-absence binary map.
#' @seealso
#' \code{\link{convert_to_pa}}, \code{\link{print.PAConversion}}
#'
#' @importFrom dplyr as_tibble
#' @import ggplot2
#' @import patchwork
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
#'   slice('band', c(1, 5, 12, 16))
#'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' # Threshold conversion
#' pa_thred <- convert_to_pa(mod$prediction,
#'   method = 'threshold', beta = 0.5)
#' plot(pa_thred)
#'}
#'
plot.PAConversion <- function(x, ...) {
  g1 <- ggplot() +
    geom_raster(data = as_tibble(x$suitability),
                aes(x = .data$x, y = .data$y,
                    fill = .data$prediction)) +
    ggtitle('Suitability') +
    scale_fill_viridis_c('Value', na.value = "transparent") +
    coord_equal() +
    theme_classic() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5))
  g2 <- ggplot() +
    geom_raster(data = as_tibble(x$probability_of_occurrence),
                aes(x = .data$x, y = .data$y,
                    fill = .data$prediction)) +
    ggtitle('Probability of ocurrence') +
    scale_fill_viridis_c('Value', na.value = "transparent") +
    coord_equal() +
    theme_classic() +
    theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5))
  g3 <- ggplot() +
    geom_raster(data = as_tibble(x$pa_map),
                aes(x = .data$x, y = .data$y,
                    fill = .data$prediction)) +
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

#' @title Exhibit suspicious outliers in an observation dataset.
#' @description Display observations and potential outliers diagnosed by
#' function \code{\link{suspicious_env_outliers}} in a dataset.
#' @param x (`EnvironmentalOutlier`) The PAConversion object to plot.
#' It could be the return of function \code{\link{suspicious_env_outliers}}.
#' @param overlay_raster (`RasterLayer` or `stars`) The environmental raster to plot
#' together with points.
#' @param pts_alpha (`numeric`) The `alpha` used by \code{\link{geom_sf}} to show points.
#' @param ... Not used.
#' @return A `ggplot2` figure of outliers distribution among all observations.
#' @seealso
#' \code{\link{suspicious_env_outliers}}, \code{\link{print.EnvironmentalOutlier}}
#'
#' @import ggplot2
#' @importFrom stars st_as_stars
#' @importFrom methods is
#' @importFrom stars geom_stars
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
#'   z_outlier = 3.5, outliers_print = 4L)
#'
#' plot(occ_outliers)
#' plot(occ_outliers,
#'   overlay_raster = env_vars %>% slice('band', 1))
#'}
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

  if (!is.null(x$outliers)) {
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
  } else {
    if (is.null(overlay_raster)) {
      ggplot() +
        geom_sf(data = x$pts_occ, aes(color = 'Normal'), size = 0.8) +
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
        scale_color_manual('',
                           values = c('Normal' = 'blue', 'Outlier' = 'red')) +
        ggtitle('Environmental outliers') +
        theme_classic() +
        theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5))
    }
  }

}

#' @title Display the figure and map of the `EnviChange` object.
#' @description Show the response curve and the map of contribution change from
#' \code{\link{detect_envi_change}}.
#' @param x (`EnviChange`) A `EnviChange` object to be messaged.
#' It could be the return of function \code{\link{detect_envi_change}}.
#' @param ... Not used.
#' @return The same object that was passed as input.
#' @seealso
#' \code{\link{detect_envi_change}}
#'
#' @import patchwork
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
#' #'
#' # Prepare data
#' data("occ_virtual_species")
#' obs_df <- occ_virtual_species %>% filter(usage == "train")
#' eval_df <- occ_virtual_species %>% filter(usage == "eval")
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#' #'
#' # Format the observations
#' obs_train_eval <- format_observation(
#'   obs_df = obs_df, eval_df = eval_df,
#'   x_col = x_col, y_col = y_col, obs_col = obs_col,
#'   obs_type = "presence_only")
#' #'
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12))
#' #'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 10,
#'   sample_size = 0.8, ndim = 1L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' # Use a fixed value
#' bio1_changes <- detect_envi_change(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables,
#'   shap_nsim = 1,
#'   target_var = "bio1",
#'   var_future = 5)
#'
#' plot(bio1_changes)
#'}
#'
plot.EnviChange <- function(x, ...) {
  x$p_curve + x$p_map
}

#' @title Display Shapley values-based spatial variable dependence maps.
#' @description Plot Shapley values-based spatial variable dependence maps
#' using ggplot2 by optionally setting target variable(s). This only works for
#' `SHAPSpatial` even though it is part of `SpatialResponse`.
#' @param x (`SHAPSpatial`) The spatial variable dependence object to plot.
#' It could be the return of function \code{\link{shap_spatial_response}}.
#' @param target_var (`vector` of `character`) The target variable to plot.
#' It could be `NA`. If it is `NA`, all variables will be plotted.
#' @param ... Not used.
#' @return `ggplot2` figure of dependent maps
#' @seealso
#' \code{\link{spatial_response}}
#'
#' @import ggplot2
#' @importFrom stars st_set_dimensions
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
#'   variables = env_vars, ntrees = 20,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' shap_spatial <- shap_spatial_response(
#'  model = mod$model,
#'  target_vars = c("bio1", "bio12"),
#'  var_occ = mod$vars_train,
#'  variables = mod$variables,
#'  shap_nsim = 1)
#'
#' plot(shap_spatial)
#' plot(shap_spatial, target_var = "bio1")
#'}
#'
plot.SHAPSpatial <- function(x,
                             target_var = NA,
                             ...) {
  # Checking
  nms <- names(x)
  checkmate::assert_character(target_var, min.len = 0)
  if (!all(is.na(target_var))) {
    stopifnot(all(target_var %in% nms))}

  # Subset
  if (!all(is.na(target_var))) x <- x[target_var]

  if (length(names(x)) > 1) {
    # SHAP plot
    ggplot() +
      geom_stars(
        data = do.call(c, x) %>%
          merge(name = 'band') %>%
          st_set_dimensions('band',
                            values = names(x))) +
      scale_fill_distiller('Value', palette = "RdYlBu",
                           na.value = "transparent") +
      ggtitle(sprintf('SHAP-based effect of %s', 'variables')) +
      coord_equal() +
      facet_wrap(~band) +
      theme_void() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            strip.text.x = element_text(size = 12, face = "bold.italic"))
  } else {
    # SHAP plot
    ggplot() +
      geom_stars(
        data = do.call(c, x)) +
      scale_fill_distiller('Value', palette = "RdYlBu",
                           na.value = "transparent") +
      ggtitle(sprintf('SHAP-based effect of %s',
                      names(x))) +
      coord_equal() +
      theme_void() +
      theme(plot.title = element_text(face = 'bold.italic', hjust = 0.5),
            strip.text.x = element_text(size = 12, face = "bold.italic"))
  }
}
# end of plot.SHAPSpatial
