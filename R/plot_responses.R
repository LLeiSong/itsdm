#' A function to plot the response curves.
#' @description This function allows you to plot response curves using ggplot2.
#' @param response_list (marginal_response or independent_response)
#' The environmental_response object obtained from modeling.
#' @return ggplot2 figure of response curves
#' @import ggplot2
#' @export
#' @examples
#' plot_responses(response_list = responses)
#'
plot_responses <- function(response_list){
  # Check inputs
  checkmate::assert_multi_class(
    response_list, c('marginal_response', 'independent_response'))

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
  ggplot(response_df, aes(x = x, y = y)) +
    geom_line(colour = "black", size = 1) +
    scale_y_continuous(limits = c(0, 1.0)) +
    facet_grid(~variable, scales = 'free') +
    xlab('Value') + ylab('Suitability') +
    theme(axis.text = element_text(size = rel(cex.axis)),
          axis.title = element_text(size = rel(cex.lab)),
          plot.title = element_text(hjust = 0.5)) +
    theme_linedraw()
}

# plot_response end
