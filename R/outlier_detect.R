#' A function to detect outliers from occurrence.
#' @description This function allows you to detect outliers from occurrence dataset.
#' @param occ (data.frame) The occurrence data frame for training.
#' There must be column x and y for coordinates. This argument is required
#' @param remove_outlier (logical) the option to remove detected outliers or not.
#' The default is FALSE.
#' @importFrom outliertree outlier.tree
#' @export
#' @examples
#' outlier_detect()
#'
outlier_detect <- function(env_vars,
                           remove_outlier = FALSE){
  # Check function arguments ----
  checkmate::assert_class(
    env_vars,
    c('data.frame', 'matrix', 'data.table', 'tibble'))

}

outlier_remove <- function(outlier_detecter){

}

# Big z_outlier value is suitable to detect wrong record over large area,
# such as global
# small z_outlier value is suitable to target critical region for a species
# (not very suitable, but have to be there)
