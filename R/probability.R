#' @title Estimate suitability on stars object using trained `isolation.forest` model.
#' @description Apply an `isolation.forest` model on a stars object to calculate
#' environmental suitability and do quantile stretch to `[0, 1]`.
#' @param x (`isolation_forest`). It could
#' be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' @param vars (`stars`) The stack of environmental variables. More specifically,
#' make sure it has x and y dimensions only, and distribute variables to
#' attributes of this `stars`. Otherwise, the function would stop.
#' @param offset (`numeric`) The offset to adjust fitted suitability. The default
#' is zero. Highly recommend to leave it as default.
#'
#' @return a `stars` of predicted habitat suitability
#' @seealso \code{\link{isotree_po}}
#' @import checkmate
#' @importFrom stats predict sd
#' @export
#' @examples
#' \dontrun{
#' # Using a pseudo presence-only occurrence dataset of
#' # virtual species provided in this package
#' library(dplyr)
#' library(sf)
#' library(stars)
#' library(itsdm)
#'
#' # Prepare data
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
#' # Modeling
#' mod_virtual_species <- isotree_po(occ = occ, occ_test = occ_test,
#'   variables = env_vars, ntrees = 50, sample_size = 0.8, ndim = 2L,
#'   seed = 123L)
#' suit <- probability(mod_virtual_species$model,
#'   mod_virtual_species$variables)
#'}
#'
probability <- function(x,
                        vars,
                        offset = 0) {
  # For later extension
  convert <- 'linear'
  pts_train <- NULL

  # Check inputs
  checkmate::assert_class(x, 'isolation_forest')
  checkmate::assert_class(vars, 'stars')
  checkmate::assert_choice(convert, c('linear', 'unify'))
  checkmate::assert_multi_class(pts_train, 'sf', null.ok = TRUE)
  if (!identical(names(dim(vars)), c('x', 'y'))) {
    stop(paste0('Please format inputs to stars object with ',
                'x and y dimensions only, and distribute variables to',
                ' attributes.'))}

  # Check data that fit x and vars
  ## model
  cols_cont <- x$metadata$cols_num
  cols_cat <- x$metadata$cols_cat

  ## vars
  bands <- names(vars)
  isfacor <- as.vector(sapply(vars, is.factor))
  bands_cont <- bands[!isfacor]
  bands_cat <- bands[isfacor]
  if (length(bands_cont) == 0) bands_cont <- NULL
  if (length(bands_cat) == 0) bands_cat <- NULL
  if (!identical(cols_cont, bands_cont) | !identical(cols_cat, bands_cat)) {
    stop('vars has different bands with data used to fit the model.')
  }; rm(cols_cat, cols_cont, bands, isfacor, bands_cont, bands_cat)

  # Check missing args
  if (convert == 'unify') {
    if (is.null(pts_train)) stop('Must set pts_train unify convertion.')
  }

  # Prediction
  ## Sample background
  var_pred <- predict(vars, x)

  if (convert == 'unify') {
    pred_train <- st_extract(var_pred, pts_train)$prediction
    # Convert score to probability using unify
    mu <- mean(pred_train)
    sigma <- sd(pred_train)
    pred_erf <- (var_pred - mu) / (sigma * sqrt(2))
    pred_erf <- .erf_stars(pred_erf)
    pred_erf[pred_erf < 0] <- 0
    pred_erf[pred_erf > 1] <- 1

    1 - pred_erf
    # linear conversion with a defined offset
  } else if (convert == 'linear') {
    var_pred <- 1 - var_pred
    .stars_stretch(var_pred, minv = 0, maxv = 1,
                   minq = offset, maxq = 1)
  }
}
