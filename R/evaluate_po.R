#' @title Evaluate the model based on presence-only data.
#' @description This function will calculate Contrast Validation Index (CVI),
#' continuous Boyce index (CBI), ROC_ratio and ROC_background.
#' @param model (`isolation_forest`) The extended isolation forest SDM. It could be
#' the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' @param occ_pred (`vector` of `numeric`) A `vector` contains predicted values
#' at occurrence locations.
#' @param bg_pred (`vector` of `numeric`) the vector contains predicted values
#' with same number of background points.
#' @param var_pred (`vector` of `numeric`) the vector contains predicted values
#' of the whole area.
#' @param visualize (`logical`) If `TRUE`, plot the evaluation figures.
#' The default is `FALSE`.
#' @return (`POEvaluation`) A list of
#' \itemize{
#' \item{cvi (`list`) A list of CVI with 0.25, 0.5, and 0.75 as threshold}
#' \item{boyce (`list`) A list of items related to continuous Boyce index (CBI)}
#' \item{roc_ratio (`list`) A list of ROC ratio and AUC ratio}
#' \item{roc_background (`list`) A list of ROC background and AUC background}}
#'
#' @seealso
#' \code{\link{print.POEvaluation}}, \code{\link{plot.POEvaluation}}
#'
#' @details
#' \bold{CVI} is the proportion of presence points falling in cells having
#' a threshold (`0.5` for example) habitat suitability index minus the proportion of
#' cells within this range of threshold of the model.
#' Here we used varied thresholds: `0.25`, `0.5`, and `0.75`.
#' \bold{ROC_background} uses background values as pseudo-absence values to generate
#' ROC curve. \bold{AUC_background} is the AUC calculated according to ROC_background.
#' \bold{ROC_ratio} curve plots the proportion of presences falling above a
#' range of thresholds against the proportion of cells falling
#' above the range of thresholds. The area under the modified
#' ROC curve was then called \bold{AUC_ratio}.
#'
#' @references
#' \itemize{
#' \item{\href{https://doi.org/10.1016/j.ecolmodel.2007.11.008}{Peterson,
#' A. Townsend, Monica Papeş, and Jorge Soberón. "Rethinking receiver operating
#' characteristic analysis applications in ecological niche modeling."
#' \emph{Ecological modelling} 213.1 (2008): 63-72.}}
#' \item{\href{https://doi.org/10.1016/j.ecolmodel.2006.05.017}{Hirzel,
#' Alexandre H., et al. "Evaluating the ability of habitat suitability models
#' to predict species presences." \emph{Ecological modelling}
#' 199.2 (2006): 142-152.}}
#' \item{\href{https://doi.org/10.1007/s00267-003-0040-3}{Hirzel, Alexandre
#' H., and Raphaël Arlettaz. "Modeling habitat suitability for complex species
#' distributions by environmental-distance geometric mean."
#' \emph{Environmental management} 32.5 (2003): 614-623.}}
#' }
#'
#' @importFrom dplyr arrange tibble
#' @importFrom ecospat ecospat.boyce
#' @importFrom ROCit rocit ciAUC
#' @export
#' @examples
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
#'   'extdata/bioclim_africa_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12, 16))
#'
#' mod <- isotree_po(
#'   occ = occ, occ_test = occ_test,
#'   variables = env_vars, ntrees = 200,
#'   sample_rate = 0.8, ndim = 2L,
#'   seed = 123L, response = FALSE,
#'   check_variable = FALSE)
#'
#' eval_train <- evaluate_po(mod$model,
#'   occ_pred = mod$pred_train$prediction,
#'   var_pred = na.omit(as.vector(mod$prediction[[1]])))
#'
evaluate_po <- function(model,
                        occ_pred,
                        bg_pred = NULL, # If NULL, skip AUC_bg
                        var_pred,
                        visualize = FALSE){
  # Check inputs
  checkmate::assert_vector(occ_pred)
  checkmate::assert_vector(var_pred)

  # CVIs
  ## CVI 0.25
  avi_test <- sum(occ_pred >= 0.25) / length(occ_pred)
  avi_all <- sum(var_pred >= 0.25) / length(var_pred)
  cvi25 <- avi_test - avi_all

  ## CVI 0.5
  avi_test <- sum(occ_pred >= 0.5)/length(occ_pred)
  avi_all <- sum(var_pred >= 0.5)/length(var_pred)
  cvi05 <- avi_test - avi_all

  ## CVI 0.75
  avi_test <- sum(occ_pred >= 0.75)/length(occ_pred)
  avi_all <- sum(var_pred >= 0.75)/length(var_pred)
  cvi75 <- avi_test - avi_all

  # CBI
  boy <- ecospat.boyce(fit = var_pred,
                       obs = occ_pred,
                       PEplot = F)
  cbi <- boy$Spearman.cor

  # AUC_ratio
  roc_r <- .roc_ratio(occ_pred, var_pred)
  auc_r <- .auc_ratio(occ_pred, var_pred)

  # AUC_background
  if (!is.null(bg_pred)) {
    occ_plus_bg <- tibble(prediction = c(occ_pred, bg_pred),
                          label = c(rep(1, length(occ_pred)),
                                    rep(0, length(bg_pred))))
    roc_b <- rocit(score = occ_plus_bg$prediction,
                   class = occ_plus_bg$label)
    auc_b <- ciAUC(roc_b)$AUC
  } else {
    message('Not set background samples, skip AUC_background.')
    roc_b <- NULL
    auc_b <- NULL
  }

  # Output
  out <- list(cvi = list(`cvi with 0.25` = cvi25,
                         `cvi with 0.5` = cvi05,
                         `cvi with 0.75` = cvi75),
              boyce = boy,
              roc_ratio = list(`roc_ratio` = roc_r,
                               `auc_ratio` = auc_r),
              roc_background = list(`roc_background` = roc_b,
                                    `auc_background` = auc_b))
  class(out) <- append("POEvaluation", class(out))

  # Visualize
  if (visualize) {
    plot(out)
  }

  # Return
  return(out)
}
# evaluate_po end
