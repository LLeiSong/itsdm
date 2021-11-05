#' A function to evaluate the model based on presence-only data.
#' @description This function will calculate Contrast Validation Index (CVI),
#' continuous Boyce index (CBI) and AUC_ratio.
#' CVI is the proportion of presence points falling in cells having
#' a threshold (0.5 here) habitat suitability index minus the proportion of
#' cells within this range of threshold of the model.
#' Here we used varied thresholds: 0.25, 0.5, and 0.75.
#' AUC_background uses background values as pseudo-absence values
#' ROC curve plots the proportion of presences falling above a
#' range of thresholds against the proportion of cells falling
#' above the range of thresholds. The area under the modified
#' ROC curve was then called AUC_ratio.
#' @param model (isolation_forest) The extended isolation forest SDM.
#' @param occ_pred (vector) the vector contains predicted values at occurrence locations.
#' @param bg_pred (vector) the vector contains predicted values with same number
#' of background points.
#' @param var_pred (vector) the vector contains predicted values of the whole area.
#' @param visualize (logical) If TRUE, print message and plot the curves.
#' The default is FALSE.
#' @return (POEvaluation) a list of evaluations
#' @importFrom dplyr arrange tibble
#' @importFrom ecospat ecospat.boyce
#' @importFrom ROCit rocit ciAUC
#' @export
#' @examples
#' evaluate_po(mod, occ_pred = occ_pred, var_pred = var_pred)
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
