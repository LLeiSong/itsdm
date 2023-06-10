#' @title Evaluate the model based on presence-only data.
#' @description This function will calculate two major types of evaluation
#' metrics in terms of presence-only data. The first type is presence-only
#' customized metrics, such as Contrast Validation Index (CVI),
#' continuous Boyce index (CBI), and ROC_ratio.
#' The second type is presence-background evaluation metrics by
#' extracting background points as pseudo absence observations.
#' @param model (`isolation_forest`) The extended isolation forest SDM.
#' It could be the item `model` of `POIsotree` made by
#' function \code{\link{isotree_po}}.
#' @param occ_pred (`vector` of `numeric`) A `vector` contains predicted values
#' at occurrence locations.
#' @param bg_pred (`vector` of `numeric`) the vector contains predicted values
#' with same number of background points.
#' @param var_pred (`vector` of `numeric`) the vector contains predicted values
#' of the whole area. The reason to take a vector is to keep this function
#' flexible for multiple types of output.
#' @param threshold (`numeric` or `NULL`) The threshold to calculate
#' threshold-based evaluation metrics. If `NULL`, a recommended threshold
#' will be calculated based on optimal TSS value. The default is `NULL`.
#' @param visualize (`logical`) If `TRUE`, plot the evaluation figures.
#' The default is `FALSE`.
#' @return (`POEvaluation`) A list of
#' \itemize{
#' \item{
#' \bold{po_evaluation} is presence-only evaluation metrics. It is a list of
#' \itemize{
#' \item{cvi (`list`) A list of CVI with 0.25, 0.5, and 0.75 as threshold}
#' \item{boyce (`list`) A list of items related to continuous Boyce index (CBI)}
#' \item{roc_ratio (`list`) A list of ROC ratio and AUC ratio}
#' }}
#' \item{
#' \bold{pb_evaluation} is presence-background evaluation metrics.
#' It is a list of
#' \itemize{
#' \item{confusion matrix (`table`) A table of confusion matrix. The columns are
#' true values, and the rows are predicted values.}
#' \item{sensitivity (`numeric`) The sensitivity or TPR}
#' \item{specificity (`numeric`) The specificity or TNR}
#' \item{TSS (`list`) A list of info related to true skill statistic (TSS)
#' \itemize{
#' \item{cutoff (`vector` of `numeric`) A vector of cutoff threshold values}
#' \item{tss (`vector` of `numeric`) A vector of TSS for each cutoff threshold}
#' \item{Recommended threshold (`numeric`) A recommended threshold
#' according to TSS}
#' \item{Optimal TSS (`numeric`) The best TSS value}}}
#' \item{roc (`list`) A list of ROC values and AUC value}
#' \item{Jaccard's similarity index (`numeric`) The Jaccard's similarity index}
#' \item{Sørensen's similarity index (`numeric`) The Sørensen's similarity index
#' or F-measure}
#' \item{Overprediction rate (`numeric`) The Overprediction rate}
#' \item{Underprediction rate (`numeric`) The Underprediction rate}
#' }}}
#'
#' @seealso
#' \code{\link{print.POEvaluation}}, \code{\link{plot.POEvaluation}}
#'
#' @details
#' \itemize{
#' \item{
#' \bold{CVI} is the proportion of presence points falling in cells having
#' a threshold (`0.5` for example) habitat suitability index minus
#' the proportion of cells within this range of threshold of the model.
#' Here we used varied thresholds: `0.25`, `0.5`, and `0.75`.}
#' \item{\bold{continuous Boyce index (CBI)} is made with a 100 resolution of
#' moving windows and Kendall method.}
#' \item{
#' \bold{ROC_ratio} curve plots the proportion of presences falling above a
#' range of thresholds against the proportion of cells falling
#' above the range of thresholds. The area under the modified
#' ROC curve was then called \bold{AUC_ratio}.}
#' \item{\bold{Sensitivity (TPR)} = TP/(TP + FN)}
#' \item{\bold{Specificity (TNR)} = TN/(TN + FP)}
#' \item{\bold{True skill statistic (TSS)} = Sensitivity + specificity - 1}
#' \item{\bold{Jaccard's similarity index} = TP/(FN + TP + FP)}
#' \item{\bold{Sørensen's similarity index (F-measure)} = 2TP/(FN + 2TP + FP)}
#' \item{\bold{Overprediction rate} = FP/(TP + FP)}
#' \item{\bold{Underprediction rate} = FN/(TP + FN)}
#' }
#'
#' @references
#' \itemize{
#' \item{Peterson,
#' A. Townsend, Monica Papeş, and Jorge Soberón. "Rethinking receiver operating
#' characteristic analysis applications in ecological niche modeling."
#' \emph{Ecological modelling} 213.1 (2008): 63-72.
#' \doi{10.1016/j.ecolmodel.2007.11.008}}
#' \item{Hirzel,
#' Alexandre H., et al. "Evaluating the ability of habitat suitability models
#' to predict species presences." \emph{Ecological modelling}
#' 199.2 (2006): 142-152.
#' \doi{10.1016/j.ecolmodel.2006.05.017}}
#' \item{Hirzel, Alexandre
#' H., and Raphaël Arlettaz. "Modeling habitat suitability for complex species
#' distributions by environmental-distance geometric mean."
#' \emph{Environmental management} 32.5 (2003): 614-623.
#' \doi{10.1007/s00267-003-0040-3}}
#' \item{Leroy, Boris, et al.
#' "Without quality presence-absence data, discrimination metrics such as
#' TSS can be misleading measures of model performance."
#' \emph{Journal of Biogeography} 45.9 (2018): 1994-2002.
#' \doi{10.1111/jbi.13402}}
#' }
#'
#' @importFrom dplyr arrange tibble
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
#' # With perfect_presence mode,
#' # which should be very rare in reality.
#' mod <- isotree_po(
#'   obs_mode = "perfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 10,
#'   sample_size = 0.8, ndim = 2L,
#'   seed = 123L, nthreads = 1,
#'   response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' # Without background samples or absences
#' eval_train <- evaluate_po(
#'   mod$model,
#'   occ_pred = mod$pred_train$prediction,
#'   var_pred = na.omit(as.vector(mod$prediction[[1]])))
#' print(eval_train)
#'
#' # With background samples
#' bg_pred <- st_extract(
#'   mod$prediction, mod$background_samples) %>%
#'   st_drop_geometry()
#' eval_train <- evaluate_po(
#'   mod$model,
#'   occ_pred = mod$pred_train$prediction,
#'   bg_pred = bg_pred$prediction,
#'   var_pred = na.omit(as.vector(mod$prediction[[1]])))
#' plot(eval_train)
#' #'
evaluate_po <- function(model,
                        occ_pred,
                        # If NULL, skip presence-background
                        bg_pred = NULL,
                        var_pred,
                        threshold = NULL,
                        visualize = FALSE){
  # Check inputs
  checkmate::assert_vector(occ_pred)
  checkmate::assert_vector(var_pred)
  checkmate::assert_vector(bg_pred, null.ok = T)
  if (!is.null(bg_pred)){
    checkmate::assert_number(threshold, lower = 0,
                             upper = 1, null.ok = T)
  }
  checkmate::assert_logical(visualize)

  # Remove NAs
  occ_pred <- na.omit(occ_pred)
  bg_pred <- na.omit(bg_pred)
  var_pred <- na.omit(var_pred)

  ###############################################
  ############## Presence-only #################
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
  boy <- .cont_boyce(fit = var_pred,
                     obs = occ_pred,
                     method = 'kendall')

  # AUC_ratio
  roc_r <- .roc_ratio(occ_pred, var_pred)
  auc_r <- .auc_ratio(occ_pred, var_pred)

  po_eval <- list(cvi = list(`cvi with 0.25` = cvi25,
                             `cvi with 0.5` = cvi05,
                             `cvi with 0.75` = cvi75),
                  boyce = boy,
                  roc_ratio = list(`roc_ratio` = roc_r,
                                   `auc_ratio` = auc_r))

  ####################################################
  ############## Presence-background #################
  if (!is.null(bg_pred)) {
    # Table
    occ_plus_bg <- tibble(prediction = c(occ_pred, bg_pred),
                          label = c(rep(1, length(occ_pred)),
                                    rep(0, length(bg_pred))))

    # ROC and AUC_background
    roc_b <- rocit(score = occ_plus_bg$prediction,
                   class = occ_plus_bg$label)
    auc_b <- ciAUC(roc_b)$AUC

    # Threshold-based TSS curve and max TSS threshold
    cutoffs <- roc_b$Cutoff; cutoffs[1] <- 1
    TSSs <- sapply(roc_b$Cutoff, function(val) {
      binary_tb <- data.frame(
        prediction = occ_plus_bg$prediction >= val,
        label = occ_plus_bg$label)
      cm <- table(binary_tb$prediction, binary_tb$label)
      if (all(rownames(cm) == 'TRUE')) {
        cm <- rbind('FALSE' = c(0, 0), cm)
      } else if (all(rownames(cm) == 'FALSE')) {
        cm <- rbind(cm, 'TRUE' = c(0, 0))
      }
      sensitivity <- cm[4] / sum(cm[4], cm[3])
      specificity <- cm[1] / sum(cm[1], cm[2])
      sensitivity + specificity - 1
    })
    tss_df <- data.frame(Threshold = cutoffs,
                         TSS = TSSs)

    # Best or user-defined threshold
    if (is.null(threshold)){
      best_tss <- tss_df[tss_df$TSS == max(tss_df$TSS), ]
      if (nrow(best_tss) > 1) {
        threds <- best_tss$Threshold
        best_thred <- threds[abs(threds - 0.5) == min(abs(threds - 0.5))]
      } else best_thred <- best_tss$Threshold
      message(sprintf('Select %s as recommended threshold according to TSS.',
                      round(best_thred, 2)))
    } else {
      best_thred <- threshold
      message(sprintf('Set %s as threshold to calculate metrics.',
                      best_thred))
    }

    # Confusion matrix
    # TN (1) | FN (3)
    # FP (2) | TP (4)
    occ_plus_bg$prediction <- occ_plus_bg$prediction >= best_thred
    cm <- table(occ_plus_bg$prediction, occ_plus_bg$label)
    if (all(rownames(cm) == 'TRUE')) {
      cm <- rbind('FALSE' = c(0, 0), cm)
    } else if (all(rownames(cm) == 'FALSE')) {
      cm <- rbind(cm, 'TRUE' = c(0, 0))}

    # Traditional ones: sensitivity (TPR), specificity (TNR), and TSS
    # TPR  = TP / (TP + FN)
    # TNR = TN / (TN + FP)
    tpr <- cm[4] / sum(cm[4], cm[3])
    tnr <- cm[1] / sum(cm[1], cm[2])
    tss <- tpr + tnr - 1

    # Similarity indices of F-measures:
    # - Jaccard's similarity index = TP/(FN + TP + FP)
    # - Sørensen's similarity index, F-measure = 2TP/(FN + 2TP + FP)
    # - Overprediction rate = FP/(TP + FP)
    # - Underprediction rate = FN/(TP + FN)
    j_index <- cm[4] / sum(cm[3], cm[4], cm[2])
    s_index <- 2 * cm[4] / sum(cm[3], 2 * cm[4], cm[2])
    op_rate <- cm[2] / sum(cm[4], cm[2])
    up_rate <- 1 - tpr

    # Collect values
    pb_eval <- list(`confusion matrix` = cm,
                    sensitivity = tpr,
                    specificity = tnr,
                    TSS = list(cutoff = tss_df$Threshold,
                               tss = tss_df$TSS,
                               `Recommended threshold` = best_thred,
                               `Optimal TSS` = tss),
                    roc = list(roc = roc_b,
                               auc = auc_b),
                    `Jaccard's similarity index` = j_index,
                    # Sørensen's similarity index
                    `f-measure` = s_index,
                    `Overprediction rate` = op_rate,
                    `Underprediction rate` = up_rate)
  } else {
    message('Not set background samples, skip presence-background evaluation.')
    pb_eval <- NULL
  }

  # Output
  out <- list(po_evaluation = po_eval,
              pb_evaluation = pb_eval)
  class(out) <- append("POEvaluation", class(out))

  # Visualize
  if (visualize) {
    print(plot(out))
  }

  # Return
  return(out)
}
# evaluate_po end
