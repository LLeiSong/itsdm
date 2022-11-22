#' @title Build Isolation forest species distribution model and explain the
#' the model and outputs.
#' @description Call Isolation forest and its variations to do
#' species distribution modeling and optionally call a collection of other
#' functions to do model explanation.
#' @param obs_mode (`string`) The mode of observations for training. It should
#' be one of `c("perfect_presence", "imperfect_presence", "presence_absence")`.
#' "perfect_presence" means presence-only occurrences without
#' errors/uncertainties/bias, which should be rare in reality.
#' "Imperfect_presence" means presence-only occurrences with
#' errors/uncertainties/bias, which should be a most common case.
#' "presence_absence" means presence-absence observations regardless quality.
#' See details to learn how to set it. The default is "imperfect_presence".
#' @param obs (`sf`) The `sf` of observation for training. It is recommended to
#' call function \code{\link{format_observation}} to format the
#' occurrence (`obs`) before passing it here.
#' Otherwise, make sure there is a column named "observation" for observation.
#' @param obs_ind_eval (`sf` or `NULL`) Optional `sf` of observations for
#' independent test. It is recommended to call function
#' \code{\link{format_observation}} to format the occurrence (`obs`)
#' before passing it here. Otherwise, make sure there is a column named
#' "observation" for observation.
#' If `NULL`, no independent test set will be used. The default is `NULL`.
#' @param variables (`RasterStack` or `stars`) The stack of environmental variables.
#' @param categ_vars (`vector` of `character` or `NULL`) The names of categorical
#' variables. Must be the same as the names in `variables`.
#' @param contamination (`numeric`) The percentage of abnormal cases within a
#' dataset. Because `iForest` is an outlier detection algorithm. It picks up
#' abnormal cases (much fewer) from normal cases. This argument is used to set
#' how many abnormal cases should be there if the users have the power to control.
#' See details for how to set it. The value should be less than 0.5. Here we
#' constrain it in (0, 0.3]. The default value is 0.1.
#' @param ntrees (`integer`) The number of trees for the isolation forest. It must
#' be integer, which you could use function \code{\link{as.integer}} to convert to.
#' The default is `100L`.
#' @param sample_size (`numeric`) It should be a rate for sampling size in `[0, 1]`.
#' The default is `1.0`.
#' @param ndim (`integer`) ExtensionLevel for isolation forest. It must
#' be integer, which you could use function \code{\link{as.integer}} to convert
#' to. Also, it must be no smaller than the dimension of environmental variables.
#' When it is 1, the model is a traditional isolation forest, otherwise the model
#' is an extended isolation forest. The default is 1.
#' @param seed (`integer`) The random seed used in the modeling. It should be an
#' integer. The default is `10L`.
#' @param ... Other arguments that \code{\link{isolation.forest}} needs.
#' @param offset (`numeric`) The offset to adjust fitted suitability. The default
#' is zero. Highly recommend to leave it as default.
#' @param response (`logical`) If `TRUE`, generate response curves.
#' The default is `TRUE`.
#' @param spatial_response (`logical`) If `TRUE`, generate spatial response maps.
#' The default is `TRUE` because it might be slow. NOTE that here SHAP-based map
#' is not generated because it is slow. If you want it be mapped, you could call
#' function \code{\link{spatial_response}} to make it.
#' @param check_variable (`logical`) If `TRUE`, check the variable importance.
#' The default is `TRUE`.
#' @param visualize (`logical`) If `TRUE`, generate the essential figures
#' related to the model. The default is `FALSE`.
#'
#' @return (`POIsotree`) A list of
#' \itemize{
#' \item{model (\code{\link{isolation.forest}}) The threshold set in
#' function inputs}
#' \item{variables (`stars`) The formatted image stack of
#' environmental variables}
#' \item{observation (\code{\link{sf}}) A \code{\link{sf}} of training occurrence
#' dataset}
#' \item{background_samples (\code{\link{sf}}) A \code{\link{sf}} of background points
#' for training dataset evaluation or SHAP dependence plot}
#' \item{independent_test (\code{\link{sf}} or `NULL`) A \code{\link{sf}} of test
#' occurrence dataset}
#' \item{background_samples_test (\code{\link{sf}} or `NULL`) A \code{\link{sf}} of
#' background points for test dataset evaluation or SHAP dependence plot}
#' \item{vars_train (\code{\link{data.frame}}) A \code{\link{data.frame}} with values of each
#' environmental variables for training occurrence}
#' \item{pred_train (\code{\link{data.frame}}) A \code{\link{data.frame}} with values of
#' prediction for training occurrence}
#' \item{eval_train (`POEvaluation`) A list of presence-only evaluation metrics
#' based on training dataset. See details of `POEvaluation` in
#' \code{\link{evaluate_po}}}
#' \item{var_test (\code{\link{data.frame}} or `NULL`) A \code{\link{data.frame}} with values of each
#' environmental variables for test occurrence}
#' \item{pred_test (\code{\link{data.frame}} or `NULL`) A \code{\link{data.frame}} with values of
#' prediction for test occurrence}
#' \item{eval_test (`POEvaluation` or `NULL`) A list of presence-only evaluation metrics
#' based on test dataset.
#' See details of `POEvaluation` in \code{\link{evaluate_po}}}
#' \item{prediction (`stars`) The predicted environmental suitability}
#' \item{marginal_responses (`MarginalResponse` or `NULL`) A list of marginal response
#' values of each environmental variables.
#' See details in \code{\link{marginal_response}}}
#' \item{offset (`numeric`) The offset value set as inputs.}
#' \item{independent_responses (`IndependentResponse` or `NULL`) A list of independent
#' response values of each environmental variables.
#' See details in \code{\link{independent_response}}}
#' \item{shap_dependences (`ShapDependence` or `NULL`) A list of variable
#' dependence values of each environmental variables.
#' See details in \code{\link{shap_dependence}}}
#' \item{spatial_responses (`SpatialResponse` or `NULL`) A list of spatial variable
#' dependence values of each environmental variables.
#' See details in \code{\link{shap_dependence}}}
#' \item{variable_analysis (`VariableAnalysis` or `NULL`) A list of variable importance
#' analysis based on multiple metrics.
#' See details in \code{\link{variable_analysis}}}}
#'
#' @seealso
#' \code{\link{evaluate_po}}, \code{\link{marginal_response}},
#' \code{\link{independent_response}}, \code{\link{shap_dependence}},
#' \code{\link{spatial_response}}, \code{\link{variable_analysis}},
#' \code{\link{isolation.forest}}
#'
#' @references
#' \itemize{
#' \item{Liu, Fei
#' Tony, Kai Ming Ting, and Zhi-Hua Zhou. "Isolation forest."
#' \emph{2008 eighth ieee international conference on data mining}.IEEE, 2008.
#' \doi{10.1109/ICDM.2008.17}}
#' \item{Liu, Fei Tony, Kai Ming
#' Ting, and Zhi-Hua Zhou. "Isolation-based anomaly detection."
#' \emph{ACM Transactions on Knowledge Discovery from Data (TKDD)} 6.1 (2012): 1-39.
#' \doi{10.1145/2133360.2133363}}
#' \item{Liu, Fei Tony,
#' Kai Ming Ting, and Zhi-Hua Zhou. "On detecting clustered anomalies using
#' SCiForest." \emph{Joint European Conference on Machine Learning and
#' Knowledge Discovery in Databases}. Springer, Berlin, Heidelberg, 2010.
#' \doi{10.1007/978-3-642-15883-4_18}}
#' \item{Ha
#' riri, Sahand, Matias Carrasco Kind, and Robert J. Brunner. "Extended
#' isolation forest." \emph{IEEE Transactions on Knowledge and Data Engineering (2019)}.
#' \doi{10.1109/TKDE.2019.2947676}}
#' \item{\url{https://github.com/david-cortes/isotree}}
#' \item{References of related feature such as response curves and variable importance
#' will be listed under their own functions}
#' }
#'
#' @details
#' For "perfect_presence", a user-defined number (`contamination`) of samples
#' will be taken from background to let `iForest` function normally.
#'
#' If "imperfect_presence", no further actions is required.
#'
#' If the \bold{obs_mode} is "presence_absence", a `contamination` percent
#' of absences will be randomly selected and work together with all presences
#' to train the model.
#'
#' NOTE: \bold{obs_mode} and \bold{mode} only works for `obs`. `obs_ind_eval`
#' will follow its own structure.
#'
#' Please read details of algorithm \code{\link{isolation.forest}} on
#' \url{https://github.com/david-cortes/isotree}, and
#' the R documentation of function \code{\link{isolation.forest}}.
#'
#' @import checkmate
#' @importFrom raster nlayers
#' @importFrom dplyr select slice mutate sample_n
#' @importFrom stars st_as_stars st_extract st_xy2sfc st_rasterize st_apply
#' @importFrom sf st_as_sf st_crs st_transform st_drop_geometry
#' @importFrom isotree isolation.forest
#' @importFrom rlang := .data
#' @importFrom stats na.omit predict
#' @importFrom methods is
#' @export
#' @examples
#' \donttest{
#' ########### Presence-absence mode #################
#' library(dplyr)
#' library(sf)
#' library(stars)
#' library(itsdm)
#'
#' # Load example dataset
#' data("occ_virtual_species")
#' obs_df <- occ_virtual_species %>% filter(usage == "train")
#' eval_df <- occ_virtual_species %>% filter(usage == "eval")
#' x_col <- "x"
#' y_col <- "y"
#' obs_col <- "observation"
#' obs_type <- "presence_absence"
#'
#' # Format the observations
#' obs_train_eval <- format_observation(
#'   obs_df = obs_df, eval_df = eval_df,
#'   x_col = x_col, y_col = y_col, obs_col = obs_col,
#'   obs_type = obs_type)
#'
#' # Load variables
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12))
#'
#' # Modeling
#' mod_virtual_species <- isotree_po(
#'   obs_mode = "presence_absence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 10,
#'   sample_size = 0.6, ndim = 1L,
#'   seed = 123L)
#'
#' # Check results
#' ## Evaluation based on training dataset
#' print(mod_virtual_species$eval_train)
#' plot(mod_virtual_species$eval_train)
#'
#' ## Response curves
#' plot(mod_virtual_species$marginal_responses)
#' plot(mod_virtual_species$independent_responses,
#'      target_var = c('bio1', 'bio5'))
#' plot(mod_virtual_species$shap_dependence)
#'
#' ## Relationships between target var and related var
#' plot(mod_virtual_species$shap_dependence,
#'      target_var = c('bio1', 'bio5'),
#'      related_var = 'bio12', smooth_span = 0)
#'
#' # Variable importance
#' mod_virtual_species$variable_analysis
#' plot(mod_virtual_species$variable_analysis)
#'
#' ########### Presence-absence mode ##################
#' # Load example dataset
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
#' # Modeling with perfect_presence mode
#' mod_perfect_pres <- isotree_po(
#'   obs_mode = "perfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 10,
#'   sample_size = 0.6, ndim = 1L,
#'   seed = 123L)
#'
#' # Modeling with imperfect_presence mode
#' mod_imperfect_pres <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 10,
#'   sample_size = 0.6, ndim = 1L,
#'   seed = 123L)
#' }
#'
isotree_po <- function(
  # SDM-related inputs
  # perfect_presence, imperfect_presence, presence_absence
  obs_mode = "imperfect_presence",
  obs, # sf
  obs_ind_eval = NULL, # sf
  variables, # rasterstack or stars. if stars, must have dimension band
  categ_vars = NULL, # categorical variables
  contamination = 0.1, # contamination level
  # Isotree-related inputs
  ntrees = 100L, # number of trees
  sample_size = 1.0, # sample size
  ndim = 1L, # extension level
  seed = 10L, # random seed, must be integer
  # Other arguments of isolation.forest
  ...,
  # Probability conversion
  offset = 0.0,
  # Other general inputs
  response = TRUE, # if generate response curves
  spatial_response = TRUE, # if generate spatial response maps
  check_variable = TRUE, # if generate var importance
  visualize = FALSE){ # if plot the curves

  # Check inputs - level 1
  checkmate::assert_choice(
    obs_mode,
    choices = c("perfect_presence", "imperfect_presence", "presence_absence"))
  checkmate::assert_multi_class(
    obs, c('sf', null.ok = F))
  checkmate::assert_multi_class(
    obs_ind_eval, c('sf', null.ok = T))
  if (!"observation" %in% names(obs)) {
    stop("No observation column in obs.")
  }

  if (!"observation" %in% names(obs_ind_eval)) {
    stop("No observation column in obs_ind_eval.")
  }

  checkmate::assert_multi_class(
    variables, c('RasterStack', 'stars'))
  checkmate::assert_vector(categ_vars, null.ok = T)
  if (obs_mode %in% c("perfect_presence", "presence_absence")) {
    checkmate::assert_number(
      contamination, lower = 0,
      upper = 0.3, na.ok = T)
  }
  checkmate::assert_int(ntrees)
  checkmate::assert_number(
    sample_size, lower = 0,
    upper = 1, na.ok = T)
  if (is(variables, 'RasterStack')) dim_max <- nlayers(variables)
  if (is(variables, 'stars')) {
    if (length(dim(variables)) == 2){
      dim_max <- length(variables)
    } else if (length(dim(variables)) == 3) {
      if (!is.null(categ_vars)){
        warning(paste0('Categorical layers are merged to an dimension.',
                       ' Be careful if they are original values.'))
      }
      dim_max <- length(st_get_dimension_values(variables, 3))
    } else {
      stop('variables has more than 3 dimensions, do not know which one to use.')
    }
  }
  checkmate::assert_int(
    ndim, lower = 1,
    upper = dim_max, na.ok = T)
  checkmate::assert_logical(response)
  checkmate::assert_logical(spatial_response)
  checkmate::assert_logical(check_variable)
  checkmate::assert_logical(visualize)

  # Check observations
  if (grepl("presence", obs_mode) &
      sum(obs$observation == 1) == 0) {
    stop(sprintf("Presences do not exist in observations for mode %s.",
                 obs_mode))}
  if (grepl("absence", obs_mode) &
      sum(obs$observation == 0) == 0) {
    stop(sprintf("Absences do not exist in observations for mode %s.",
                 obs_mode))}

  # Check categ_cols
  if (exists('categ_cols')) {
    stop('Set categ_vars instead.')
  }

  # Check inputs - level 2
  ## Check related columns
  if(!"observation" %in% colnames(obs)){
      stop("There must be observation column in occ.")}
  if((!is.null(obs_ind_eval)) &
     !("observation" %in% colnames(obs_ind_eval))){
      stop("There must be observation column in non-null occ_ind_eval.")}

  # Reformat the inputs
  # Variables
  ## -- RasterStack - resume cat original values, convert to stars with multiple
  ##                  attributes, check if categ_vars are existing cat vars.
  ## -- Stars with 3 dims - split it, no need to check cat vars
  ## -- Stars with multiple attributes - check categ_vars are existing cat vars
  if (is(variables, 'RasterStack')){
    # Check categ_vars
    ## step 1
    stopifnot(all(categ_vars %in% names(variables)))
    ## step 2
    if (!identical(names(variables) %in% categ_vars,
                   is.factor(variables))) {
      warning('Categorical layers detected in RasterStack do not match with categ_vars.')
    }
    variables <- .remove_cats(variables)
    variables <- st_as_stars(variables) %>% split('band')
  } else {
    if (length(dim(variables)) == 3) {
      variables <- variables %>% split(3)
      # Check categ_vars
      stopifnot(all(categ_vars %in% names(variables)))
    } else {
      # Check categ_vars
      ## step 1
      stopifnot(all(categ_vars %in% names(variables)))
      # step 2
      isfacor <- as.vector(sapply(variables, is.factor))
      if (!identical(names(variables) %in% categ_vars,
                     isfacor)) {
        warning('Categorical layers detected in RasterStack do not match with categ_vars.')
      }
    }
  }

  # Convert to factors
  for (nm in categ_vars) {
    if (!is.factor(variables[[nm]])) {
      variables <- variables %>%
        mutate(!!nm := as.factor(variables[[nm]]))
    }
  }

  # Match observation and spatial variables
  if (st_crs(variables) != st_crs(obs)){
    obs <- st_transform(obs, st_crs(variables))
  }
  if(!is.null(obs_ind_eval)) {
    if (st_crs(variables) != st_crs(obs_ind_eval)){
      obs_ind_eval <- st_transform(obs_ind_eval, st_crs(variables))
    }}

  ############## Reformat observations #############
  # TODO
  # Reformat observations based on obs_mode, mode,
  # contamination
  ##################################################

  # Collect results
  ## Fill out the train and evaluate dataset if necessary
  if (obs_mode != "imperfect_presence") {
    if (obs_mode == "perfect_presence") {
      # Make a template
      rst_template <- st_apply(merge(variables), c("x", "y"),
                               "sum", na.rm = FALSE)
      rst_template[[1]][!is.na(rst_template[[1]])] <- 1

      # Select a contamination percent of background samples
      obs <- .bg_sampling(
        rst_template, obs, seed, floor(contamination * nrow(obs))) %>%
        mutate("observation" = 0) %>%
        rbind(obs)
      rm(rst_template)

      # In case there is very limited background to sample
      contamination <- sum(obs$observation == 0) /
        sum(obs$observation == 1)
    } else if (obs_mode == "presence_absence") {
      obs_full <- obs
      # Select a contamination percent of absences
      obs_frd <- obs %>% filter(.data$observation == 1)
      set.seed(seed)
      obs_frd <- obs_frd %>%
        rbind(obs %>% filter(.data$observation == 0) %>%
                sample_n(min((nrow(obs_frd) * contamination),
                             nrow(.))))
      obs <- obs_frd; rm(obs_frd)
      contamination <- sum(obs$observation == 0) /
        sum(obs$observation == 1)
    }
  }

  ################ Preparation is done ################
  ############### Time to run the model ###############

  # Extract variable values
  obs_vars_mat <- st_extract(x = variables, at = obs) %>%
    st_drop_geometry()

  # Remove NAs from both variable data.frame and spatial points
  ids <- rowSums(is.na(obs_vars_mat)) > 0
  obs_vars_mat <- obs_vars_mat[!ids, ]
  obs <- obs[!ids, ]; rm(ids)

  # Do the same thing to test dataset if there is one
  if (!is.null(obs_ind_eval)) {
    obs_eval_vars_mat <- st_extract(x = variables, at = obs_ind_eval) %>%
      st_drop_geometry()

    # Remove NAs from both variable data.frame and spatial points
    ids <- rowSums(is.na(obs_eval_vars_mat)) > 0
    obs_eval_vars_mat <- obs_eval_vars_mat[!ids, ]
    obs_ind_eval <- obs_ind_eval[!ids, ]; rm(ids)
  } else {
    obs_eval_vars_mat <- NULL}

  # Train the model
  isotree_mod <- isolation.forest(
    obs_vars_mat,
    ntrees = ntrees,
    sample_size = sample_size,
    ndim = ndim,
    seed = seed,
    ...)

  # Do prediction
  ## Raster
  var_pred <- probability(isotree_mod, variables, offset = offset)

  # Generate response curves
  if (response) {
    # marginal variable dependence
    marginal_responses <- marginal_response(
      model = isotree_mod,
      var_occ = obs_vars_mat,
      variables = variables,
      visualize = visualize)

    # independent variable dependence
    independent_responses <- independent_response(
      model = isotree_mod,
      var_occ = obs_vars_mat,
      variables = variables,
      visualize = visualize)

    # Shapley value based variable dependence
    shap_dependences <- shap_dependence(
      model = isotree_mod,
      var_occ = obs_vars_mat,
      variables = variables,
      visualize = visualize,
      seed = seed)
  } else {
    marginal_responses <- NULL
    independent_responses <- NULL
    shap_dependences <- NULL}

  # Spatial dependence
  if (spatial_response) {
    spatial_responses <- spatial_response(
      model = isotree_mod,
      var_occ = obs_vars_mat,
      variables = variables,
      shap_nsim = 0,
      seed = seed,
      visualize = visualize)
  } else {
    spatial_responses <- NULL}

  # Check variable importance
  if (check_variable) {
    vimp <- variable_analysis(
      model = isotree_mod,
      pts_occ = obs,
      pts_occ_test = obs_ind_eval,
      variables = variables,
      visualize = visualize,
      seed = seed)
  } else {vimp <- NULL}

  # Make evaluation using presence-only
  #######################################
  ##         Check the inputs          ##
  #######################################
  # If input obs is presence_absence, then
  # use it as the train evaluation.
  # Check the obs_ind_eval, if it includes
  # absence, then use the absence to do
  # the evaluation.
  # NOTE: if the number of absence is way
  # smaller than number of absence, still
  # use the background samples to eval.

  # Training
  if (obs_mode == "presence_absence") {
    nums <- obs_full %>%
      st_drop_geometry() %>%
      group_by(.data$observation) %>%
      summarise(n = n()) %>%
      arrange(.data$observation)

    # Check the balance of the dataset
    # Just give a warning to let the user deal with
    if ((max(nums$n) - min(nums$n)) / (sum(nums$n) - 2) > 0.5) {
      warning(sprintf("Unbalanced presence-absence obs: %s-%s",
                      nums$n[2], nums$n[1]))
    }

    # Get the prediction
    obs_pred <- st_extract(x = var_pred, at = obs_full) %>%
      st_drop_geometry()
    obs_pred <- data.frame('prediction' = obs_pred$prediction,
                           "observation" = obs_full$observation) %>%
      na.omit()

    ## Calculate
    eval_train <- evaluate_po(
      isotree_mod,
      obs_pred[obs_pred$observation == 1, "prediction"],
      obs_pred[obs_pred$observation == 0, "prediction"],
      na.omit(as.vector(var_pred[[1]])),
      visualize = visualize)

    obs_bg <- NULL
  } else {
    obs_pred <- st_extract(x = var_pred, at = obs) %>%
      st_drop_geometry()

    # Make a template
    rst_template <- st_apply(merge(variables), c("x", "y"),
                             "mean", na.rm = FALSE)
    rst_template[[1]][!is.na(rst_template[[1]])] <- 1

    # Observations
    obs_bg <- .bg_sampling(rst_template, obs, seed, nrow(obs))
    obs_bg_pred <- st_extract(x = var_pred, at = obs_bg) %>%
      st_drop_geometry()
    rm(rst_template)

    ## Calculate
    eval_train <- evaluate_po(
      isotree_mod,
      obs_pred$prediction,
      obs_bg_pred$prediction,
      na.omit(as.vector(var_pred[[1]])),
      visualize = visualize)
  }

  # evaluation
  if (!is.null(obs_ind_eval)) {
    if (0 %in% obs_ind_eval$observation) {
      nums <- obs_ind_eval %>%
        st_drop_geometry() %>%
        group_by(.data$observation) %>%
        summarise(n = n()) %>%
        arrange(.data$observation)

      # Check the balance of the dataset
      # Just give a warning to let the user deal with
      if ((max(nums$n) - min(nums$n)) / (sum(nums$n) - 2) > 0.5) {
        warning(sprintf("Unbalanced presence-absence in obs_ind_eval: %s-%s",
                        nums$n[2], nums$n[1]))
      }

      # Get the prediction
      obs_eval_pred <- st_extract(x = var_pred, at = obs_ind_eval) %>%
        st_drop_geometry()
      obs_eval_pred <- data.frame(
        'prediction' = obs_eval_pred$prediction,
        "observation" = obs_ind_eval$observation) %>%
        na.omit()

      ## Calculate
      eval_test <- evaluate_po(
        isotree_mod,
        obs_eval_pred[obs_eval_pred$observation == 1, "prediction"],
        obs_eval_pred[obs_eval_pred$observation == 0, "prediction"],
        na.omit(as.vector(var_pred[[1]])),
        visualize = visualize)

      obs_eval_bg <- NULL
    } else {
      obs_eval_pred <- st_extract(x = var_pred, at = obs_ind_eval) %>%
        st_drop_geometry()

      # Make a template
      rst_template <- st_apply(merge(variables), c("x", "y"),
                               "mean", na.rm = FALSE)
      rst_template[[1]][!is.na(rst_template[[1]])] <- 1

      # Observations
      obs_eval_bg <- .bg_sampling(rst_template, obs, seed, nrow(obs_ind_eval))
      obs_eval_bg_pred <- st_extract(x = var_pred, at = obs_eval_bg) %>%
        st_drop_geometry()
      rm(rst_template)

      ## Calculate
      eval_test <- evaluate_po(
        isotree_mod,
        obs_eval_pred$prediction,
        obs_eval_bg_pred$prediction,
        na.omit(as.vector(var_pred[[1]])),
        visualize = visualize)
    }
  } else eval_test <- NULL

  # Return
  out <- list(
    model = isotree_mod, # model
    variables = variables, # the formatted stars
    observation = obs, # clean occurrence points
    background_samples = obs_bg, # used background samples paired with training
    independent_test = obs_ind_eval, # clean test points
    # used background samples paired with test
    background_samples_test = obs_eval_bg,
    vars_train = obs_vars_mat, # table of training features
    pred_train = obs_pred, # prediction of training dataset
    eval_train = eval_train, # evaluation on training dataset
    vars_test = obs_eval_vars_mat, # table of test dataset
    pred_test = obs_eval_pred, # prediction of test dataset
    eval_test = eval_test, # evaluation on test dataset
    prediction = var_pred, # prediction map
    offset = offset, # the offset to adjust result
    marginal_responses = marginal_responses, # marginal response plot
    independent_responses = independent_responses, # independent response plot
    shap_dependences = shap_dependences, # Shapley value based response plot
    spatial_responses = spatial_responses, # spatial response maps
    variable_analysis = vimp) # variable evaluation
  class(out) <- append("POIsotree", class(out))

  # Return
  out
}

# isotree_po end
