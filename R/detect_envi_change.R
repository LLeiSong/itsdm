#' @title Detect areas influenced by a changing environment variable.
#' @description Use shapley values to detect the potential areas that will
#' impact the species distribution. It only works on continuous variables.
#' @param model (`isolation_forest` or other model). It could
#' be the item `model` of `POIsotree` made by function \code{\link{isotree_po}}.
#' It also could be other user-fitted models as long as the `pfun` can work on it.
#' @param var_occ (`data.frame`, `tibble`) The `data.frame` style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (`stars`) The `stars` of environmental variables.
#' It should have multiple `attributes` instead of `dims`.
#' If you have `raster` object instead, you
#' could use \code{\link[stars:st_as_stars]{st_as_stars}} to convert it to `stars` or use
#' \code{\link[stars:read_stars]{read_stars}} directly read source data as a `stars`.
#' You also could use item `variables` of `POIsotree` made by function
#' \code{\link{isotree_po}}.
#' @param target_var (`character`) The selected variable to process.
#' @param bins (`integer`) The bin to cut the target variable for the analysis.
#' If it is `NULL`, no cut to apply. The default is `NULL`.
#' @param shap_nsim (`integer`) The number of Monte Carlo repetitions in SHAP
#' method to use for estimating each Shapley value. See details in documentation
#' of function \code{\link[fastshap:explain]{explain}} in package `fastshap`.
#' When the number of variables is large, a smaller shap_nsim could be used.
#' Be cautious that making SHAP-based spatial dependence will be slow
#' because of Monte-Carlo computation for all pixels.
#' But it is worth the time because it is much more
#' informative. See details in documentation of function \code{\link[fastshap:explain]{explain}}.
#' in package `fastshap`. The default is 10. Usually a value 10 - 20 is enough.
#' @param seed (`integer`) The seed for any random progress. The default is `10L`.
#' @param var_future (`numeric` or `stars`) A number to apply to the current
#' variable or a `stars` layer as the future variable. It can be `NULL` if
#' `variables_future` is set.
#' @param variables_future (`stars`) A `stars` raster stack for future variables.
#' It could be `NULL` if `var_future` is set.
#' @param pfun (`function`) The predict function that requires two arguments,
#' `object` and `newdata`.
#' It is only required when `model` is not \code{isolation_forest}.
#' The default is the wrapper function designed for iForest model in `itsdm`.
#' @param method Argument passed on to \code{\link[ggplot2:geom_smooth]{geom_smooth}}
#' to fit the line.
#' Note that the same arguments will be used for all target variables.
#' User could set variable one by one to set the arguments separately.
#' Default value is "gam".
#' @param formula Argument passed on to \code{\link[ggplot2:geom_smooth]{geom_smooth}}
#' to fit the line.
#' Note that the same arguments will be used for all target variables.
#' User could set variable one by one to set the arguments separately.
#' The default is y ~ s(x).
#'
#' @return (`EnviChange`) A list of
#' \itemize{
#' \item{A figure of fitted variable curve}
#' \item{A map of variable contribiution change}
#' \item{Tipping points of variable contribution}
#' \item{A `stars` of variable contribution under current and future condition,
#' and the detected changes}
#' }
#'
#' @details
#' The values show how changes in environmental variable affects the modeling
#' prediction in space. These maps could help to answer questions of where will
#' be affected by a changing variable.
#'
#' @seealso
#' \code{\link{shap_spatial_response}}
#'
#' @importFrom isotree isolation.forest
#' @importFrom dplyr select slice mutate as_tibble pull n %>% group_by
#' @importFrom stars st_as_stars st_dimensions geom_stars st_dimensions<-
#' st_warp
#' @importFrom stats predict setNames uniroot
#' @importFrom tidyselect all_of
#' @importFrom fastshap explain
#' @importFrom ggplot2 geom_point geom_smooth geom_vline geom_text
#' scale_fill_viridis_d
#' @importFrom mgcv gam
#' @importFrom raster reclassify
#' @export
#' @examples
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
#'   slice('band', c(1, 12))
#' #'
#' # With imperfect_presence mode,
#' mod <- isotree_po(
#'   obs_mode = "imperfect_presence",
#'   obs = obs_train_eval$obs,
#'   obs_ind_eval = obs_train_eval$eval,
#'   variables = env_vars, ntrees = 5,
#'   sample_size = 0.8, ndim = 1L,
#'   nthreads = 1,
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
#' \dontrun{
#' # Use a future layer
#' ## Read the future Worldclim variables
#' future_vars <- system.file(
#'   'extdata/future_bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   split() %>% select(bioc1, bioc12)
#' # Rename the bands
#' names(future_vars) <- paste0("bio", c(1, 12))
#'
#' ## Just use the target future variable
#' climate_changes <- detect_envi_change(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables,
#'   shap_nsim = 1,
#'   target_var = "bio1",
#'   var_future = future_vars %>% select("bio1"))
#'
#' ## Use the whole future variable tack
#' bio12_changes <- detect_envi_change(
#'   model = mod$model,
#'   var_occ = mod$vars_train,
#'   variables = mod$variables,
#'   shap_nsim = 1,
#'   target_var = "bio12",
#'   variables_future = future_vars)
#'
#' print(bio12_changes)
#'
#' ##### Use Random Forest model as an external model ########
#' library(randomForest)
#'
#' # Prepare data
#' data("occ_virtual_species")
#' obs_df <- occ_virtual_species %>%
#'   filter(usage == "train")
#'
#' env_vars <- system.file(
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12)) %>%
#'   split()
#'
#' model_data <- stars::st_extract(
#'   env_vars, at = as.matrix(obs_df %>% select(x, y))) %>%
#'   as.data.frame()
#' names(model_data) <- names(env_vars)
#' model_data <- model_data %>%
#'   mutate(occ = obs_df[['observation']])
#' model_data$occ <- as.factor(model_data$occ)
#'
#' mod_rf <- randomForest(
#'   occ ~ .,
#'   data = model_data,
#'   ntree = 200)
#'
#' pfun <- function(X.model, newdata) {
#'   # for data.frame
#'   predict(X.model, newdata, type = "prob")[, "1"]
#' }
#'
#' # Use a fixed value
#' bio5_changes <- detect_envi_change(
#'   model = mod_rf,
#'   var_occ = model_data %>% select(-occ),
#'   variables = env_vars,
#'   target_var = "bio5",
#'   bins = 20,
#'   var_future = 5,
#'   pfun = pfun)
#'
#' plot(bio5_changes)
#'}
#'
# Now just for continuous variables
detect_envi_change <- function(model,
                               var_occ,
                               variables,
                               target_var,
                               bins = NULL,
                               shap_nsim = 10,
                               seed = 10,
                               # a number or a stars
                               var_future = NULL,
                               variables_future = NULL,
                               pfun = .pfun_shap,
                               method = 'gam',
                               formula = y ~ s(x)){
  # Check inputs
  ## required
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_int(bins, null.ok = TRUE)
  checkmate::assert_class(variables, 'stars')
  stopifnot(length(dim(variables)) <= 2)
  checkmate::assert_character(target_var)
  checkmate::assert_int(shap_nsim)
  checkmate::assert_int(seed)
  bands <- names(variables)
  stopifnot(all(bands %in% colnames(var_occ)))

  ## optional
  checkmate::assert_class(
    variables_future, 'stars',
    null.ok = TRUE)
  if (!is.null(variables_future)) {
    message("Set variables_future, will use all future variables.")
    stopifnot(length(dim(variables_future)) <= 2)
    stopifnot(all(names(variables_future) %in% bands))
  } else {
    message("Just set the single future variable.")
    checkmate::assert_multi_class(
      var_future, classes = c("numeric", "stars"))

    if (inherits(var_future, "numeric")) {
      message(sprintf("Change current %s with %s.", target_var, var_future))
    } else {
      message(sprintf("Change current %s to a new one.", target_var))
    }
  }

  ## Split numeric and categorical
  bands <- names(variables)
  isfacor <- as.vector(sapply(variables, is.factor))
  bands_cont <- bands[!isfacor]
  bands_cat <- bands[isfacor]

  # response curve analysis
  # detect the cross points between response curve and y = 0.
  shap_dependences <- shap_dependence(
    model = model,
    var_occ = var_occ,
    variables = variables,
    visualize = FALSE,
    seed = seed,
    pfun = pfun)

  # Make the plot
  p_curve <- plot(shap_dependences,
                  sample_prop = 1.0,
                  target_var = target_var,
                  smooth_line = FALSE) +
    geom_smooth(color = 'red',
                alpha = 0,
                # Use GAM
                method = method, formula = formula,
                linewidth = 1)

  # Get the intersection with y = 0
  ## In case curve failed to make
  tryCatch(
    expr = {
      xy <- ggplot_build(p_curve)$data[[1]][, c("x", "y")]
      fit_line <- gam(formula = y ~ s(x), data = xy)
      fit_fun <- function(x) predict(fit_line, data.frame(x = x))
    },
    error = function(e) {
      stop(paste0("Failed to generate response curve. ",
                  "Please check the inputs, such as ",
                  "variables or pfun."))
    })

  # In case no curve generated above
  if(inherits(try(
    roots <- .find_intersects(fit_fun, xy$x)
  ), "try-error")){

  }

  # Extend the plot
  if (length(roots) > 0) {
    p_curve <- p_curve +
      geom_point(data = data.frame(x = roots, y = rep(0, length(roots))),
                 aes(x = .data$x, y = .data$y), size = 3, color = "blue") +
      geom_vline(xintercept = min(roots), linetype = "dotted") +
      geom_vline(xintercept = max(roots), linetype = "dotted") +
      geom_text(data = data.frame(x = roots, y = rep(min(xy$y), length(roots))),
                aes(x = .data$x, y = .data$y, label = round(.data$x, 2))) +
      ylab("Shapley value")
  }

  # Spatial map analysis
  # Bin cut the maps, and then calculate the spatial map
  # and then make the change map: 0 -> 1, 1 -> 0
  if (is.null(variables_future)) {
    ## Extract the layers
    var_old <- variables[target_var]
    if (inherits(var_future, "numeric")){
      var_future <- variables[target_var] + var_future}

    if (!is.null(bins)) {
      dims <- st_dimensions(var_old)
      # Cut the maps
      breaks <- seq(min(var_old[[1]], na.rm = TRUE),
                    max(var_old[[1]], na.rm = TRUE),
                    floor(max(var_old[[1]], na.rm = TRUE) -
                            min(var_old[[1]], na.rm = TRUE)) / bins)
      var_old <- reclassify(
        as(var_old, "Raster"),
        rcl = cbind(breaks[-length(breaks)], breaks[-1],
                    (breaks[-length(breaks)] - breaks[-1]) / 2 + breaks[-1])) %>%
        st_as_stars()
      names(var_old) <- target_var; st_dimensions(var_old) <- dims

      breaks <- seq(min(var_future[[1]], na.rm = TRUE),
                    max(var_future[[1]], na.rm = TRUE),
                    floor(max(var_future[[1]], na.rm = TRUE) -
                            min(var_future[[1]], na.rm = TRUE)) / bins)
      var_future <- reclassify(
        as(var_future, "Raster"),
        rcl = cbind(breaks[-length(breaks)], breaks[-1],
                    (breaks[-length(breaks)] - breaks[-1]) / 2 + breaks[-1])) %>%
        st_as_stars()
      names(var_future) <- target_var; st_dimensions(var_future) <- dims
    } else {
      var_future <- st_warp(var_future, var_old, no_data_value = 0,
                            use_gdal = TRUE, method = "bilinear")
      names(var_future) <- target_var
    }

    # change the new variable
    variables <- variables %>% select(-all_of(target_var))
    variables <- c(variables, var_old) %>%
      select(all_of(bands))

    variables_future <- variables %>% select(-all_of(target_var))
    variables_future <- c(variables_future, var_future) %>%
      select(all_of(bands))
  }

  shap_old <- shap_spatial_response(
    model,
    var_occ,
    variables,
    target_vars = target_var,
    shap_nsim,
    seed,
    pfun = pfun)

  shap_future <- shap_spatial_response(
    model,
    var_occ,
    variables_future,
    target_vars = target_var,
    shap_nsim,
    seed,
    pfun = pfun)

  # Collect result
  if (!identical(st_dimensions(shap_old[[1]]),
                st_dimensions(shap_future[[1]]))) {
    shap_future[[1]] <- st_warp(
      shap_future[[1]], shap_old[[1]], no_data_value = 0,
      use_gdal = TRUE, method = "near")
  }
  shap_change <- c(shap_old[[1]], shap_future[[1]])
  names(shap_change) <- c("current", "future")
  shap_change <- shap_change %>% mutate(
    change = ifelse(.data$current > 0 & .data$future > 0, 1,
                    ifelse(.data$current > 0 & .data$future <= 0, 2,
                    ifelse(.data$current <= 0 & .data$future > 0, 3, 4)))) %>%
    mutate(change = factor(
      .data$change, levels = 1:4,
      labels = c("Positive-Positive", "Positive-Negative",
                 "Negative-Positive", "Negative-Negative")))

  p_map <- ggplot() +
    geom_stars(data = shap_change %>% select(.data$change),
               na.action = na.omit) + coord_equal() +
    scale_fill_viridis_d("Contribution change") + theme_void()

  # Out
  out <- list("p_curve" = p_curve,
              "p_map" = p_map,
              "Tipping_points" = roots,
              "variable_contribution_change" = shap_change)
  class(out) <- append("EnviChange", class(out))

  # Return
  return(out)
}
