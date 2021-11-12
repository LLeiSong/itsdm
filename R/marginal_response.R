#' @title A function to calculate marginal responses of each variables.
#' @description This function allows you to calculate the marginal responses of
#' each variables within the model.
#' @param model Any predictive model. In this package, it could be isolation_forest.
#' In this package, it is isolation_forest.
#' @param var_occ (data.frame, tibble) the data.frame style table that
#' include values of environmental variables at occurrence locations.
#' @param variables (stars) The stars of environmental variables. It should have
#' multiple attributes instead of dims. If you have `raster` object instead, you
#' could use `st_as_stars` to convert it to `stars` or use `read_stars` directly
#' read source data as a `stars`.
#' @param si (integer) the number of samples to generate response curves.
#' If it is too small, the response curves might be biased.
#' The default value is 1000.
#' @param visualize (logical) if TRUE, plot the response curves.
#' The default is FALSE.
#' @details The curves show how each environmental variable affects the modeling
#' prediction. The curves show how the predicted result changes as each environmental
#' variable is varied while keeping all other environmental variables at average
#' sample value. The curves might be hard to interpret if there are strongly correlated
#' variables. The users could use `dim_reduce` function to remove the strong correlation
#' from original environmental variable stack.
#' @return (MarginalResponse) a nested list of data.frame of response and variable values.
#' The response values correspond to suitability of this single variable.
#' @importFrom dplyr select slice as_tibble pull
#' @importFrom stars st_as_stars
#' @export
#' @examples
#' marginal_response(model = mod$model, variables = env_vars)
marginal_response <- function(model,
                              var_occ,
                              variables,
                              si = 1000,
                              visualize = FALSE) {

  # Check inputs
  checkmate::assert_data_frame(var_occ)
  checkmate::assert_int(si)
  checkmate::assert_class(variables, 'stars')
  stopifnot(length(dim(variables)) == 2)
  checkmate::assert_logical(visualize)
  bands <- names(variables)
  stopifnot(all(bands %in% colnames(var_occ)))

  # Reformat data
  ## Split numeric and categorical
  isfacor <- as.vector(sapply(variables, is.factor))
  bands_cont <- bands[!isfacor]
  bands_cat <- bands[isfacor]

  # Numeric variables
  ## Not limited to data volume, could generate as many as possible
  ## pseudo observations in [min, max], so the function could make
  ## smoother curves.
  if (length(bands_cont) > 0) {
    mins <- sapply(bands_cont, function(nm) {
      min(variables %>% select(nm) %>% pull, na.rm = T)})
    maxs <- sapply(bands_cont, function(nm) {
      max(variables %>% select(nm) %>% pull, na.rm = T)})
    vals_cont <- do.call(cbind, lapply(1:length(mins), function(n) {
      seq(from = mins[n], to = maxs[n],
          length.out = si)})) %>%
      data.frame() %>%
      setNames(bands_cont)

    # Get means and expand
    means <- lapply(bands, function(nm) {
      if (nm %in% bands_cont) {
        mean(var_occ %>% pull(nm), na.rm = T)
      } else {
        .mode(var_occ %>% pull(nm))
      }}) %>%
      setNames(bands) %>%
      as_tibble() %>%
      slice(rep(1:n(), each = nrow(vals_cont)))
    responses_cont <- lapply(names(vals_cont), function(nm) {
      vals_tmp <- means
      vals_tmp <- vals_tmp %>%
        mutate('{nm}' := vals_cont %>% pull(nm))
      pred_tmp <- predict(model, vals_tmp)
      pred_tmp <- 1 - pred_tmp
      pred_tmp <- .norm(pred_tmp)
      data.frame(y = pred_tmp,
                 x = vals_cont %>% pull(nm)) %>%
        setNames(c("y", "x"))
    })
    names(responses_cont) <- bands_cont
    rm(mins, maxs, vals_cont, means)
  } else responses_cont <- NULL

  # Categorical variables
  ## The number of pseudo observations is limited to factor levels
  ## So categorical variables should generate one by one
  if (length(bands_cat) > 0) {
    responses_cat <- lapply(bands_cat, function(nm) {
      vals_this <- var_occ %>% pull(nm) %>%
        levels() %>% as.factor()

      means_this <- lapply(bands, function(nm) {
        if (nm %in% bands_cont) {
          mean(var_occ %>% pull(nm), na.rm = T)
        } else {
          .mode(var_occ %>% pull(nm))
        }}) %>%
        setNames(bands) %>%
        as_tibble() %>%
        slice(rep(1:n(), each = length(vals_this)))

      means_this <- means_this %>%
        mutate('{nm}' := vals_this)
      pred_tmp <- predict(model, means_this)
      pred_tmp <- 1 - pred_tmp
      pred_tmp <- .norm(pred_tmp)
      data.frame(y = pred_tmp,
                 x = vals_this) %>%
        setNames(c("y", "x"))
    }) %>% setNames(bands_cat)
  } else responses_cat <- NULL

  # Put together
  responses <- list(responses_cont = responses_cont,
                    responses_cat = responses_cat)
  class(responses) <- append("MarginalResponse", class(responses))

  # Visualize
  if (visualize) {
    plot(responses)
  }

  # Return
  return(responses)
}

# marginal_response end
