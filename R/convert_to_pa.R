#' @title Convert predicted suitability to presence-absence map.
#' @description Use threshold-based, logistic or linear conversion method to
#' convert predicted suitability map to presence-absence map.
#' @param suitability (`stars` or `RasterLayer`) The suitability raster.
#' @param method (`character`) The conversion method, must be one of
#' 'threshold', 'logistic', and 'linear'. The default is 'logistic'.
#' @param beta (`numeric`) Works for 'threshold' or 'logistic' method.
#' If `method` is threshold, then `beta` is the threshold value to cutoff.
#' If `method` is logistic, it is the sigmoid midpoint. The default is `0.5`.
#' @param alpha (`numeric`) Works for logistic method.
#' It is the logistic growth rate or steepness of the curve.
#' The default is `-.05`.
#' @param a (`numeric`) Works for linear method. It is the slope of the line.
#' The default is `1`.
#' @param b (`numeric`) Works for linear method.
#' It is the intercept of the line. The default is `0`.
#' @param species_prevalence (`numeric` or `NA`) Works for all three methods.
#' It is the species prevalence to classify suitability map.
#' It could be `NA`, when the will be calculated automatically
#' based on other arguments. The default is `NA`.
#' @param threshold (`numeric`) The threshold used to convert probability of
#' occurrence to presence-absence map. It ranges in `[0, 1]`. The default is 0.5.
#' @param seed (`integer`) The seed for random progress. The default is `10L`
#' @param visualize (`logical`) If `TRUE`, plot map of suitability,
#' probability of occurrence, and presence-absence together.
#' The default is `TRUE`.
#' @return (`PAConversion`) A list of
#' \itemize{
#' \item{suitability (`stars`) The input suitability map}
#' \item{probability_of_occurrence (`stars`) The map of occurrence probability}
#' \item{pa_conversion (`list`) A list of conversion arguments}
#' \item{pa_map (`stars`) The presence-absence map}}
#'
#' @seealso
#' \code{\link{plot.PAConversion}}
#'
#' @details
#' Multiple methods and arguments could be used as a combination to do
#' the conversion.
#' @references
#' \href{https://github.com/Farewe/virtualspecies/blob/master/R/convertToPA.R}{c
#' onvertToPA in package `virtualspecies`}
#'
#' @importFrom stars st_as_stars
#' @importFrom dplyr case_when
#' @importFrom stats na.omit
#' @importFrom methods is
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
#'   'extdata/bioclim_tanzania_10min.tif',
#'   package = 'itsdm') %>% read_stars() %>%
#'   slice('band', c(1, 5, 12, 16))
#'
#' mod <- isotree_po(
#'   occ = occ, occ_test = occ_test,
#'   variables = env_vars, ntrees = 50,
#'   sample_size = 0.8, ndim = 3L,
#'   seed = 123L, response = FALSE,
#'   spatial_response = FALSE,
#'   check_variable = FALSE)
#'
#' # Threshold conversion
#' pa_thred <- convert_to_pa(mod$prediction,
#'   method = 'threshold', beta = 0.5)
#' pa_thred
#' plot(pa_thred)
#'
#' # Logistic conversion
#' pa_log <- convert_to_pa(mod$prediction, method = 'logistic',
#'   beta = 0.5, alpha = -.05)
#'
#' # Linear conversion
#' pa_lin <- convert_to_pa(mod$prediction, method = 'linear',
#'   a = 1, b = 0)
#'
convert_to_pa <- function(suitability, # prediction from isotree_sdm
                          method = "logistic",
                          beta = 0.5, # could be NA, for threshold or logistic
                          alpha = -.05, # only for logistic
                          a = 1, # for linear
                          b = 0, # for linear
                          species_prevalence = NA,  # could be NA, for all
                          threshold = 0.5, # threshold to convert to PA
                          seed = 10L,
                          visualize = TRUE) {
  # Check inputs - level 1
  checkmate::assert_multi_class(
    suitability, c('RasterLayer', 'stars'))
  # Convert suitability if it is a raster
  if (is(suitability, 'RasterLayer')){
    suitability <- st_as_stars(suitability)}

  checkmate::assert_choice(
    method, c('threshold', 'logistic', 'linear'))
  suit_min <- .min_value(suitability)
  suit_max <- .max_value(suitability)
  checkmate::assert_number(beta, lower = suit_min,
                           upper = suit_max, na.ok = T)
  checkmate::assert_number(alpha, na.ok = T)
  checkmate::assert_number(a, na.ok = T)
  checkmate::assert_number(b, na.ok = T)
  checkmate::assert_number(species_prevalence, lower = 0,
                           upper = 1, na.ok = T)
  checkmate::assert_number(threshold,  lower = 0, upper = 1)
  checkmate::assert_int(seed)
  checkmate::assert_logical(visualize)

  # Check inputs - level 2
  # threshold - beta or species_prevalence
  if (method == 'threshold' & (is.na(beta) & is.na(species_prevalence))) {
    stop('Must set beta or species prevalence for threshold conversion.')
  }
  # logistic - two of alpha, beta, and species_prevalence
  # or just species_prevalence
  if (method == 'logistic') {
    if (sum(is.na(c(beta, alpha, species_prevalence))) == 3){
      stop(paste0('No beta, alpha, or species_prevalence',
                  ' is set for logistic conversion.'))
    } else if (sum(is.na(c(beta, alpha, species_prevalence))) == 2){
      if (is.na(species_prevalence)) {
        stop('species_prevalence must be set if beta and alpha is missing.')}
      }
    }
  # linear - a and b or species_prevalence, could be all NAs, so skip check.

  # Convert arguments
  ## threshold conversion: if both set, just use prevalence
  if (method == 'threshold') {
    if (!is.na(beta) & !is.na(species_prevalence)) {
      warning('Both beta and species_prevalence are set. Ignore beta.')
      beta <- NA}
    }

  ## logistic conversion
  ## if all set, set beta to NA
  ## just prevalence, choose alpha randomly and decide beta programmatically
  ## keep two of three non-NA
  if (method == 'logistic') {
    if (length(na.omit(c(beta, alpha, species_prevalence))) == 1){
      if (!is.na(species_prevalence)) {
        warning(
          paste0('Neither alpha or beta set. Randomly select alpha and',
                 ' adjust beta programmatically to the desired prevalence.'))
        set.seed(seed)
        alpha <- -sample(c(seq((suit_max - suit_min)/1000,
                               (suit_max - suit_min)/100, length = 10),
                           seq((suit_max - suit_min)/100,
                               (suit_max - suit_min)/10, length = 100),
                           seq((suit_max - suit_min)/10,
                               (suit_max - suit_min)*10, length = 10)),
                         size = 1)
      }
    } else if (length(na.omit(c(beta, alpha, species_prevalence))) == 3) {
      beta <- NA}

    # Provide necessary message
    if (!is.na(species_prevalence)) {
      if (!is.na(beta)) {
        message(paste0('species_prevalence and beta are set,',
                       ' choose alpha automatically.'))
      } else {
        message(paste0('species_prevalence and alpha are set,',
                       ' choose beta automatically.'))
      }
    } else {
      message(paste0('alpha and beta are set,',
                     ' choose species_prevalence automatically.'))
    }
  }

  ## linear conversion
  ## Nothing set, use a = 1 and b = 0 instead
  ## if species_prevalence set, just use it
  ## if not, must set a and b to get a linear conversion.
  if (method == 'linear') {
    if (is.na(a) & is.na(b) & is.na(species_prevalence)) {
      message(paste0('No arguments set, set slope = 1 and intercept = 0',
                     ' for linear conversion.'))
      a <- 1; b <- 0
    }
    if (is.na(species_prevalence)) {
      if (is.na(a) | is.na(b)) {
        stop('Must set a and b if no species_prevalence provided.')
      } else {
        message('a and b are set, choose species_prevalence automatically.')
      }
    } else {
      message(paste0('Ignore a and b if set, search for a linear conversion',
                     ' that fits species_prevalence.'))
      a <- NA; b <- NA
    }
  }

  # Start conversion
  ## threshold conversion
  if (method == 'threshold') {
    if(!is.na(species_prevalence)) {
      beta <- .quantile_stars(suitability, 1 - species_prevalence)
      names(beta) <- NULL}

    pa_map <- suitability %>%
      mutate(prediction = case_when(prediction < beta ~ 0,
                                    prediction >= beta ~ 1))
    prob_of_occurrence <- pa_map
    # Convert to binary
    pa_map <- pa_map == 1
  }

  ## logistic conversion
  if (method == 'logistic') {
    if (!is.na(species_prevalence)) {
      if (!is.na(beta)) {
        alpha_test <- do.call(rbind, lapply(
          c((.max_value(suitability) - .min_value(suitability))/1000,
            (.max_value(suitability) - .min_value(suitability)) * 10),
          function(alpha) {
            if(alpha > 0) alpha <- -alpha
            prob_of_occurrence <- .logistic(
              suitability, beta = beta, alpha = alpha)
            pa_map <- .binary_convert(prob_of_occurrence, threshold = threshold)
            c(alpha, .mean_value(pa_map))}))
        epsilon <- species_prevalence - alpha_test[, 2]
        if(all(epsilon > 0)){
          warning(
            paste0(
              'The desired species prevalence cannot be obtained, ',
              'because of the chosen beta and available environmental',
              ' conditions.\n',
              sprintf("The closest possible estimate of prevalence was %s.",
                      round(alpha_test[2, 2], 2)),
              "\nPerhaps you can try a lower beta value."))
          alpha <- alpha_test[2, 1]
        } else if (all(epsilon < 0)){
          warning(
            paste0(
              'The desired species prevalence cannot be obtained, ',
              'because of the chosen beta and available environmental',
              ' conditions.\n',
              sprintf("The closest possible estimate of prevalence was %s.",
                      round(alpha_test[1, 2], 2)),
              "\nPerhaps you can try a higher beta value."))
          alpha <- alpha_test[1, 1]
        } else {
          while (all(abs(epsilon) > 0.01)){
            alpha <- (alpha_test[which(epsilon == max(epsilon[epsilon < 0])), 1] +
                        alpha_test[which(epsilon == min(epsilon[epsilon > 0])), 1]) / 2
            prob_of_occurrence <- .logistic(suitability, beta = beta, alpha = alpha)
            pa_map <- .binary_convert(prob_of_occurrence, threshold = threshold)
            alpha_test <- rbind(alpha_test, c(alpha, .mean_value(pa_map)))

            epsilon <- species_prevalence - alpha_test[, 2]
          }
        }
      } else {
        # We define the upper and lower boundaries for beta.
        # We choose to be able to define beta values beyond the boundaries of
        # our probability of occurrence raster, so we can have a larger range
        # of prevalence
        aa <- .min_value(suitability) -
          diff(c(.min_value(suitability), .max_value(suitability))) / 2
        bb <- .max_value(suitability) +
          diff(c(.min_value(suitability), .max_value(suitability))) / 2
        beta_test <- do.call(rbind, lapply(c(aa, bb), function(beta) {
          prob_of_occurrence <- .logistic(suitability,
                                          beta = beta, alpha = alpha)
          pa_map <- .binary_convert(prob_of_occurrence, threshold = threshold)
          c(beta, .mean_value(pa_map))}))

        epsilon <- data.frame(epsi = species_prevalence - beta_test[, 2],
                              prevalence = beta_test[, 2])
        if(all(epsilon$epsi > 0)) {
          warning(
            paste0(
              'The desired species prevalence cannot be obtained, ',
              'because of the chosen alpha and available environmental',
              ' conditions.\n',
              sprintf("The closest possible estimate of prevalence was %s.",
                      round(beta_test[1, 2], 2)),
              "\nPerhaps you can try an alpha value closer to 0."))
          beta <- beta_test[1, 1]
        } else if (all(epsilon$epsi < 0)) {
          warning(
            paste0(
              'The desired species prevalence cannot be obtained, ',
              'because of the chosen alpha and available environmental',
              ' conditions.\n',
              sprintf("The closest possible estimate of prevalence was %s.",
                      round(beta_test[2, 2], 2)),
              "\nPerhaps you can try an alpha value closer to 0."))
          beta <- beta_test[2, 1]
        } else {
          while (all(apply(
            epsilon, 1, function(x) {
              ifelse(abs(x[1]) > 0.001, TRUE,
                     ifelse(x[2] == 0, TRUE, FALSE))}
            ))){
            beta <- (beta_test[which(epsilon$epsi == max(epsilon$epsi[epsilon$epsi < 0])), 1][1] +
                       beta_test[which(epsilon$epsi == min(epsilon$epsi[epsilon$epsi > 0])), 1][1]) / 2
            prob_of_occurrence <- .logistic(suitability, beta = beta, alpha = alpha)
            pa_map <- .binary_convert(prob_of_occurrence, threshold = threshold)

            beta_test <- rbind(beta_test, c(beta, .mean_value(pa_map)))
            epsilon <- data.frame(epsi = species_prevalence - beta_test[, 2],
                                  prevalence = beta_test[, 2])
          }
        }
      }
    }

    # Convert
    prob_of_occurrence <- .logistic(suitability, beta = beta, alpha = alpha)
    pa_map <- .binary_convert(prob_of_occurrence, threshold = threshold)

    # Adjust result for very low prevalence
    if (.max_value(suitability) == 0){
      while (.max_value(suitability) == 0) {
        prob_of_occurrence <- .logistic(suitability, beta = beta, alpha = alpha)
        pa_map <- .binary_convert(prob_of_occurrence, threshold = threshold)
      }
    }
  }

  # linear conversion
  if (method == 'linear') {
    if (!is.na(species_prevalence)) {
      tmp <- .find_linear_conversion(
        suitability, species_prevalence,
        threshold)
      a <- tmp$a
      b <- tmp$b
      prob_of_occurrence <- tmp$prob_of_occurrence
      pa_map <- tmp$distribution
    } else {
      prob_of_occurrence <- .linear_convert(suitability, c(a, b))
      pa_map <- .binary_convert(prob_of_occurrence, threshold = threshold)
    }
  }

  # summarize results
  species_prevalence <- .mean_value(pa_map)

  if(method == "threshold") {
    pa_conversion <-  list(cutoff = beta,
                         conversion_method = method,
                         species_prevalence = species_prevalence)
  } else if(method == "logistic") {
    pa_conversion <- list(conversion_method = method,
                         alpha = alpha,
                         beta = beta,
                         species_prevalence = species_prevalence)
  } else if(method == "linear") {
    pa_conversion <- list(conversion_method = method,
                         a = a,
                         b = b,
                         species_prevalence = species_prevalence)
  }
  out <- list(suitability = suitability,
              probability_of_occurrence = prob_of_occurrence,
              pa_conversion = pa_conversion,
              pa_map = pa_map)
  class(out) <- append("PAConversion", class(out))

  if (visualize) plot(out)
  # Return
  return(out)
}

