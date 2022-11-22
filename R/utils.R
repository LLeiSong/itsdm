# A few internal functions used in this package
# Calculate the mean of absolute values in a vector
.abs_mean <- function(v) mean(abs(v))

# Get spatial resolution of stars
#' @importFrom stars st_dimensions
.res_stars <- function(s){
  x <- abs(st_dimensions(s)$x$delta)
  y <- abs(st_dimensions(s)$y$delta)
  list(x = x,
       y = y)
}

# Get ROC_ratio
.roc_ratio <- function(occ, full) {
  roc_r <- lapply(seq(0, 1, 0.01), function(t){
    ratio_test <- sum(occ >= t) / length(occ)
    ratio_all <- sum(full >= t) / length(full)
    data.frame(presence = ratio_test,
               cell = ratio_all)
  })
  do.call(rbind, roc_r)
}

# Approximately calculate AUC_ratio
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
.auc_ratio <- function(occ, full) {
  roc_r <- .roc_ratio(occ, full)
  roc_r <- roc_r %>% arrange(.data$cell)
  sum(sapply(2:nrow(roc_r), function(n){
    (roc_r$cell[n] - roc_r$cell[n - 1]) * roc_r$presence[n - 1] +
      (roc_r$cell[n] - roc_r$cell[n - 1]) *
      (roc_r$presence[n] - roc_r$presence[n - 1]) / 2
  }))
}

# A norm function
.norm <- function(val) {
  min <- min(val)
  max <- max(val)
  (val - min) / (max - min)
}
# Logistic transfer function
.logistic <- function(orig_values,
                      beta = 0.5,
                      alpha = 0.05) {
  1 / (1 + exp((orig_values - beta) / alpha))
}

# Functions to get min max, mean and std values of a single band stars
.min_value <- function(x) min(x[[1]], na.rm = T)
.max_value <- function(x) max(x[[1]], na.rm = T)
.mean_value <- function(x) mean(x[[1]], na.rm = T)

# Function to get sd
#' @importFrom stats sd
.std_value <- function(x) sd(x[[1]], na.rm = T)

# Function to get erf of a single band stars
.erf_stars <- function(x) {
  vals <- x[[1]]
  vals <- 2 * stats::pnorm(vals * sqrt(2)) - 1
  x[[1]] <- vals
  x
}

# Linear stretch of single-band stars
#' @importFrom stats quantile
.stars_stretch <- function(x, # stars
                           new_values = NULL,
                           minv = 0,
                           maxv = 1,
                           minq = 0, # 0.5
                           maxq = 1) {
  # Check inputs
  checkmate::assert_class(x, 'stars')
  checkmate::assert_multi_class(
    new_values, c('numeric', 'RasterLayer', 'stars'),
    null.ok = T)
  checkmate::assert_number(minv)
  checkmate::assert_number(maxv)
  stopifnot(maxv > minv)

  checkmate::assert_number(minq)
  checkmate::assert_number(maxq)
  stopifnot(maxq > minq)

  # Convert values
  minq <- max(0, minq)
  maxq <- min(1, maxq)
  stopifnot(minq < maxq)

  if (minq == 0 & maxq == 1) {
    q <- cbind(.min_value(x), .max_value(x))
  } else {
    q <- quantile(x[[1]], c(minq, maxq), na.rm = TRUE)
  }

  # Stretch values
  if (is.null(new_values)) {
    x_strech <- x
  } else {
    x_strech <- new_values
  }

  mult <- maxv / (q[2] - q[1])
  x_strech <- mult * (x_strech - q[1])
  x_strech[x_strech < minv] <- minv
  x_strech[x_strech > maxv] <- maxv

  x_strech
}

# Linear stretch of a vector, paralleling to .stars_stretch
#' @importFrom stats quantile
.stretch <- function(x,
                     new_values = NULL,
                     minv = 0,
                     maxv = 1,
                     minq = 0.5,
                     maxq = 1) {
  # Convert values
  minq <- max(0, minq)
  maxq <- min(1, maxq)
  stopifnot(minq < maxq)

  if (minq == 0 & maxq == 1) {
    q <- cbind(min(x), max(x))
  } else {
    q <- quantile(x, c(minq, maxq), na.rm = TRUE)
  }

  # Stretch values
  if (is.null(new_values)) {
    x_strech <- x
  } else {
    x_strech <- new_values
  }

  mult <- maxv / (q[2] - q[1])
  x_strech <- mult * (x_strech - q[1])
  x_strech[x_strech < minv] <- minv
  x_strech[x_strech > maxv] <- maxv

  x_strech
}

# Predict_wrapper for SHAP
#' @importFrom stats predict
.pfun_shap <- function(X.model, newdata) {
  pred <- 1 - predict(X.model, newdata)
  #.stretch(pred)
}

# Functions related to convert_to_pa
# Calculate quantile of stars
#' @importFrom stats quantile
#'
.quantile_stars <- function(x, # stars with one band
                            ...,
                            na.rm = TRUE) {
  v <- try(as.vector(x[[1]]))
  return(quantile(v, ..., na.rm = na.rm))
}

# Functions to find linear conversion
# Reference: https://github.com/Farewe/virtualspecies/blob/master/R/convertToPA.R
# Get line coefficients from two points
.abcoefs <- function(x1, y1, x2, y2) {
  list(b = y1 - x1 * (y1 - y2) / (x1 - x2),
       a = (y1 - y2) / (x1 - x2))
}

# Function for a line with intercept (b) and slope (a)
.lab <- function(x, b, a) a * x + b

# Function to convert
.binary_convert <- function(prob_of_occurrence,
                            threshold = 0.5,
                            ...) {
  # Make binary
  prob_of_occurrence >= threshold
}

.find_linear_conversion <- function(suitability,
                                    target_prevalence,
                                    threshold = 0.5) {
  suit_max <- .max_value(suitability)
  suit_mean <- .mean_value(suitability)
  suit_min <- .min_value(suitability)

  xs <- c(suit_min, suit_max)
  ys <- c(0, 1)

  # Only include (0, 0) case if suitability >= 0
  if (suit_min >= 0) {
    xs <- c(0, xs)
    ys <- c(0, ys)
  }

  AB <- .abcoefs(suit_mean,
                target_prevalence,
                xs,
                ys)

  ymn <- .lab(suit_min, AB$a, AB$b)
  ymx <- .lab(suit_max, AB$a, AB$b)

  # Round to avoid very small floating point calculation errors
  ymn <- round(ymn, 6)
  ymx <- round(ymx, 6)

  I <- min(which(ymn >= 0 & ymx <= 1)) # Find first one that works

  # Calculate the resulting prevalence:
  new_suit <- AB$a[I] * suitability + AB$b[I]
  distr <- .binary_convert(new_suit, threshold = threshold)
  prev <- .mean_value(distr)

  return(list(a = AB$a[I],
              b = AB$b[I],
              prevalence = prev,
              prob_of_occurrence = new_suit,
              distribution = distr))
}

# Linear transfer
.linear_convert <- function(x, coefs) {
  x <- x * coefs[1] + coefs[2]
  if(.min_value(x) < 0 | .max_value(x) > 1) {
    if(.min_value(x) < 0 & .max_value(x) > 1) {
      message(paste0('The linear transformation resulted in probability values ',
      'below 0 and above 1, so these were respectively truncated to 0 and 1.\n'))
    } else if(.min_value(x) < 0) {
      message(paste0('The linear transformation resulted in probability values',
                     ' below 0 so these were truncated to 0\n'))
    } else {
      message(paste0('The linear transformation resulted in probability values',
                     ' above 1 so these were truncated to 1\n'))
    }
  }
  x[x < 0] <- 0
  x[x > 1] <- 1
  return(x)
}

# kruskal.test between RasterStack and RasterLayer
## Return the p values
# .kruskal.test.raster <- function(rst_stack,
#                                  cat_rst){
#   # Convert data
#   cats <- getValues(cat_rst)
#   cors <- sapply(1:nlayers(rst_stack), function(n) {
#     vals <- getValues(subset(rst_stack, n))
#     # Calculate
#     kruskal.test(x = vals, g = cats, na.action = 'na.omit')[['p.value']]
#   })
#   data.frame(cors) %>% setNames(names(cat_rst))
# }

# Convert categorical variables to numeric in RasterStack
#' @importFrom raster deratify
.remove_cats <- function(rst_stack){
  # Check
  categ_vars <- names(rst_stack)[is.factor(rst_stack)]
  if (length(categ_vars) == 0) {
    rst_stack
  } else {
    for (nm in categ_vars) {
      rst_stack[[nm]] <- deratify(rst_stack[[nm]], complete = TRUE)
    }
    rst_stack
  }
}

## Get the mode of categorical integerish vector
#' @importFrom dplyr %>%
.mode <- function(v) {
  # Get mode
  summary(v) %>% sort(decreasing = T) %>%
    `[`(1) %>% names() %>% as.factor()
}

# Get the intersects of fitted curve and y = 0
.find_intersects <- function(fun, values, bins = 100) {
  vals_start <- seq(min(values), max(values),
                    ceiling((max(values) - min(values)) / bins))
  vals_end <- c(vals_start[-1], max(values))
  roots <- lapply(1:length(vals_start), function(i){
    inv <- c(vals_start[i], vals_end[i])
    if (!inherits(try(rt <- uniroot(fun, inv)$root, silent = TRUE),
                  "try-error")) rt
  })
  unlist(roots[!sapply(roots,is.null)])
}

## A simple function to calculate Boyce Index in Hirzel et al. 2006
## The objective of this code chunk is to avoid import heavy package
## the interested users can directly use their packages.
## Reference from https://github.com/adamlilith/enmSdm/blob/master/R/contBoyce.r
## and https://github.com/ecospat/ecospat/blob/master/ecospat/R/ecospat.boyce.R
## License: GNU General Public License v3.0
#' @importFrom stats cor
.cont_boyce <- function(obs, fit,
                       num_bins = 101,
                       bin_width = 0.1,
                       auto_window = TRUE,
                       method = 'spearman',
                       drop.zero = TRUE,
                       na.rm = FALSE) {
  # Check inputs
  checkmate::check_int(num_bins)
  checkmate::check_number(bin_width)
  checkmate::check_choice(method,
                          choices = c("pearson", "kendall", "spearman"))
  checkmate::check_logical(auto_window)
  checkmate::check_logical(drop.zero)
  checkmate::check_logical(na.rm)

  # if all NAs
  if (all(is.na(obs)) | all(is.na(fit))) return(NA)

  # [val1, val2) (assumes max value is > 0)
  residual <- .Machine$double.eps
  lowest <- ifelse(auto_window, min(c(obs, fit), na.rm = na.rm), 0)
  highest <- ifelse(auto_window,
                    max(c(obs, fit), na.rm = na.rm) + residual,
                    1 + residual)

  window_width <- bin_width * (highest - lowest)

  lows <- seq(lowest, highest - window_width, length.out = num_bins)
  highs <- seq(lowest + window_width + residual,
               highest, length.out = num_bins)

  ### tally proportion of test presences/background sites in each class
  freq_obs_bg <- do.call(rbind, lapply(1:num_bins, function(count_class) {
    # number of presence predictions in this class
    obs_in <- obs >= lows[count_class] & obs < highs[count_class]
    freq_obs <- sum(obs_in, na.rm = na.rm)

    # number of background predictions in this class
    bg_in <- fit >= lows[count_class] & fit < highs[count_class]
    freq_bg <- sum(bg_in, na.rm = na.rm)

    data.frame(freq_obs = freq_obs, freq_bg = freq_bg)
  }))
  freq_obs <- freq_obs_bg$freq_obs
  freq_bg <- freq_obs_bg$freq_bg
  rm(freq_obs_bg)

  # mean bin prediction
  mean_pred <- rowMeans(cbind(lows, highs))

  # Store original values for return
  HS <- mean_pred
  F.ratio <- (freq_obs / length(obs)) / (freq_bg / length(fit))

  # add small number to each bin that has 0 background frequency
  # but does have a presence frequency > 0
  if (any(freq_obs > 0 & freq_bg == 0)) {
    freq_bg[freq_obs > 0 & freq_bg == 0] <- 0.1
  }

  # remove classes with 0 presence frequency
  if (drop.zero && (0 %in% freq_obs)) {
    zeros <- which(freq_obs == 0)
    mean_pred[zeros] <- NA
    freq_obs[zeros] <- NA
    freq_bg[zeros] <- NA
  }

  # remove classes with 0 background frequency
  if (any(0 %in% freq_bg)) {
    zeros <- which(freq_bg == 0)
    mean_pred[zeros] <- NA
    freq_obs[zeros] <- NA
    freq_bg[zeros] <- NA
  }

  P <- freq_obs / length(obs)
  E <- freq_bg / length(fit)
  PE <- P / E

  # remove NAs
  ids <- !(is.na(mean_pred) | is.na(PE))
  mean_pred <- mean_pred[ids]
  PE <- PE[ids]

  # calculate continuous Boyce index (cbi)
  cbi <- cor(x = mean_pred, y = PE, method = method)

  return(list(F.ratio = F.ratio, cor = round(cbi, 3), HS = HS))
}

# A function to thin points to plot shap dependence
#' @importFrom rlang .data
#' @importFrom dplyr slice_sample ungroup
.shap_plot_thin <- function(x_trans,
                            sample_prop,
                            sample_bin,
                            seed,
                            category = FALSE) {
  if (category) {
    # Resample within each category
    set.seed(seed)
    x_trans %>%
      group_by(.data$variable, .data$x) %>%
      slice_sample(prop = sample_prop) %>%
      ungroup()
  } else {
    # Resample within bins
    x_trans_bins <- do.call(
      rbind, lapply(unique(x_trans$variable), function(var_select) {
        x_trans_sub <- x_trans %>% filter(.data$variable == var_select)
        width <- (max(x_trans_sub$x) - min(x_trans_sub$x)) / sample_bin
        x_trans_sub %>%
          mutate(bin = cut_width(.data$x, width = width, boundary = 0))
      }))

    set.seed(seed)
    x_trans_bins %>%
      group_by(.data$variable, .data$bin) %>%
      slice_sample(prop = sample_prop) %>%
      ungroup()
  }
}

# Function to sample background
#' @importFrom rlang .data
#' @importFrom dplyr select mutate sample_n
#' @importFrom stars st_rasterize
#' @importFrom sf st_as_sf
#' @importFrom stars st_xy2sfc
.bg_sampling <- function(rst_template, obs, seed, num) {
  if (is.null(obs)) {
    bg_samples <- rst_template
  } else {
    # Remove observations from background
    bg_samples <- obs %>% select() %>%
      mutate(value = 0) %>%
      st_rasterize(rst_template)
    bg_samples[bg_samples == 0] <- NA
  }

  # Sampling
  set.seed(seed)
  bg_samples %>%
    st_xy2sfc(as_points = T) %>%
    st_as_sf() %>%
    sample_n(min(num, nrow(.))) %>%
    select()
}

# TODO: let the model have high flexibility to take care of the
# input dataset.
# Function to mode the observations
#' @importFrom rlang .data
#' @importFrom dplyr select mutate sample_n filter
.obs_moding <- function(mode, obs_mode, rst_template,
                        obs, seed, contamination) {
  # Sampling
  if (mode == "normal") {
    if (obs_mode == "perfect_presence") {
      # Sample background
      obs <- .bg_sampling(
        rst_template, obs, seed, nrow(obs) * contamination) %>%
        mutate("observation" = 0) %>%
        rbind(obs)
      contamination <- sum(obs$observation == 0) /
        sum(obs$observation == 1)
    } else if (obs_mode == "presence_absence") {
      # Sample absence
      obs_frd <- obs %>% filter(.data$observation == 1)
      set.seed(seed)
      obs_frd <- obs_frd %>%
        rbind(obs %>% filter(.data$observation == 0) %>%
                sample_n(min((nrow(obs_frd) * contamination),
                             nrow(.))))
      obs <- obs_frd; rm(obs_frd)
      contamination <- sum(obs$observation == 0) /
        sum(obs$observation == 1)
    } else {
      contamination <- 0.001
    }
  } else if (mode == "reverse") {
    # Extract presences
    obs <- obs %>% filter(.data$observation == 1)
    # Sample background
    obs <- .bg_sampling(
      rst_template, obs, seed, nrow(obs) * (1 / contamination)) %>%
      mutate("observation" = 0) %>%
      rbind(obs)

    # Double check the contamination
    if (sum(obs$observation == 1) /
        sum(obs$observation == 0) > 0.3) {
      # Sample absence
      obs_frd <- obs %>% filter(.data$observation == 1)
      set.seed(seed)
      obs_frd <- obs_frd %>%
        rbind(obs %>% filter(.data$observation == 0) %>%
                sample_n((nrow(obs_frd) * 0.3)))
      obs <- obs_frd; rm(obs_frd)
      contamination <- 0.3
    } else {
      contamination <- sum(obs$observation == 1) /
        sum(obs$observation == 0)
    }

  } else { # joint
    if (obs_mode == "presence_absence") {
      ## For forward
      obs_frd <- obs %>% filter(.data$observation == 1)
      set.seed(seed)
      obs_frd <- obs_frd %>%
        rbind(obs %>% filter(.data$observation == 0) %>%
                sample_n((nrow(obs_frd) * contamination)))
      contamination_frd <- sum(obs_frd$observation == 0) /
        sum(obs_frd$observation == 1)

      ## For reverse
      obs_rev <- obs %>% filter(.data$observation == 0)
      set.seed(seed)
      obs_rev <- obs_rev %>%
        rbind(obs %>% filter(.data$observation == 1) %>%
                sample_n((nrow(obs_rev) * contamination)))
      contamination_rev <- sum(obs_rev$observation == 1) /
        sum(obs_rev$observation == 0)
    } else {
      # For forward, use presence
      if (obs_mode == "perfect_presence") {
        # Sample background
        obs_frd <- .bg_sampling(
          rst_template, obs, seed, nrow(obs) * contamination) %>%
          mutate("observation" = 0) %>%
          rbind(obs)
        contamination_frd <- sum(obs_frd$observation == 0) /
          sum(obs_frd$observation == 1)
      } else {
        obs_frd <- obs
        contamination_frd <- 0.001
      }

      # For reserve, use background
      # Extract presences
      obs <- obs %>% filter(.data$observation == 1)
      # Sample background
      obs_rev <- .bg_sampling(
        rst_template, obs, seed, nrow(obs) * (1 / contamination)) %>%
        mutate("observation" = 0) %>%
        rbind(obs)

      # Double check the contamination
      if (sum(obs_rev$observation == 1) /
          sum(obs_rev$observation == 0) > 0.3) {
        # Sample absence
        obs_rev_frd <- obs_rev %>% filter(.data$observation == 1)
        set.seed(seed)
        obs_rev_frd <- obs_rev_frd %>%
          rbind(obs_rev %>% filter(.data$observation == 0) %>%
                  sample_n((nrow(obs_rev_frd) * 0.3)))
        obs_rev <- obs_rev_frd; rm(obs_rev_frd)
        contamination <- rev <- 0.3
      } else {
        contamination_rev <- sum(obs_rev$observation == 1) /
          sum(obs_rev$observation == 0)
      }
    }
  }

  # Return the result according to mode
  if (mode == "joint") {
    return(list("obs_frd" = obs_frd,
                "obs_rev" = obs_rev,
                "contamination_frd" = contamination_frd,
                "contamination_rev" = contamination_rev))
  } else {
    return(list("obs" = obs,
                "contamination" = contamination))
  }
}
