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
      rst_stack[[nm]] <- deratify(rst_stack[[nm]])
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
