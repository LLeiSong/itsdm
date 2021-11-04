# A few internal functions used in this package
# Calculate the mean of absolute values in a vector
.abs_mean <- function(v) mean(abs(v))

# Get spatial resolution of stars
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
.auc_ratio <- function(occ, full) {
  roc_r <- .roc_ratio(occ, full)
  roc_r <- roc_r %>% arrange(cell)
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

# Function to get min and max values of a single band stars
.min_value <- function(x) min(x[[1]], na.rm = T)
.max_value <- function(x) max(x[[1]], na.rm = T)

# Linear stretch of single-band stars
.stars_stretch <- function(x, # stars
                           new_values = NULL,
                           minv = 0,
                           maxv = 1,
                           minq = 0.5,
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
    q <- cbind(.min_value(x), .max_value(x))
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
