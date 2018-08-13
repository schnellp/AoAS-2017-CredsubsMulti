quantile.band <- function(sample, alpha, verbose=FALSE, low.mem=FALSE) {
  # Constructs quantile-based simultaneous credible bands from
  # a sample from the joint distribution of the variables.
  #
  # Args:
  #   sample: A matrix whose rows are draws from the
  #           joint distribution
  #   alpha: the (simultaneous) level of the credible band
  #   verbose: (logical) print status updates
  #   low.mem: (logical) use slower, low memory operations
  #
  # Returns:
  #   A list containing:
  #     upper: a vector corresponding to the upper bounds,
  #     lower: a vector corresponding to the lower bounds,
  #     alpha: the significance level,
  #     W.crit: the critical value (adjusted alpha / 2).
  #   If the columns of the sample are named, those names
  #   are retained in the bounding vectors.
  
  if (low.mem) {
    return(quantile.band.low.mem(sample, alpha, verbose))
  }
  
  sample <- as.matrix(sample)
  
  # 1-G_{Y}(y) = F_{-Y}(-y)
  double.sample <- cbind(sample, -sample)
  
  p <- ncol(double.sample)
  
  if (verbose) {
    print("Computing empirical CDF...")
  }
  
  mecdf <- apply(double.sample, 2, ecdf)
  
  if (verbose) {
    print("Finding tails of distribution...")
  }
  
  tails <- apply(as.matrix(1:p), 1, function(i) {
    mecdf[[i]](double.sample[, i])
  })
  
  # reclaim a bit of memory
  rm(mecdf, double.sample)
  
  if (verbose) {
    print("Computing suprema...")
  }
  
  W <- do.call(pmin, as.data.frame(tails))
  
  # type=1 uses the true (non-interpolated) definition,
  # the default interpolation (type=7) yields anticonservative bands
  W.crit <- quantile(W, alpha, type=1)
  
  if (verbose) {
    print("Computing bounds...")
  }
  
  upper <- -apply(-sample, 2, quantile, prob=W.crit, type=1)
  lower <-  apply( sample, 2, quantile, prob=W.crit, type=1)
  
  names(upper) <- names(lower) <- colnames(sample)
  
  list(upper=upper, lower=lower, alpha=alpha, W.crit=W.crit)
}

quantile.band.low.mem <- function(sample, alpha, verbose=FALSE) {
  # Constructs quantile-based simultaneous credible bands from
  # a sample from the joint distribution of the variables.
  # Uses less memory than quantile.band but is slower when
  # enough memory is available for quantile.band.
  #
  # Args:
  #   sample: A matrix whose rows are draws from the
  #           joint distribution
  #   alpha: the (simultaneous) level of the credible band
  #   verbose: (logical) print status updates
  #
  # Returns:
  #   A list containing:
  #     upper: a vector corresponding to the upper bounds,
  #     lower: a vector corresponding to the lower bounds,
  #     alpha: the significance level,
  #     W.crit: the critical value (adjusted alpha / 2).
  #   If the columns of the sample are named, those names
  #   are retained in the bounding vectors.
  
  sample <- as.matrix(sample)
  
  # 1-G_{Y}(y) = F_{-Y}(-y)
  double.sample <- cbind(sample, -sample)
  
  n <- nrow(double.sample)
  p <- ncol(double.sample)
  
  if (verbose) {
    print("Finding tails of distribution...")
  }
  
  for (j in 1:p) {
    double.sample[, j] <- ecdf(double.sample[, j])(double.sample[, j])
  }
  
  if (verbose) {
    print("Computing suprema...")
  }
  
  W <- numeric(n)
  
  for (i in 1:n) {
    W[i] <- min(double.sample[i, ])
  }
  
  rm(double.sample)
  
  # type=1 uses the true (non-interpolated) definition,
  # the default interpolation (type=7) yields anticonservative bands
  W.crit <- quantile(W, alpha, type=1)
  
  upper <- -apply(-sample, 2, quantile, prob=W.crit, type=1)
  lower <-  apply( sample, 2, quantile, prob=W.crit, type=1)
  
  names(upper) <- names(lower) <- colnames(sample)
  
  list(upper=upper, lower=lower, alpha=alpha, W.crit=W.crit)
}

multi.quantile.band <- function(sample, alpha, verbose=FALSE, low.mem=FALSE) {
  # A wrapper to quantile.band for vector-valued functions.
  #
  # Args:
  #   sample: An array of parameter values with at least two dimensions.
  #           Dims:
  #             - First: samples.
  #             - Second and (optionally) beyond: multivariate coordinates
  #   alpha: the (simultaneous) level of the credible band
  #   verbose: (logical) print status updates
  #
  # Returns:
  #   A list containing:
  #     upper: a vector corresponding to the upper bounds,
  #     lower: a vector corresponding to the lower bounds,
  #     alpha: the significance level,
  #     W.crit: the critical value (adjusted alpha / 2).
  
  n.draws <- dim(sample)[1]
  n.coords <- prod(dim(sample)[-1])
  
  sample.MAT <- sample
  dim(sample.MAT) <- c(n.draws, n.coords)
  
  sim.cred.band <- quantile.band(sample.MAT, alpha=alpha,
                                 verbose=verbose, low.mem=low.mem)
  
  multi.sim.cred.band <- sim.cred.band
  multi.sim.cred.band$upper <- array(sim.cred.band$upper, dim=dim(sample)[-1])
  multi.sim.cred.band$lower <- array(sim.cred.band$lower, dim=dim(sample)[-1])
  multi.sim.cred.band
}
