# from https://github.com/mohuangx/SAVER/blob/master/R/optimize_variance.R
# and https://github.com/mohuangx/SAVER/blob/master/R/calc_loglik.R
calculate_fano_parameter <- function(y, mu, sf,
                                     max_val = stats::median(sf),
                                     min_val = 1) {
  stopifnot(length(y) == length(mu), length(y) == length(sf))

  n <- length(y)
  res <- stats::optimize(calc.loglik.b,
                         maximum = F,
                         interval = c(min_val, max_val),
                         y = y, mu = mu,
                         sf = sf)

  res$minimum
}

calc.loglik.b <- function(b, y, mu, sf) {
  func1 <- sum(mu*b*log(b))
  func2 <- -sum(lgamma(mu*b))
  func3 <- sum(lgamma(y+mu*b))
  func4 <- -sum((y+mu*b)*log(sf+b))
  return(-sum(func1, func2, func3, func4))
}
