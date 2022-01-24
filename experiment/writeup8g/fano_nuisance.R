# from https://github.com/mohuangx/SAVER/
calculate_fano_parameter <- function(y, mu, sf,
                                     max_val = stats::median(sf),
                                     min_val = min(sf)) {
  stopifnot(length(y) == length(mu), length(y) == length(sf))

  n <- length(y)

  # determine initial estimate
  res <- Rmpfr::optimizeR(calc.loglik.b, lower = 1/max_val,
                          upper = 1/min_val,
                          maximum = F, y = y, mu = mu, sf = sf,
                          maxiter = 100)

  1/res$minimum
}

calc.loglik.b <- function(b, y, mu, sf, digits = 5) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  func1 <- sum(Rmpfr::mpfr(mu/b*log(1/b), digits))
  func2 <- -sum(lgamma(Rmpfr::mpfr(mu/b, digits)))
  func3 <- sum(lgamma(Rmpfr::mpfr(y+mu/b, digits)))
  func4 <- -sum(Rmpfr::mpfr(y+mu/b, digits)*log(Rmpfr::mpfr(sf+1/b, digits)))
  return(-sum(func1, func2, func3, func4))
}
