# from https://github.com/mohuangx/SAVER/
calculate_fano_parameter <- function(y, mu, sf,
                                     max_val = 1e5,
                                     min_val = mean(y/sf)/var(y/sf)) {
  stopifnot(length(y) == length(mu), length(y) == length(sf))

  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  if (sum(mu) == 0) {
    return(c(0, 0))
  }

  # determine initial estimate
  b.vec <- stats::optimize(calc.loglik.b, interval = c(1/max_val, 1/min_val),
                           y = y, mu = mu, sf = sf)
  b <- b.vec$minimum
  b.loglik <- b.vec$objective
  min.b <- -calc.loglik.b(1/max_val, y, mu, sf)
  mle.b <- -calc.loglik.b(b, y, mu, sf)
  max.b <- -calc.loglik.b(1/min_val, y, mu, sf)

  # determine if the estimate is numerically stable
  if (mle.b - min.b < 0.5) {
    if (mle.b-max.b > 10) {
      b.max <- stats::uniroot(function(x) calc.loglik.b(x, y, mu, sf) + mle.b - 10,
                       c(1/max_val, 1/min_val))$root
    } else {
      b.max <- 1/min_val
    }
    samp <- (exp((exp(stats::ppoints(100))-1)/2)-1)
    samp <- samp/max(samp)*(max_val - min_val) + min_val
    loglik <- mapply(calc.loglik.b, b = samp,
                     MoreArgs = list(y = y, mu = mu, sf = sf))
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    b <- mean(sample(samp, 10000, replace = TRUE, prob = loglik3))
    b.loglik <- calc.loglik.b(b, y, mu, sf)
  }
  return(c(1/b, b.loglik))
}

calc.loglik.b <- function(b, y, mu, sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  func1 <- sum(mu/b*log(1/b))
  func2 <- -sum(lgamma(mu/b))
  func3 <- sum(lgamma(y+mu/b))
  func4 <- -sum((y+mu/b)*log(sf+1/b))
  return(-sum(func1, func2, func3, func4))
}
