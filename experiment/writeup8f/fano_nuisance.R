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
  res <- stats::optimize(calc.loglik.b, interval = c(1/max_val, 1/min_val),
                           maximum = F, y = y, mu = mu, sf = sf)
  # b <- res$minimum
  # b_loglik <- res$objective
  # min_loglik <- calc.loglik.b(1/max_val, y, mu, sf)
  # max_loglik <- calc.loglik.b(1/min_val, y, mu, sf)
  #
  # # determine if the estimate is numerically stable
  # if (abs(b_loglik - min_loglik) < 0.5) {
  #   if (abs(b_loglik-max_loglik) > 10) {
  #     res_max <- stats::uniroot(function(x) calc.loglik.b(x, y, mu, sf) + max_loglik - 10,
  #                      c(1/max_val, 1/min_val))$root
  #   } else {
  #     b.max <- 1/min_val
  #   }
  #   samp <- (exp((exp(stats::ppoints(100))-1)/2)-1)
  #   samp <- samp/max(samp)*(max_val - min_val) + min_val
  #   loglik <- mapply(calc.loglik.b, b = samp,
  #                    MoreArgs = list(y = y, mu = mu, sf = sf))
  #   loglik2 <- exp(-loglik-min(-loglik))
  #   loglik3 <- loglik2/sum(loglik2)
  #   b <- mean(sample(samp, 10000, replace = TRUE, prob = loglik3))
  #   b.loglik <- calc.loglik.b(b, y, mu, sf)
  # }
  return(c(1/res$minimum, res$objective))
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
