rm(list=ls())

load("../../out/tmp.RData")

teststat_vec <- gaussian_teststat
observed_quantile <- c(0.05,0.95)

z_vec <- teststat_vec
fn <- function(param_vec,
               lb,
               N,
               N0,
               z0_vec,
               ub){
  delta0 <- param_vec[1]
  sigma0 <- param_vec[2]
  theta <-  param_vec[3]
  if(theta > 0.99 | theta < 0.01) return(Inf)
  denom <- stats::pnorm(ub, mean = delta0, sd = sigma0) - stats::pnorm(lb, mean = delta0, sd = sigma0)

  # compute each sample's likelihood
  sample_llvec <- sapply(z0_vec, function(z){
    stats::dnorm(z, mean = delta0, sd = sigma0, log = T)
  })

  # compute full log-likelihood, Equation 4.12
  loglik <- N0*log(theta) + (N-N0)*log(1-theta) + sum(sample_llvec) - N0*log(denom)
  -loglik
}

N <- length(z_vec)
tmp <- stats::quantile(z_vec, probs = observed_quantile)
lb <- tmp[1]; ub <- tmp[2]
idx <- intersect(which(z_vec >= lb), which(z_vec <= ub))
N0 <- length(idx)
z0_vec <- z_vec[idx]

init_delta0 <- mean(z0_vec)
init_sigma0 <- stats::sd(z0_vec)
init_theta <- 0.95*(stats::pnorm(ub, mean = init_delta0, sd = init_sigma0) - stats::pnorm(lb, mean = init_delta0, sd = init_sigma0))
res <- stats::optim(par = c(init_delta0, init_sigma0, init_theta),
                    fn = fn,
                    method = "Nelder-Mead",
                    lb = lb,
                    N = N,
                    N0 = N0,
                    z0_vec = z0_vec,
                    ub = ub)
res$par

###################################
null_mean2 <- mean(gaussian_teststat[-c(1:10)])
null_sd2 <- sd(gaussian_teststat[-c(1:10)])

null_mean2; null_sd2
