rm(list=ls())
dispersion <- 0.1
nat <- 4

trials <- 10000
set.seed(10)
gamma_vec <- stats::rgamma(trials, shape = dispersion, scale = exp(nat)/dispersion)
poisson_vec <- stats::rpois(trials, lambda = gamma_vec)

set.seed(10)
nbinom_vec <- stats::rnbinom(trials, size = dispersion, mu = exp(nat))

plot(sort(poisson_vec), sort(nbinom_vec), asp = T)
mean(poisson_vec); sd(poisson_vec)
mean(nbinom_vec); sd(nbinom_vec)
