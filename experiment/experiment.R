rm(list=ls())

mu = 0.05
gamma = 1
size = 1000
n = 10000

set.seed(10)
lambda <- stats::rgamma(n, shape = gamma*mu, rate = gamma)
x <- stats::rpois(n, lambda = size*lambda)
par(mar = c(4, 0.5, 0.5, 0.5))
hist(x, breaks = 50)
