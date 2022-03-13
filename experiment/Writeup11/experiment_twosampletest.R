rm(list=ls())
source("../eSVD2_examples/experiment/Writeup11/twosampletest.R")

set.seed(10)
n_each <- 15
factor_vec <- factor(c(rep("a", n_each), rep("b", n_each)))
n <- length(factor_vec)

trials <- 1000
pval_vec <- sapply(1:trials, function(i){
  set.seed(i)
  var_vec <- rep(1, n)
  mean_vec <- stats::rnorm(n)
  compute_twosample_pvalue(factor_vec = factor_vec,
                           mean_vec = mean_vec,
                           var_vec = var_vec,
                           verbose = 0)
})

plot(sort(pval_vec), seq(0, 1, length.out = length(pval_vec)), asp = T)

#############################

trials <- 1000
pval_vec <- sapply(1:trials, function(i){
  set.seed(i)
  var_vec <- rep(1, n)
  mean_vec <- c(stats::rnorm(n/2), stats::rnorm(n/2)+2)
  compute_twosample_pvalue(factor_vec = factor_vec,
                           mean_vec = mean_vec,
                           var_vec = var_vec,
                           verbose = 0)
})

plot(sort(pval_vec), seq(0, 1, length.out = length(pval_vec)), asp = T)

##############

trials <- 1000
pval_vec <- sapply(1:trials, function(i){
  set.seed(i)
  var_vec <- c(stats::runif(n/2), stats::runif(n/2)+1)
  mean_vec <- rep(0, n)
  compute_twosample_pvalue(factor_vec = factor_vec,
                           mean_vec = mean_vec,
                           var_vec = var_vec,
                           verbose = 0)
})

plot(sort(pval_vec), seq(0, 1, length.out = length(pval_vec)), asp = T)

##############

trials <- 1000
pval_vec <- sapply(1:trials, function(i){
  set.seed(i)
  var_vec <- rep(1, n)
  mean_vec <- c(stats::rnorm(n/2), stats::rnorm(7)-2, stats::rnorm(8)+2)
  compute_twosample_pvalue(factor_vec = factor_vec,
                           mean_vec = mean_vec,
                           var_vec = var_vec,
                           verbose = 0)
})

plot(sort(pval_vec), seq(0, 1, length.out = length(pval_vec)), asp = T)



