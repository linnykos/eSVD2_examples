rm(list=ls())
library(eSVD2)

# example for negative binomial
set.seed(5)
nat_mat <- matrix(-1/5, 100, 5)
dat <- generate_data(nat_mat, family = "neg_binom", nuisance_param_vec = 10, library_size_vec = 1)

res <- initialize_esvd(dat, k = 1, family = "poisson", library_size_vec = 1)
res2 <- opt_esvd(res$x_mat, res$y_mat, dat, family = "poisson",
                 library_size_vec = 1)

nuisance_vec <- initialize_nuisance_param(dat, res2$x %*% t(res2$y), family = "neg_binom",
                                          library_size_vec = res$library_size_vec)
nuisance_vec

# example for curved gaussian
set.seed(5)
nat_mat <- matrix(1/5, 100, 5)
dat <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2, library_size_vec = 1)

res <- initialize_esvd(dat, k = 1, family = "exponential", library_size_vec = 1)
res2 <- opt_esvd(res$x_mat, res$y_mat, dat, family = "exponential",
                 library_size_vec = 1)

nuisance_vec <- initialize_nuisance_param(dat, res2$x %*% t(res2$y), family = "curved_gaussian",
                                          library_size_vec = res$library_size_vec)
nuisance_vec
