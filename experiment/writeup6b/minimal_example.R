# x_mat <- init$x_mat; y_mat <- init$y_mat; b_mat <- init$b_mat
# save(mat, x_mat, y_mat, b_mat, covariates, nuisance_vec, file = "../../../../out/writeup6b/tmp.RData")

rm(list=ls())
load("../../../../out/writeup6b/tmp.Rdata")
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(x_mat, y_mat, mat, family = "neg_binom",
                             nuisance_param_vec = nuisance_vec, library_size_vec = 1,
                             b_init = b_mat, covariates = covariates,
                             max_iter = 100, verbose = 1)
