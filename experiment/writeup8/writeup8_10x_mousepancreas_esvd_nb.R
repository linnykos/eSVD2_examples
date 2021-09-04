rm(list=ls())

load("../../../../out/writeup6b/writeup6b_10x_mousepancreas_esvd.RData")

library(eSVD2)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
print("Starting NB nuisance")
time_start3 <- Sys.time()
nuisance_vec <- eSVD2::initialize_nuisance_param(mat, nat_mat, family = "neg_binom",
                                                 library_size_vec = 1)
time_end3 <- Sys.time()
rm(list = "nat_mat")
save.image("../../../../out/writeup8/writeup8_10x_mousepancreas_esvd_nb.RData")

print("Starting NB initialization")
time_start4 <- Sys.time()
set.seed(10)
init2 <- eSVD2::initialize_esvd(mat, k = K, family = "neg_binom",
                                nuisance_param_vec = nuisance_vec,
                                library_size_vec = 1,
                                covariates = covariates,
                                check_rank = F,
                                config = eSVD2::initialization_options(), verbose = 1)
time_end4 <- Sys.time()
save.image("../../../../out/writeup8/writeup8_10x_mousepancreas_esvd_nb.RData")

print("Starting NB estimation")
time_start5 <- Sys.time()
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(init2$x_mat, init2$y_mat, mat, family = "neg_binom",
                             nuisance_param_vec = nuisance_vec,
                             library_size_vec = 1,
                             b_init = init2$b_mat,
                             covariates = covariates,
                             max_iter = 100,
                             verbose = 1)
time_end5 <- Sys.time()
print("Finished")
save.image("../../../../out/writeup8/writeup8_10x_mousepancreas_esvd_nb.RData")

