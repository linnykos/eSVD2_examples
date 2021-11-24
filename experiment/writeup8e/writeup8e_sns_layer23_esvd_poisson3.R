rm(list=ls())
load("../../../../out/writeup8e/writeup8e_sns_layer23_processed_onlymat.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

K <- 30
n <- nrow(mat)
p <- ncol(mat)

time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat,
                                   k = K,
                                   family = "poisson",
                                   covariates = covariates,
                                   column_set_to_one = NULL,
                                   offset_vec = rep(0, nrow(mat)),
                                   verbose = 1)
time_end1 <- Sys.time()

print("Starting Poisson fit")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
                            init_res$y_mat,
                            mat,
                            family = "poisson",
                            nuisance_param_vec = NA,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = init_res$b_mat,
                            covariates = init_res$covariates,
                            offset_vec = init_res$offset_vec,
                            reestimate_nuisance = F,
                            global_estimate = F,
                            reparameterize = F,
                            max_iter = 50,
                            l2pen = 0.01,
                            verbose = 1)
time_end2 <- Sys.time()

save.image("../../../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson3.RData")

quantile(esvd_res$b_mat[,"Log_UMI"])
quantile(esvd_res$b_mat[,"diagnosis_ASD"])

