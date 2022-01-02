rm(list=ls())
load("../../../../out/writeup8f/writeup8f_pseudoreal_data.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

K <- min(50, round(min(dim(mat))*.75))
n <- nrow(mat)
p <- ncol(mat)

time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat,
                                   k = K,
                                   family = "poisson",
                                   covariates = covariates,
                                   column_set_to_one = "Log_UMI",
                                   offset_vec = rep(0, nrow(mat)),
                                   verbose = 1)
time_end1 <- Sys.time()

print("Starting large Poisson fit")
time_start3 <- Sys.time()
set.seed(10)
esvd_res_full <- eSVD2::opt_esvd(init_res$x_mat,
                                 init_res$y_mat,
                                 mat,
                                 family = "poisson",
                                 nuisance_param_vec = NA,
                                 library_size_vec = 1,
                                 method = "newton",
                                 b_init = init_res$b_mat,
                                 covariates = init_res$covariates,
                                 offset_vec = init_res$offset_vec,
                                 global_estimate = F,
                                 l2pen = 0.01,
                                 max_iter = 50,
                                 reparameterize = T,
                                 reestimate_nuisance = F,
                                 verbose = 1)
time_end3 <- Sys.time()

true_objects = list(autism_gene_idx = autism_gene_idx,
                    up_idx = up_idx, down_idx = down_idx,
                    nat_mat = nat_mat, lambda_mat = lambda_mat,
                    true_esvd = true_esvd)

save(date_of_run, session_info,
     mat, esvd_res_full, time_start3, time_end3,
     covariates, metadata,
     true_objects,
     file = "../../../../out/writeup8f/writeup8f_sns_pseudoreal_esvd4.RData")

