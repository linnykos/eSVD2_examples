rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_layer23_processed.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))

K <- min(50, round(min(dim(mat))*.5))
n <- nrow(mat)
p1 <- ncol(mat)

covariates_nolib <- covariates[,which(colnames(covariates) != "Log_UMI")]

time_start1 <- Sys.time()
init_res <- initialize_esvd2(mat,
                             k = K,
                             family = "poisson",
                             covariates = covariates_nolib,
                             column_set_to_one = NULL,
                             offset_vec = covariates[,"Log_UMI"],
                             verbose = 1)
time_end1 <- Sys.time()

print("Starting large Poisson fit, with library coef fixed at 1")
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
                            global_estimate = F,
                            l2pen = 0.01,
                            max_iter = 50,
                            reparameterize = T,
                            reestimate_nuisance = F,
                            verbose = 1)
time_end2 <- Sys.time()

save(date_of_run, session_info,
     sns, esvd_res, time_start2, time_end2,
     file = "../../../../out/Writeup10/Writeup10_sns_layer23_esvd.RData")

#################

print("Starting final fit, where library size coef can change")
time_start3 <- Sys.time()
set.seed(10)
tmp_b <- matrix(1, nrow = nrow(esvd_res$b_mat), ncol = ncol(covariates))
colnames(tmp_b) <- colnames(covariates)
for(i in colnames(tmp_b)){
  idx <- which(colnames(esvd_res$b_mat) == i)
  if(length(idx) == 1){
    tmp_b[,i] <- esvd_res$b_mat[,idx]
  }
}

esvd_res_full <- eSVD2::opt_esvd(esvd_res$x_mat,
                                 esvd_res$y_mat,
                                 mat,
                                 family = "poisson",
                                 nuisance_param_vec = NA,
                                 library_size_vec = 1,
                                 method = "newton",
                                 b_init = tmp_b,
                                 covariates = covariates,
                                 offset_vec = rep(0, nrow(mat)),
                                 global_estimate = F,
                                 l2pen = 0.01,
                                 max_iter = 50,
                                 reparameterize = T,
                                 reestimate_nuisance = F,
                                 verbose = 1)
time_end3 <- Sys.time()

save(date_of_run, session_info,
     sns, esvd_res, time_start2, time_end2,
     esvd_res_full, time_start3, time_end3,
     file = "../../../../out/Writeup10/Writeup10_sns_layer23_esvd.RData")

