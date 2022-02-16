rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_layer4_processed2.RData")
source("initialization.R")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))

K <- min(50, round(min(dim(mat))*.5))
n <- nrow(mat)
p <- ncol(mat)

time_start1 <- Sys.time()
init_res <- initialize_esvd2(mat,
                             k = K,
                             family = "poisson",
                             covariates = covariates,
                             column_set_to_one = NULL,
                             offset_vec = rep(0, n),
                             verbose = 2)
time_end1 <- Sys.time()

save(date_of_run, session_info,
     sns, init_res, time_start1, time_end1,
     file = "../../../../out/Writeup10/Writeup10_sns_layer4_esvd2.RData")


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
     sns, init_res, time_start1, time_end1,
     esvd_res, time_start2, time_end2,
     file = "../../../../out/Writeup10/Writeup10_sns_layer4_esvd2.RData")

#################### 

print("Starting final fit, where library size coef can change")

covariates <- cbind(esvd_res$offset_vec, esvd_res$covariates)
colnames(covariates)[1] <- "Log_UMI"
b_mat <- cbind(1, init_res$b_mat)
colnames(b_mat)[1] <- "Log_UMI"

time_start3 <- Sys.time()
esvd_res_full <- eSVD2::opt_esvd(esvd_res$x_mat,
                                 esvd_res$y_mat,
                                 mat,
                                 family = "poisson",
                                 nuisance_param_vec = NA,
                                 library_size_vec = 1,
                                 method = "newton",
                                 b_init = b_mat,
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
     sns, init_res, time_start1, time_end1,
     esvd_res, time_start2, time_end2,
     esvd_res_full, time_start3, time_end3,
     file = "../../../../out/Writeup10/Writeup10_sns_layer4_esvd2.RData")


###########################

print("Starting nuisance parameter estimation")
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
diagnosis_idx <- which(colnames(esvd_res_full$covariates) == "diagnosis_ASD")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,diagnosis_idx,drop=F], esvd_res_full$b_mat[,diagnosis_idx,drop=F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,-diagnosis_idx],
  esvd_res_full$b_mat[,-diagnosis_idx]
))

time_start4 <- Sys.time()
nuisance_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')
  val <- tryCatch(eSVD2:::gamma_rate(x = mat[,j],
                                     mu = mean_mat_nolib[,j],
                                     s = library_mat[,j]),
                  error = function(c) 0)
  val
})
time_end4 <- Sys.time()

save(date_of_run, session_info,
     sns, init_res, time_start1, time_end1,
     esvd_res, time_start2, time_end2,
     esvd_res_full, time_start3, time_end3,
     nuisance_vec, time_start4, time_end4,
     file = "../../../../out/Writeup10/Writeup10_sns_layer4_esvd2.RData")

