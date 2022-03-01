rm(list=ls())
load("../../../../out/Writeup11/Writeup11_sns_invip_esvd_coef.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##################

print("Starting first eSVD fit")

case_control_variable = "diagnosis_ASD"
offset_var <- setdiff(colnames(esvd_init$covariates), case_control_variable)
offset_mat <- tcrossprod(esvd_init$covariates[,offset_var], esvd_init$b_mat[,offset_var])
covariate_init <- esvd_init$covariates[,case_control_variable,drop = F]
b_init <- esvd_init$b_mat[,case_control_variable,drop = F]

time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(esvd_init$x_mat,
                            esvd_init$y_mat,
                            mat,
                            family = "poisson",
                            nuisance_param_vec = NA,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = b_init,
                            covariates = covariate_init,
                            offset_vec = NULL,
                            offset_mat = offset_mat,
                            global_estimate = F,
                            l2pen = 0.01,
                            max_iter = 50,
                            run_cpp = F,
                            reparameterize = F,
                            reestimate_nuisance = F,
                            verbose = 1)
time_end2 <- Sys.time()

save(date_of_run, session_info,
     esvd_init,
     esvd_res, time_start2, time_end2,
     file = "../../../../out/Writeup11/Writeup11_sns_invip_esvd_fit.RData")

##################

print("Starting final fit, where library size coef can change")

covariates <- esvd_init$covariates
b_mat <- esvd_init$b_mat
b_mat[,case_control_variable] <- esvd_res$b_mat[,case_control_variable]

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
                                 offset_mat = NULL,
                                 global_estimate = F,
                                 l2pen = 0.01,
                                 max_iter = 50,
                                 run_cpp = T,
                                 reparameterize = F,
                                 reestimate_nuisance = F,
                                 verbose = 1)
time_end3 <- Sys.time()

save(date_of_run, session_info,
     esvd_init,
     esvd_res, time_start2, time_end2,
     esvd_res_full, time_start3, time_end3,
     file = "../../../../out/Writeup11/Writeup11_sns_invip_esvd_fit.RData")

###########

print("Starting nuisance parameter estimation")
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,case_control_variable,drop = F],
                       esvd_res_full$b_mat[,case_control_variable,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,offset_var],
  esvd_res_full$b_mat[,offset_var]
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
     esvd_init,
     esvd_res, time_start2, time_end2,
     esvd_res_full, time_start3, time_end3,
     nuisance_vec, time_start4, time_end4,
     file = "../../../../out/Writeup11/Writeup11_sns_invip_esvd_fit.RData")

print("Finished")
