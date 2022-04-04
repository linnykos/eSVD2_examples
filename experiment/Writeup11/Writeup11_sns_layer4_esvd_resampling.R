rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_layer4_processed2.RData")

library(Seurat)
library(eSVD2)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

replication_total <- 10
mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))

covariate_dat <- sns@meta.data[,c("percent.mt", "individual", "region", "age", "sex",
                                  "RNA.Integrity.Number", "post.mortem.hours",
                                  "diagnosis", "Seqbatch")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"individual"] <- as.factor(covariate_df[,"individual"])
covariate_df[,"region"] <- as.factor(covariate_df[,"region"])
covariate_df[,"diagnosis"] <- factor(covariate_df[,"diagnosis"], levels = c("Control", "ASD"))
covariate_df[,"sex"] <- as.factor(covariate_df[,"sex"])
covariate_df[,"Seqbatch"] <- as.factor(covariate_df[,"Seqbatch"])

idx_list <- lapply(unique(covariate_df[,"individual"]), function(indiv){
  which(covariate_df[,"individual"] == indiv)
})
names(idx_list) <- unique(covariate_df[,"individual"])

for(replication_idx in 1:replication_total){
  print("======")
  print(paste0("Starting replication index ", replication_idx))
  print("======")

  covariate_df[,"diagnosis"] <- factor(sns@meta.data[,"diagnosis"], levels = c("Control", "ASD"))
  indiv_status <- t(sapply(unique(covariate_df[,"individual"]), function(indiv){
    diagnosis <- covariate_df[which(covariate_df[,"individual"] == indiv)[1],"diagnosis"]
    c(as.character(indiv), as.character(diagnosis))
  }))
  set.seed(10*replication_idx)
  indiv_status[,2] <- sample(indiv_status[,2])
  diagnosis_vec <- as.character(covariate_df[,"diagnosis"])
  for(i in 1:length(idx_list)){
    diagnosis_vec[idx_list[[i]]] <- indiv_status[i,2]
  }
  covariate_df[,"diagnosis"] <- factor(diagnosis_vec, levels = c("Control", "ASD"))

  covariates <- eSVD2:::format_covariates(dat = mat,
                                          covariate_df = covariate_df,
                                          mixed_effect_variables = c("individual", "Seqbatch"))

  #####################
  mixed_effect_variables <- c(colnames(covariates)[grep("^individual", colnames(covariates))],
                              colnames(covariates)[grep("^Seqbatch", colnames(covariates))])

  time_start1 <- Sys.time()
  esvd_init <- eSVD2:::initialize_esvd(dat = mat,
                                       covariates = covariates,
                                       case_control_variable = "diagnosis_ASD",
                                       k = 30,
                                       lambda = 0.1,
                                       mixed_effect_variables = mixed_effect_variables,
                                       offset_variables = "Log_UMI",
                                       verbose = 2,
                                       tmp_path = paste0("../../../../out/Writeup11/Writeup11_sns_layer4_esvd_resample_tmp_",
                                                         replication_idx, ".RData"))
  time_end1 <- Sys.time()

  save(date_of_run, session_info, sns, covariate_df,
       esvd_init, time_start1, time_end1,
       file = paste0("../../../../out/Writeup11/Writeup11_sns_layer4_esvd_resample_", replication_idx, ".RData"))

  #############

  print("Starting first eSVD fit")

  mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
  case_control_variable <- "diagnosis_ASD"
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
                              l2pen = 0.1,
                              max_iter = 50,
                              run_cpp = F,
                              reparameterize = F,
                              reestimate_nuisance = F,
                              verbose = 1)
  time_end2 <- Sys.time()

  save(date_of_run, session_info, sns, covariate_df,
       esvd_init, time_start1, time_end1,
       esvd_res, time_start2, time_end2,
       file = paste0("../../../../out/Writeup11/Writeup11_sns_layer4_esvd_resample_", replication_idx, ".RData"))

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
                                   l2pen = 0.1,
                                   max_iter = 50,
                                   run_cpp = T,
                                   reparameterize = F,
                                   reestimate_nuisance = F,
                                   verbose = 1)
  time_end3 <- Sys.time()

  save(date_of_run, session_info, sns, covariate_df,
       esvd_init, time_start1, time_end1,
       esvd_res, time_start2, time_end2,
       esvd_res_full, time_start3, time_end3,
       file = paste0("../../../../out/Writeup11/Writeup11_sns_layer4_esvd_resample_", replication_idx, ".RData"))

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

  save(date_of_run, session_info, sns, covariate_df,
       esvd_init, time_start1, time_end1,
       esvd_res, time_start2, time_end2,
       esvd_res_full, time_start3, time_end3,
       nuisance_vec, time_start4, time_end4,
       file = paste0("../../../../out/Writeup11/Writeup11_sns_layer4_esvd_resample_", replication_idx, ".RData"))

  print("Finished")
}

