rm(list=ls())
library(Seurat)
load("../../../../out/Writeup11/Writeup11_adams_ciliated_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(Matrix::t(adams[["RNA"]]@counts[adams[["RNA"]]@var.features,]))
case_control_variable <- "Disease_Identity_IPF"
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,case_control_variable,drop = F],
                       esvd_res_full$b_mat[,case_control_variable,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
offset_var <- setdiff(colnames(esvd_init$covariates), case_control_variable)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,offset_var],
  esvd_res_full$b_mat[,offset_var]
))

idx <- which(nuisance_vec == 0)
nuisance_error_list <- vector("list", length = 4)
names(nuisance_error_list) <- c("error", "nan", "large", "normal")
nuisance_error_list$error <- lapply(idx, function(j){
  data.frame(x = mat[,j],
             mu = mean_mat_nolib[,j],
             s = library_mat[,j])
})

j <- 8
eSVD2:::gamma_rate(x = nuisance_error_list$error[[j]][,"x"],
                   mu = nuisance_error_list$error[[j]][,"mu"],
                   s = nuisance_error_list$error[[j]][,"s"])

save(nuisance_error_list, session_info, date_of_run,
     file = "../../error_reports/2022-03-29-nuisance.RData")

##########

load("../../../../out/Writeup11/Writeup11_habermann_ciliated_esvd.RData")

mat <- as.matrix(Matrix::t(habermann[["RNA"]]@counts[habermann[["RNA"]]@var.features,]))
case_control_variable <- "Diagnosis_IPF"
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,case_control_variable,drop = F],
                       esvd_res_full$b_mat[,case_control_variable,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
offset_var <- setdiff(colnames(esvd_init$covariates), case_control_variable)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,offset_var],
  esvd_res_full$b_mat[,offset_var]
))

idx <- which(nuisance_vec == 0)[1:10]
tmp <- lapply(idx, function(j){
  data.frame(x = mat[,j],
             mu = mean_mat_nolib[,j],
             s = library_mat[,j])
})
nuisance_error_list$error <- c(nuisance_error_list$error, tmp)

idx <- which(is.na(nuisance_vec))
nuisance_error_list$na <- lapply(idx, function(j){
  data.frame(x = mat[,j],
             mu = mean_mat_nolib[,j],
             s = library_mat[,j])
})

idx <- which(nuisance_vec > 1000)
nuisance_error_list$large <- lapply(idx, function(j){
  data.frame(x = mat[,j],
             mu = mean_mat_nolib[,j],
             s = library_mat[,j])
})

save(nuisance_error_list, session_info, date_of_run,
     file = "../../error_reports/2022-03-29-nuisance.RData")

######################

val_vec <- c(0.01, 0.5, 1, 3, 10)
nuisance_vec[is.na(nuisance_vec)] <- 5000
nuisance_error_list$normal <- lapply(val_vec, function(val){
  j <- which.min(abs(val - nuisance_vec))
  data.frame(x = mat[,j],
             mu = mean_mat_nolib[,j],
             s = library_mat[,j])
})
quantile(nuisance_error_list$normal[[1]]$x)
names(nuisance_error_list$normal) <- paste0("value_close_to_", val_vec)

nuisance_error_list$error <- nuisance_error_list$error[c(1:8,10:18)]

save(nuisance_error_list, session_info, date_of_run,
     file = "../../error_reports/2022-03-29-nuisance.RData")

##############################3

sapply(nuisance_error_list, length)
k <- 4
j <- 4
eSVD2:::gamma_rate(x = nuisance_error_list[[k]][[j]][,"x"],
                   mu = nuisance_error_list[[k]][[j]][,"mu"],
                   s = nuisance_error_list[[k]][[j]][,"s"])


