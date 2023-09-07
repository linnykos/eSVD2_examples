rm(list=ls())
library(Seurat)
library(eSVD2)
library(npregfast)

load("../../../out/main/habermann_T_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

lib_vec <- eSVD_obj$covariates[,"Log_UMI"]
dat <- as.matrix(eSVD_obj$dat)
p <- ncol(dat)
before_npreg_list <- sapply(1:p, function(j){
  print(j)

  tmp_df <- data.frame(y = dat[,j], x = lib_vec)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- lib_vec
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})

save(before_npreg_list, lib_vec,
     date_of_run, session_info,
     file = "../../../out/main/habermann_T_esvd_illustration-npreg.RData")

nat_mat1 <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat)
case_control_var <- eSVD2:::.get_object(eSVD_obj = eSVD_obj, what_obj = "init_case_control_variable", which_fit = "param")
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,case_control_var], eSVD_obj$fit_Second$z_mat[,case_control_var])
mean_mat_nolib <- exp(nat_mat1 + nat_mat2)
after_npreg_list <- sapply(1:p, function(j){
  print(j)

  tmp_df <- data.frame(y = mean_mat_nolib[,j], x = lib_vec)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- lib_vec
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})

save(before_npreg_list, after_npreg_list, lib_vec,
     date_of_run, session_info,
     file = "../../../out/main/habermann_T_esvd_illustration-npreg.RData")

