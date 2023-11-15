rm(list=ls())
library(Seurat)
library(eSVD2)
library(npregfast)

load("../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

lib_vec <- eSVD_obj$covariates[,"Log_UMI"]
dat <- Matrix::t(sns[["RNA"]]@data[Seurat::VariableFeatures(sns),])
p <- ncol(dat)
lognorm_npreg_list <- sapply(1:p, function(j){
  print(j)

  tmp_df <- data.frame(y = dat[,j], x = lib_vec)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- lib_vec
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})

save(before_npreg_list, after_npreg_list, lib_vec,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_esvd_illustration-npreg2.RData")

