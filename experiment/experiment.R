rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_invip_esvd.RData")
eSVD_obj2 <- eSVD_obj

load("../../../out/Writeup11f/Writeup11f_sns_invip_esvd.RData")

colnames(eSVD_obj2$covariates)[which(!colnames(eSVD_obj2$covariates) %in% colnames(eSVD_obj$covariates))]
colnames(eSVD_obj$covariates)[which(!colnames(eSVD_obj$covariates) %in% colnames(eSVD_obj2$covariates))]

var_names <- colnames(eSVD_obj$covariates)[which(colnames(eSVD_obj$covariates) %in% colnames(eSVD_obj2$covariates))]
for(var_name in var_names){
  print(var_name)
  print(sum(abs(eSVD_obj$covariates[,var_name] - eSVD_obj2$covariates[,var_name])))
}
