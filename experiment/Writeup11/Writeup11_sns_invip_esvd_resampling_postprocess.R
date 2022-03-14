rm(list=ls())
library(Seurat)
load("../../../../out/Writeup10/Writeup10_sns_invip_processed.RData")
load("../../../../out/Writeup11/Writeup11_sns_invip_esvd_resample.RData")

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
case_control_variable = "diagnosis_ASD"
offset_var <- setdiff(colnames(esvd_init$covariates), case_control_variable)

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,case_control_variable,drop = F],
                       esvd_res_full$b_mat[,case_control_variable,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,offset_var],
  esvd_res_full$b_mat[,offset_var]
))

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
AplusAlpha <- mat + Alpha
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2
