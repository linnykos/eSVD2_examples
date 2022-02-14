rm(list=ls())
library(Seurat)
load("../../out/Writeup10/Writeup10_sns_layer23_esvd.RData")
load("../../out/Writeup10/Writeup10_sns_layer23_esvd_nuisance.RData")

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

library_idx <- which(!colnames(esvd_res_full$covariates) %in% c("Intercept", "diagnosis_ASD"))
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,-library_idx], esvd_res_full$b_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,library_idx],
  esvd_res_full$b_mat[,library_idx]
))

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
AplusAlpha <- mat + Alpha
quantile(library_mat)
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2

###############

gene_names <- c("MT-ND1", "MT-ND2", "MT-CO2", "MT-CO3",
                "MT-ND3", "MT-ND4L", "MT-ND4", "MT-CYB")
idx <- which(colnames(mat) == gene_names[1])
plot(mean_mat[,idx], mat[,idx], asp = T)
hist(esvd_res_full$b_mat[,"percent.mt"])
rug(esvd_res_full$b_mat[gene_names,"percent.mt"], col = 2, lwd = 2)
hist(esvd_res_full$b_mat[,"diagnosis_ASD"])
rug(esvd_res_full$b_mat[gene_names,"diagnosis_ASD"], col = 2, lwd = 2)


plot(mean_mat_nolib[,idx], mat[,idx])
plot(posterior_mean_mat[,idx], mat[,idx])

metadata <- sns@meta.data
case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
case_idx <- which(metadata[,"diagnosis"] == "ASD")
control_idx <- which(metadata[,"diagnosis"] == "Control")

par(mfrow = c(1,2))
plot(posterior_mean_mat[case_idx,idx], mat[case_idx,idx], pch = 16,
     col = rgb(0.5,0.5,0.5,0.1),
     xlim = range(posterior_mean_mat[,idx]),
     ylim = range(mat[,idx]))
plot(posterior_mean_mat[control_idx,idx], mat[control_idx,idx], pch = 16,
     col = rgb(0.5,0.5,0.5,0.1),
     xlim = range(posterior_mean_mat[,idx]),
     ylim = range(mat[,idx]))

par(mfrow = c(1,2))
hist(sns$percent.mt[case_idx])
hist(sns$percent.mt[control_idx])


