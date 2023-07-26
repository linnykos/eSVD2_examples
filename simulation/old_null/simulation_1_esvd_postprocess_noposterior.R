rm(list=ls())
library(Seurat)

load("../../out/simulation/simulation_1_esvd.RData")

nat_mat1 <- tcrossprod(eSVD_obj[["fit_Second"]]$x_mat, eSVD_obj[["fit_Second"]]$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,"cc_1"], eSVD_obj[["fit_Second"]]$z_mat[,"cc_1"])
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
case_idx <- which(individual_vec %in% case_individuals)
control_idx <- which(individual_vec %in% control_individuals)

p <- ncol(mean_mat)
vanilla_pvalue <-  apply(mean_mat, 2, function(y){
  res <- stats::wilcox.test(x = y[case_idx],
                            y = y[control_idx])
  res$p.value
})
