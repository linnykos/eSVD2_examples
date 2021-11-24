rm(list=ls())
load("../../../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson.RData")

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)
quantile(angle_vec[gene_idx])
