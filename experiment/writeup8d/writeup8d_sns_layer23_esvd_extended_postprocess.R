rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended.RData")

nat_mat1 <- tcrossprod(esvd_res_nb2$x_mat, esvd_res_nb2$y_mat)
nat_mat2 <- tcrossprod(esvd_res_nb2$covariates, esvd_res_nb2$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)
quantile(angle_vec[gene_idx])

###################

nat_mat1 <- tcrossprod(esvd_res_nb1$x_mat, esvd_res_nb1$y_mat)
nat_mat2 <- tcrossprod(esvd_res_nb1$covariates, esvd_res_nb1$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)
quantile(angle_vec[gene_idx])

cor_vec <- sapply(1:p, function(j){
  stats::cor(mat[,j], mean_mat[,j])
})
quantile(cor_vec)
quantile(cor_vec[gene_idx])
quantile(cor_vec[which(colnames(mat) %in% de_genes)])

cor_vec_nonzero <- sapply(1:p, function(j){
  idx <- which(mat[,j] != 0)
  stats::cor(mat[idx,j], mean_mat[idx,j])
})
quantile(cor_vec_nonzero)
quantile(cor_vec_nonzero[gene_idx])
quantile(cor_vec_nonzero[which(colnames(mat) %in% de_genes)])

####################

round(quantile(esvd_res_nb2$b_mat[,"diagnosis_ASD"]),2)
round(quantile(esvd_res_nb2$b_mat[de_genes,"diagnosis_ASD"]),2)
round(quantile(esvd_res_nb2$b_mat[gene_idx,"diagnosis_ASD"]),2)

round(quantile(esvd_res_nb1$b_mat[,"diagnosis_ASD"]),2)
round(quantile(esvd_res_nb1$b_mat[de_genes,"diagnosis_ASD"]),2)
round(quantile(esvd_res_nb1$b_mat[gene_idx,"diagnosis_ASD"]),2)

