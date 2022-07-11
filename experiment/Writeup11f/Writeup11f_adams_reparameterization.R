rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/adams_T_esvd.RData")

covariate_mat <- eSVD_obj$covariates
covariate_mat <- covariate_mat[,-grep("Subject_Identity", colnames(covariate_mat))]
rank_val <- Matrix::rankMatrix(covariate_mat)
dim(covariate_mat)
rank_val

cor_mat <- stats::cor(covariate_mat, eSVD_obj$fit_Second$x_mat)
round(cor_mat, 2)

nat_mat1 <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat) + tcrossprod(eSVD_obj$covariates, eSVD_obj$fit_Second$z_mat)

x_mat <- eSVD_obj$fit_Second$x_mat
y_mat <- eSVD_obj$fit_Second$y_mat
z_mat <- eSVD_obj$fit_Second$z_mat
k <- ncol(x_mat)
p <- nrow(y_mat)
covariate_mat2 <- covariate_mat[,which(colnames(covariate_mat) != "Intercept")]
for(ell in 1:k){
  tmp_df <- cbind(x_mat[,ell], covariate_mat2)
  colnames(tmp_df)[1] <- "x"
  tmp_df <- as.data.frame(tmp_df)

  lm_res <- stats::lm(x ~ . , data = tmp_df)
  coef_vec <- stats::coef(lm_res)
  names(coef_vec)[1] <- "Intercept"
  print(paste0(ell, ": R2 of ", round(summary(lm_res)$r.squared, 2)))

  for(j in 1:p){
    z_mat[j,names(coef_vec)] <- z_mat[j,names(coef_vec)] + coef_vec*y_mat[j,ell]
  }

  x_mat[,ell] <- stats::residuals(lm_res)
}

nat_mat2 <- tcrossprod(x_mat, y_mat) + tcrossprod(eSVD_obj$covariates, z_mat)
sum(abs(nat_mat1 - nat_mat2))

#########

res <- eSVD2:::.reparameterize(x_mat, y_mat, equal_covariance = T)
x_mat <- res$x_mat; y_mat <- res$y_mat

cor_mat <- stats::cor(covariate_mat, x_mat)
print(round(cor_mat,2))
apply(x_mat, 2, quantile)
round(apply(z_mat, 2, quantile),2)
