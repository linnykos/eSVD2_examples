plot(eSVD_obj$fit_First$z_mat[,"covariate_1"])

nat_mat1 <- tcrossprod(eSVD_obj$fit_First$x_mat, eSVD_obj$fit_First$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates, eSVD_obj$fit_First$z_mat)
pred_mat <- exp(nat_mat1 + nat_mat2)
image(pred_mat)

esvd_res <- eSVD_obj$fit_First
covariates <- eSVD_obj$covariates
library_size_variable <- "Log_UMI"
library_idx <- which(colnames(covariates) %in% c("Intercept", library_size_variable))
nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(covariates[,-library_idx],
                       esvd_res$z_mat[,-library_idx])
# image(nat_mat1 + nat_mat2)

pred_mat <- exp(nat_mat1 + nat_mat2)
image(pred_mat)
image(log(pred_mat))

plot(esvd_res$z_mat[,"covariate_1"])
plot(esvd_res$z_mat[,"covariate_2"])
plot(esvd_res$z_mat[,"covariate_3"])
plot(esvd_res$z_mat[,"covariate_4"])
plot(esvd_res$z_mat[,"Intercept"])
plot(esvd_res$z_mat[,"Log_UMI"])

esvd_res <- eSVD_obj$fit_Init
plot(eSVD_obj$initial_Reg$z_mat1[,"covariate_1"])
plot(eSVD_obj$initial_Reg$log_pval)
plot(esvd_res$z_mat[,"covariate_1"])
plot(esvd_res$z_mat[,"covariate_2"])
plot(esvd_res$z_mat[,"covariate_3"])
plot(esvd_res$z_mat[,"covariate_4"])
plot(esvd_res$z_mat[,"Intercept"])
plot(esvd_res$z_mat[,"Log_UMI"])


###################

plot(eSVD_obj$fit_First$nuisance_vec, jitter(nuisance_vec))

zz <- eSVD_obj$fit_First$posterior_mean_mat
plot(zz[,1])
plot(zz[,100])
plot(zz[,150])
plot(zz[,90])

image(eSVD_obj$fit_First$posterior_mean_mat)
image(eSVD_obj$fit_First$posterior_var_mat)
image(log(eSVD_obj$fit_First$posterior_mean_mat))
image(log(eSVD_obj$fit_First$posterior_var_mat))
zz <- eSVD_obj$fit_First$posterior_mean_mat
zz <- pmin(zz, 1)
image(zz)

eSVD_obj2 <- compute_posterior(input_obj = eSVD_obj,
                               bool_adjust_covariates = F)
image(eSVD_obj2$fit_First$posterior_mean_mat)
zz <- eSVD_obj2$fit_First$posterior_mean_mat
zz <- pmin(zz, 1)
image(zz)

######################
case_control_variable <- "covariate_1"
case_control_idx <- which(colnames(covariates) == case_control_variable)
esvd_res <- eSVD_obj2$fit_First
covariates <- eSVD_obj2$covariates
library_size_variable <- "Log_UMI"
library_idx <- which(colnames(covariates) %in% c("Intercept", library_size_variable))

hist(covariates[,library_size_variable])
hist(esvd_res$z_mat[,library_size_variable])
hist(esvd_res$z_mat[,library_size_variable])

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(covariates[,-library_idx],
                       esvd_res$z_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_tmp <- tcrossprod(
  covariates[,library_idx], esvd_res$z_mat[,library_idx]
)
library_mat <- exp(library_tmp)
library_mat <- pmin(library_mat, 500)

plot(as.numeric(library_mat), rep(Matrix::rowSums(dat), times = p), asp = T)
plot(as.numeric(library_mat), rep(Matrix::rowSums(dat), times = p))

#########

library_tmp <- tcrossprod(
  covariates[,"Log_UMI"], esvd_res$z_mat[,"Log_UMI"]
)
library_mat <- exp(library_tmp)
library_mat <- pmin(library_mat, 500)
plot(as.numeric(library_mat), rep(Matrix::rowSums(dat), times = p))

