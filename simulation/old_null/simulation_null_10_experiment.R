rm(list=ls())
load("../eSVD2_examples/simulation/simulation_null_10_esvd.RData")

dispersion_est_vec <- eSVD_obj$fit_Second$nuisance_vec
posterior_mat <- eSVD_obj$fit_Second$posterior_mean_mat
p <- ncol(posterior_mat)
denoised_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat,
                           eSVD_obj$fit_Second$y_mat) +
  tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
denoised_mat <- exp(denoised_mat)
sparsity_vec <- sapply(1:p, function(j){
  length(which(eSVD_obj$dat[,j] == 0))/nrow(eSVD_obj$dat)
})

cor_vec <- sapply(1:p, function(j){
  stats::cor(denoised_mat[,j], posterior_mat[,j])
})
quantile(cor_vec)

plot(jitter(1/dispersion_est_vec), cor_vec, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
plot(sparsity_vec, cor_vec, pch = 16, col = rgb(0.5,0.5,0.5,0.5))

