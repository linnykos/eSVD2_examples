rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../out/simulation/simulation_1.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

true_coefficient_mat <- t(apply(nat_mat, 2, function(y){
  df <- data.frame(cbind(y, covariates))
  colnames(df)[1] <- "y"
  lm_res <- stats::lm(y ~ . - 1, data = df)
  stats::coef(lm_res)
}))

load("../../out/simulation/simulation_1_esvd.RData")
r <- ncol(eSVD_obj$fit_Second$z_mat)
quantile_list <- vector("list", length(seq(1, 0.6, by = -.05)))
quantile_list[[1]] <- sapply(1:r, function(j){
  stats::cor(eSVD_obj$fit_Second$z_mat[,j], true_coefficient_mat[,j])
})

# true_residual_mat <- apply(nat_mat, 2, function(y){
#   df <- data.frame(cbind(y, covariates))
#   colnames(df)[1] <- "y"
#   lm_res <- stats::lm(y ~ . - 1, data = df)
#   stats::residuals(lm_res) + stats::coef(lm_res)["cc"]*covariates[,"cc"]
# })
# nat_mat_1 <- tcrossprod(eSVD_obj[["fit_Second"]]$x_mat, eSVD_obj[["fit_Second"]]$y_mat)
# nat_mat_2 <- tcrossprod(eSVD_obj[["covariates"]][,"cc_1"], eSVD_obj[["fit_Second"]]$z_mat[,"cc_1"])
# nat_mat <- nat_mat_1 + nat_mat_2
# n <- nrow(nat_mat)
# quantile_list[[1]] <- sapply(1:n, function(i){
#   stats::cor(true_residual_mat[i,], nat_mat[i,])
# })

downsample_values <- seq(0.95, 0.6, by = -0.05)
for(k in 1:length(downsample_values)){
  downsample_value <- downsample_values[k]
  print(paste0("Working on ", downsample_value))
  load(paste0("../../out/simulation/simulation_1_esvd_downsampled-", downsample_value, ".RData"))

  quantile_list[[k+1]] <- sapply(1:r, function(j){
    stats::cor(eSVD_obj$fit_Second$z_mat[,j], true_coefficient_mat[,j])
  })

  # nat_mat_1 <- tcrossprod(eSVD_obj[["fit_Second"]]$x_mat, eSVD_obj[["fit_Second"]]$y_mat)
  # nat_mat_2 <- tcrossprod(eSVD_obj[["covariates"]][,"cc_1"], eSVD_obj[["fit_Second"]]$z_mat[,"cc_1"])
  # nat_mat <- nat_mat_1 + nat_mat_2
  #
  # quantile_list[[k+1]] <- sapply(1:n, function(i){
  #   stats::cor(true_residual_mat[i,], nat_mat[i,])
  # })
}

################

load("../../out/simulation/simulation_1_nbreg_downsampled.RData")
nb_quantile_list <- vector("list", length(seq(1, 0.6, by = -.05)))
for(k in 1:length(downsample_coef_list)){
  downsample_value <- downsample_values[k]
  nb_quantile_list[[k]] <- sapply(1:r, function(j){
    stats::cor(downsample_coef_list[[k]][,j], true_coefficient_mat[,j])
  })
}
