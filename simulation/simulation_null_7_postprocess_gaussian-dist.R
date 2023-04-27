rm(list=ls())
load("../eSVD2_examples/simulation/simulation_null_7_esvd.RData")

denoised_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat,
                           eSVD_obj$fit_Second$y_mat) +
  tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
denoised_mat <- exp(denoised_mat)

# first do it without any posterior correction
cell_individual_list <- lapply(1:max(covariate[,"Individual"]), function(i){which(covariate[,"Individual"] == i)})
mean_denoised <- t(sapply(cell_individual_list, function(idx){
  Matrix::colMeans(denoised_mat[idx,])
}))
var_denoised <- t(sapply(cell_individual_list, function(idx){
  matrixStats::colVars(denoised_mat[idx,])
}))
case_df_idx <- which(df[,"CC"] == 1)
control_df_idx <- which(df[,"CC"] == 0)
# j <- 1; plot(mean_denoised[,j], sqrt(var_denoised[,j]), col = c(rep(1,10),rep(2,10)), pch = 16, asp = T)

# https://djalil.chafai.net/blog/2010/04/30/wasserstein-distance-between-two-gaussians/
.gaussian_distances <- function(case_idx,
                                control_idx,
                                mean_vec,
                                var_vec){
  perm_mat1 <- cbind(rep(1:length(case_idx), each = length(control_idx)),
                     rep(1:length(control_idx), times = length(case_idx)))
  dist_vec1 <- sapply(1:nrow(perm_mat), function(j){
    i1 <- case_idx[perm_mat1[j,1]]; i2 <- control_idx[perm_mat1[j,2]]

    m1 <- mean_vec[i1]; m2 <- mean_vec[i2]
    v1 <- var_vec[i1]; v2 <- var_vec[i2]
    (m1-m2)^2 + (sqrt(v1)-sqrt(v2))^2
  })

  perm_mat2 <- utils::combn(1:length(case_idx),2)
  dist_vec2 <- sapply(1:nrow(perm_mat2), function(j){
    i1 <- case_idx[perm_mat2[j,1]]; i2 <- case_idx[perm_mat2[j,2]]

    m1 <- mean_vec[i1]; m2 <- mean_vec[i2]
    v1 <- var_vec[i1]; v2 <- var_vec[i2]
    (m1-m2)^2 + (sqrt(v1)-sqrt(v2))^2
  })

  perm_mat3 <- utils::combn(1:length(control_idx),2)
  dist_vec3 <- sapply(1:nrow(perm_mat3), function(j){
    i1 <- control_idx[perm_mat3[j,1]]; i2 <- control_idx[perm_mat3[j,2]]

    m1 <- mean_vec[i1]; m2 <- mean_vec[i2]
    v1 <- var_vec[i1]; v2 <- var_vec[i2]
    (m1-m2)^2 + (sqrt(v1)-sqrt(v2))^2
  })

  stats::median(dist_vec1)/c(stats::median(c(dist_vec2, dist_vec3)))
}

p <- ncol(mean_denoised)
dist_vec <- sapply(1:p, function(j){
  .gaussian_distances(case_idx = case_df_idx,
                      control_idx = control_df_idx,
                      mean_vec = mean_denoised[,j],
                      var_vec = var_denoised[,j])
})
quantile(dist_vec)
plot(log(dist_vec)); points(1:10, log(dist_vec[1:10]), pch = 16, col = 2)


#########################

posterior_mean_mat <- eSVD_obj$fit_Second$posterior_mean_mat
posterior_var_mat <- eSVD_obj$fit_Second$posterior_var_mat
avg_posterior_mean_mat <- t(sapply(cell_individual_list, function(idx){
  Matrix::colMeans(posterior_mean_mat[idx,])
}))
avg_posterior_var_mat <- t(sapply(cell_individual_list, function(idx){
  Matrix::colMeans(posterior_var_mat[idx,])
}))
case_df_idx <- which(df[,"CC"] == 1)
control_df_idx <- which(df[,"CC"] == 0)
case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_df_idx,])
control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_df_idx,])
case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[case_df_idx,],
  avg_posterior_var_mat = avg_posterior_var_mat[case_df_idx,]
)
control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[control_df_idx,],
  avg_posterior_var_mat = avg_posterior_var_mat[control_df_idx,]
)
# j <- 2; plot(avg_posterior_mean_mat[,j], sqrt(avg_posterior_var_mat[,j]), col = c(rep(1,10),rep(2,10)), pch = 16, asp = T)

dist_vec2 <- sapply(1:p, function(j){
  .gaussian_distances(case_idx = case_df_idx,
                      control_idx = control_df_idx,
                      mean_vec = avg_posterior_mean_mat[,j],
                      var_vec = avg_posterior_var_mat[,j])
})
quantile(dist_vec2)
plot(log(dist_vec2)); points(1:10, log(dist_vec2[1:10]), pch = 16, col = 2)
