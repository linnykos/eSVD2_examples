rm(list=ls())
load("../eSVD2_examples/simulation/simulation_null_6_esvd.RData")

gene_plot_idx <- c(which(y_block_assignment == 1), which(y_block_assignment == 2), which(y_block_assignment == 3))

plot(eSVD_obj$fit_Second$z_mat[,"CC"])

# let's see what the heatmap looks like without posterior
denoised_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat,
                           eSVD_obj$fit_Second$y_mat) +
  tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
denoised_mat <- exp(denoised_mat)
cor_mat <- cor(denoised_mat)
image(cor_mat[gene_plot_idx, gene_plot_idx])

# the a few naive methods:
# 1) A Wilcoxon test
p <- ncol(denoised_mat)
case_idx <- which(seurat_obj$CC == 1)
control_idx <- which(seurat_obj$CC == 0)
wilcoxon_pval <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  set.seed(j)
  res <- stats::wilcox.test(x = denoised_mat[case_idx,j],
                            y = denoised_mat[control_idx,j])
  res$p.value
})
plot(sort(wilcoxon_pval[-c(1:10)]),
     seq(0,1,length.out = length(wilcoxon_pval[-c(1:10)])), asp = T,
     main = "Wilcoxon")
lines(c(0,1), c(0,1), col = 2, lty = 2)

# 2) Collapse onto pseudo-bulk, and then test
cell_individual_list <- lapply(1:max(covariate[,"Individual"]), function(i){which(covariate[,"Individual"] == i)})
mean_denoised <- t(sapply(cell_individual_list, function(idx){
  Matrix::colMeans(denoised_mat[idx,])
}))
case_df_idx <- which(df[,"CC"] == 1)
control_df_idx <- which(df[,"CC"] == 0)
wilcoxon_df_pval <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  set.seed(j)
  res <- stats::wilcox.test(x = mean_denoised[case_df_idx,j],
                            y = mean_denoised[control_df_idx,j])
  res$p.value
})
plot(sort(wilcoxon_df_pval[-c(1:10)]),
     seq(0,1,length.out = length(wilcoxon_df_pval[-c(1:10)])), asp = T,
     main = "Wilcoxon Pseudobulk")
lines(c(0,1), c(0,1), col = 2, lty = 2)

plot(wilcoxon_pval, wilcoxon_df_pval)

# 3) A Gaussian individual test, but no posterior
cell_individual_list <- lapply(1:max(covariate[,"Individual"]), function(i){which(covariate[,"Individual"] == i)})
mean_denoised <- t(sapply(cell_individual_list, function(idx){
  Matrix::colMeans(denoised_mat[idx,])
}))
var_denoised <- t(sapply(cell_individual_list, function(idx){
  matrixStats::colVars(denoised_mat[idx,])
}))
# j <- 15; plot(mean_denoised[,j], sqrt(var_denoised[,j]), col = c(rep(1,10),rep(2,10)), pch = 16, asp = T)
case_df_idx <- which(df[,"CC"] == 1)
control_df_idx <- which(df[,"CC"] == 0)
n1 <- length(case_df_idx)
n2 <- length(control_df_idx)
case_gaussian_mean <- Matrix::colMeans(mean_denoised[case_df_idx,,drop = F])
control_gaussian_mean <- Matrix::colMeans(mean_denoised[control_df_idx,,drop = F])
case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = mean_denoised[case_df_idx,,drop = F],
  avg_posterior_var_mat = var_denoised[case_df_idx,,drop = F]
)
control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = mean_denoised[control_df_idx,,drop = F],
  avg_posterior_var_mat = var_denoised[control_df_idx,,drop = F]
)
teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
  (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))
# plot(teststat_vec)
# oh... it didn't seem that important
