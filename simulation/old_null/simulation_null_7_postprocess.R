rm(list=ls())
load("../eSVD2_examples/simulation/simulation_null_7_esvd.RData")

plot(eSVD_obj$fit_Second$z_mat[,"CC"])
plot(eSVD_obj$teststat_vec)

plot(sort(multtest_res$pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(multtest_res$pvalue_vec[-c(1:10)])), asp = T)
lines(c(0,1), c(0,1), col = 2, lty = 2)

# plot the "naive" p-values
df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
teststat_vec <- eSVD_obj$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})
null_mean <- 0; null_sd <- 1
pvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    2*Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = F)
  } else {
    2*Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = F)
  }
})
plot(sort(pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(pvalue_vec[-c(1:10)])), asp = T,
     main = "Mixture-Gaussian Pseudobulk")
lines(c(0,1), c(0,1), col = 2, lty = 2)


# let's see what the heatmap looks like without posterior
denoised_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat,
                           eSVD_obj$fit_Second$y_mat) +
  tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
denoised_mat <- exp(denoised_mat)
cor_mat <- cor(denoised_mat)
image(cor_mat[gene_plot_idx, gene_plot_idx])

##################################################################
##################################################################

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

##################################################################
##################################################################

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

gaussian_teststat <- sapply(1:length(wilcoxon_df_pval), function(j){
  if(mean(mean_denoised[case_df_idx,j]) > mean(mean_denoised[control_df_idx,j])){
    stats::qnorm(wilcoxon_df_pval[j]/2)
  } else {
    stats::qnorm(1-wilcoxon_df_pval[j]/2)
  }
})
res2 <- eSVD2:::multtest(gaussian_teststat)

plot(sort(res2$pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(res2$pvalue_vec[-c(1:10)])), asp = T,
     main = "Wilcoxon Pseudobulk, Multtest")
lines(c(0,1), c(0,1), col = 2, lty = 2)
# uhoh. nonsense


##################################################################
##################################################################

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
plot(teststat_vec)

numerator_vec <- (case_gaussian_var/n1 + control_gaussian_var/n2)^2
denominator_vec <- (case_gaussian_var/n1)^2/(n1-1) + (control_gaussian_var/n2)^2/(n2-1)
df_vec <- numerator_vec/denominator_vec
p <- ncol(mean_denoised)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})
null_mean <- 0; null_sd <- 1
pvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    2*Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = F)
  } else {
    2*Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = F)
  }
})
res2 <- eSVD2:::multtest(gaussian_teststat)

plot(sort(pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(pvalue_vec[-c(1:10)])), asp = T,
     main = "Mixture-Gaussian Pseudobulk")
lines(c(0,1), c(0,1), col = 2, lty = 2)

plot(sort(res2$pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(res2$pvalue_vec[-c(1:10)])), asp = T,
     main = "Mixture-Gaussian Pseudobulk, Multtest")
lines(c(0,1), c(0,1), col = 2, lty = 2)
