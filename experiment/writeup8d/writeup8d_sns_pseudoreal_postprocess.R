rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

stats::cor(esvd_res2$nuisance_param_vec, true_esvd$nuisance_param_vec)
quantile(esvd_res2$nuisance_param_vec)
for(j in 1:ncol(esvd_res2$covariates)){
  print(paste0(colnames(esvd_res2$covariates)[j], ": ",j))
  print(stats::cor(esvd_res2$b_mat[,j], true_esvd$b_mat[,j]))
}

png("../../../../out/fig/writeup8d/sns_pseudoreal_nuisance_scatter.png",
    height = 1200, width = 1200, units = "px", res = 300)
plot(true_esvd$nuisance_param_vec,
     esvd_res2$nuisance_param_vec,
     xlab = "True nuisance", ylab = "Estimated nuisance",
     asp = T, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
graphics.off()

########
quantile(esvd_res2$b_mat[autism_gene_idx,"diagnosis_ASD"])
quantile(esvd_res2$b_mat[-autism_gene_idx,"diagnosis_ASD"])

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), p)
col_vec[autism_gene_idx] <- 2
cex_vec <- rep(0.5, p)
cex_vec[autism_gene_idx] <- 1
png("../../../../out/fig/writeup8d/sns_pseudoreal_asdcoef_scatter.png",
    height = 1200, width = 1200, units = "px", res = 300)
plot(true_esvd$b_mat[,"diagnosis_ASD"],
     esvd_res2$b_mat[,"diagnosis_ASD"],
     xlab = "True ASD coefficient", ylab = "Estimated ASD coefficient",
     asp = T, pch = 16, col = col_vec, cex = cex_vec)
graphics.off()

png("../../../../out/fig/writeup8d/sns_pseudoreal_hist_diagnosis.png",
    width = 1500, height = 1500, units = "px", res = 300)
hist(esvd_res2$b_mat[,"diagnosis_ASD"], main = "Pseudoreal,\nCoefficient for diagnosis (ASD)",
     col = "gray", xlab = "ASD coefficient", breaks = 50,
     xlim = range(esvd_res2$b_mat[,"diagnosis_ASD"]))
rug(esvd_res2$b_mat[autism_gene_idx,"diagnosis_ASD"], col = "red")
graphics.off()

###############################

est_nat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(esvd_res2$covariates, esvd_res2$b_mat)
est_mean <- exp(est_nat)

max_value <- 30
png("../../../../out/fig/writeup8d/sns_pseudoreal_est_scatterplot.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat,
                            mean_mat = est_mean,
                            size_vec = esvd_res2$nuisance_param_vec,
                            quantile_shoulder = 0.5,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Pseudoreal: Estimated",
                            included_col = rgb(0.5, 0.5, 0.5, 0.5),
                            excluded_col = rgb(0.5, 0, 0, 0.5),
                            include_percentage_in_main = T,
                            verbose = T)
graphics.off()

xlim <- c(0, 30)
png("../../../../out/fig/writeup8d/sns_pseudoreal_est_vs_true_mean.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
idx <- which(mat > 0)
idx <- sample(idx, 1e5)
plot(NA, xlim = xlim, ylim = xlim, xlab = "True lambda",
     ylab = "Estimated mean", asp = T)
points(x = eSVD2:::.mult_vec_mat(s_vec, lambda_mat)[idx], y = est_mean[idx], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.5))
graphics.off()

#####################
#####################
#####################
#####################
#####################

# what about the version where all genes share the same nuisance?

quantile(esvd_res$b_mat[,"Log-UMI"])
for(j in 1:ncol(esvd_res$covariates)){
  print(paste0(colnames(esvd_res$covariates)[j], ": ",j))
  print(stats::cor(esvd_res$b_mat[,j], true_esvd$b_mat[,j]))
}


col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), p)
col_vec[autism_gene_idx] <- 2
cex_vec <- rep(0.5, p)
cex_vec[autism_gene_idx] <- 1
png("../../../../out/fig/writeup8d/sns_pseudoreal_asdcoef_scatter_samenuisance.png",
    height = 1200, width = 1200, units = "px", res = 300)
plot(true_esvd$b_mat[,"diagnosis_ASD"],
     esvd_res$b_mat[,"diagnosis_ASD"],
     xlab = "True ASD coefficient", ylab = "Estimated ASD coefficient",
     asp = T, pch = 16, col = col_vec, cex = cex_vec)
graphics.off()

png("../../../../out/fig/writeup8d/sns_pseudoreal_hist_diagnosis_samenuisance.png",
    width = 1500, height = 1500, units = "px", res = 300)
hist(esvd_res$b_mat[,"diagnosis_ASD"], main = "Pseudoreal,\nCoefficient for diagnosis (ASD)",
     col = "gray", xlab = "ASD coefficient", breaks = 50,
     xlim = range(esvd_res$b_mat[,"diagnosis_ASD"]))
rug(esvd_res$b_mat[autism_gene_idx,"diagnosis_ASD"], col = "red")
graphics.off()

###############################

est_nat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(esvd_res$covariates, esvd_res$b_mat)
est_mean <- exp(est_nat)
true_nat <- tcrossprod(true_esvd$x_mat, true_esvd$y_mat) + tcrossprod(true_esvd$covariates, true_esvd$b_mat)
true_mean <- exp(true_nat)
true_mean[true_mean > 100] <- 100
mean((true_mean - mean(true_mean))*as.numeric(est_mean - mean(est_mean)))/(sd(true_mean) * sd(est_mean))

max_value <- 30
png("../../../../out/fig/writeup8d/sns_pseudoreal_est_scatterplot_samenuisance.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat,
                            mean_mat = est_mean,
                            size_vec = esvd_res$nuisance_param_vec,
                            quantile_shoulder = 0.5,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Pseudoreal: Estimated",
                            included_col = rgb(0.5, 0.5, 0.5, 0.5),
                            excluded_col = rgb(0.5, 0, 0, 0.5),
                            include_percentage_in_main = T,
                            verbose = T)
graphics.off()

xlim <- c(0, 30)
png("../../../../out/fig/writeup8d/sns_pseudoreal_est_vs_true_mean_samenuisance.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
idx <- which(mat > 0)
idx <- sample(idx, 1e5)
plot(NA, xlim = xlim, ylim = xlim, xlab = "True lambda and library size",
     ylab = "Estimated mean", asp = T)
points(x = eSVD2:::.mult_vec_mat(s_vec, lambda_mat)[idx], y = est_mean[idx], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.5))
graphics.off()

######################################



