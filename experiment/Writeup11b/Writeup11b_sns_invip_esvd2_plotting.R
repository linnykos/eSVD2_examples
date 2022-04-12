rm(list=ls())
load("../../../../out/Writeup11/Writeup11_sns_invip_esvd2.RData")


mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
tmp <- compute_posterior(mat = mat,
                         esvd_res = esvd_res_full,
                         nuisance_vec = nuisance_vec,
                         case_control_variable = "diagnosis_ASD",
                         alpha_max = 50,
                         nuisance_lower_quantile = 0.01)

set.seed(10)
png("../../../../out/fig/Writeup11b/Writeup11b_sns_invip_scatterplot_mean.png",
    height = 1800, width = 3500,
    units = "px", res = 300)
par(mfrow = c(1,2))
plot_scatterplot_mean(mat = mat,
                      esvd_res = esvd_res_full,
                      nuisance_vec = nuisance_vec,
                      case_control_variable = "diagnosis_ASD",
                      mean_type = "predicted",
                      main = "IN-VIP (Observed vs. predicted)",
                      xlab = "Predicted value",
                      ylab = "Observed value")

plot_scatterplot_mean(mat = mat,
                      esvd_res = esvd_res_full,
                      nuisance_vec = nuisance_vec,
                      case_control_variable = "diagnosis_ASD",
                      bool_logscale = T,
                      mean_type = "posterior",
                      main = "IN-VIP (Observed vs. posterior)",
                      xlab = "Posterior value (Log-scale)",
                      ylab = "Depth-normalized observed value (Log-scale)")
graphics.off()

###########

res <- plot_gene_variability(mat = mat,
                             library_size_vec = Matrix::rowSums(mat),
                             ngroups_cells = 5,
                             ngroups_genes = 6)
plot1 <- ggplot2::ggplot(res, ggplot2::aes(fill=cell_cluster, y=variance, x=gene_cluster)) +
  ggplot2::geom_bar(position="fill", stat="identity")
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/Writeup11b_sns_invip_geneVariability.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

other_covariates <- colnames(esvd_res_full$covariates)[colnames(esvd_res_full$covariates) != "Log_UMI"]
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,other_covariates,drop = F],
                       esvd_res_full$b_mat[,other_covariates,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
res <- plot_gene_variability(mat = mean_mat_nolib,
                             library_size_vec = Matrix::rowSums(mat),
                             ngroups_cells = 5,
                             ngroups_genes = 6)
plot1 <- ggplot2::ggplot(res, ggplot2::aes(fill=cell_cluster, y=variance, x=gene_cluster)) +
  ggplot2::geom_bar(position="fill", stat="identity")
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/Writeup11b_sns_invip_geneVariability_esvd.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

#########################

avg_expression <- Matrix::colMeans(mat)
set.seed(10)
res <- compute_gene_librarysize(mat = mat,
                                avg_expression = avg_expression,
                                library_size_vec = Matrix::rowSums(mat),
                                verbose = 1)
png("../../../../out/fig/Writeup11b/Writeup11b_sns_invip_gene_librarysize.png",
    height = 2000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(2,3))
plot_gene_librarysize(res,
                      xlab = "Library size",
                      ylab = "Observed gene count")
graphics.off()


other_covariates <- colnames(esvd_res_full$covariates)[colnames(esvd_res_full$covariates) != "Log_UMI"]
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,other_covariates,drop = F],
                       esvd_res_full$b_mat[,other_covariates,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
set.seed(10)
res2 <- compute_gene_librarysize(mat = mean_mat_nolib,
                                 avg_expression = avg_expression,
                                library_size_vec = Matrix::rowSums(mat),
                                verbose = 1)
png("../../../../out/fig/Writeup11b/Writeup11b_sns_invip_gene_librarysize_esvd.png",
    height = 2000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(2,3))
plot_gene_librarysize(res2,
                      xlab = "Library size",
                      ylab = "Denoised gene expression")
graphics.off()
