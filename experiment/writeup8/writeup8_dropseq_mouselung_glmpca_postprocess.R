rm(list=ls())
load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")
mat <- as.matrix(Matrix::t(mat))

png("../../../../out/fig/writeup8/dropseq_mouselung_glmpca_factor_heatmap.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(glmpca_res$factors,
                            membership_vec = as.factor(lung@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2, initialized glmGamPoi) Log-scale")
graphics.off()

max_value <- quantile(mat[mat > 2], probs = 0.999)
set.seed(10)
png("../../../../out/fig/writeup8/dropseq_mouselung_glmpca_scatter_onlymean.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scatterplot_poisson(mat, mean_mat = t(pred_mat),
                                 xlab = "Predicted mean",
                                 ylab = "Observed value",
                                 xlim = c(0, max_value),
                                 main = "Mouse Lung (Dropseq)\n(GLM-PCA, NB, only mean)")
graphics.off()

png("../../../../out/fig/writeup8/dropseq_mouselung_glmpca_scatter.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat,
                            mean_mat = t(pred_mat),
                            size_vec = rep(glmpca_res$glmpca_family$nb_theta[1], ncol(mat)),
                            log_scale = F,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Mouse Lung (Dropseq)\neSVD (GLM-PCA, NB):",
                            include_percentage_in_main = T)
graphics.off()
