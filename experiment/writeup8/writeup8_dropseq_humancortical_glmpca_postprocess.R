rm(list=ls())

library(Seurat)

load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_glmpca_nb.RData")
quantile(glmpca_res$dev) # this number is huge. what are these?

png("../../../../out/fig/writeup8/dropseq_humancortical_glmpca_factor_heatmap.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(glmpca_res$factors,
                            membership_vec = as.factor(cortical@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Human cortical (Dropseq)\nGLM-PCA (NB)")
graphics.off()

max_value <- quantile(mat@x[mat@x > 2], probs = 0.999)
set.seed(10)
png("../../../../out/fig/writeup8/dropseq_humancortical_glmpca_scatter_nb_onlymean.png",
    height = 2000, width = 2000, units = "px", res = 300)
plot_scatterplot_poisson(as.matrix(mat), mean_mat = pred_mat,
                         xlab = "Predicted mean",
                         ylab = "Observed value",
                         xlim = c(0, max_value),
                         main = "Human cortical (Dropseq)\nGLM-PCA (NB, but visualizing only mean)")
graphics.off()

set.seed(10)
png("../../../../out/fig/writeup8/dropseq_humancortical_glmpca_scatter_nb.png",
    height = 2000, width = 2000, units = "px", res = 300)
plot_scatterplot_nb(as.matrix(mat),
                    mean_mat = pred_mat,
                    size_vec = rep(glmpca_res$glmpca_family$nb_theta, ncol(pred_mat)),
                    log_scale = F,
                    xlab = "Predicted mean",
                    ylab = "Observed value",
                    xlim = c(0, max_value),
                    main = "Human cortical (Dropseq)\nGLM-PCA (NB)")
graphics.off()
