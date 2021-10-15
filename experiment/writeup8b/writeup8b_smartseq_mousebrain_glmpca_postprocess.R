rm(list=ls())

load("../../../../out/writeup6/writeup6_smartseq_mousebrain_glmpca_nb.RData")
mat2 <- as.matrix(Matrix::t(mat))
source("../writeup8b/plotting.R")

plotting_func_glmpca(seurat_obj = brain,
                     mat = mat2,
                     glmpca_res,
                     group.by = "subclass",
                     max_value = quantile(mat2[mat2 > 2], probs = 0.999),
                     factor_title = "Mouse Brain (Smartseq)\nGLM-PCA, Factor",
                     factor_filename = "../../../../out/fig/writeup8b/smartseq_mousebrain_glmpca_factor.png",
                     heatmap_title = "Heatmap for Mouse Brain (Smartseq)\nGLM-PCA",
                     heatmap_filename = "../../../../out/fig/writeup8b/smartseq_mousebrain_glmpca_heatmap.png",
                     scatter_title = "Mouse Brain (Smartseq)\nGLM-PCA:",
                     scatter_filename = "../../../../out/fig/writeup8b/smartseq_mousebrain_glmpca_scatter.png")

##########################

nat_mat <- tcrossprod(as.matrix(glmpca_res$factors), as.matrix(glmpca_res$loadings)) +
  tcrossprod(as.matrix(glmpca_res$X), as.matrix(glmpca_res$coefX))
nat_mat <- sweep(nat_mat, 1, glmpca_res$offsets, "+")
mean_mat <- exp(nat_mat)

max_value = quantile(mat2[mat2 > 2], probs = 0.999)
png("../../../../out/fig/writeup8b/smartseq_mousebrain_glmpca_scatter.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat2,
                            mean_mat = mean_mat,
                            size_vec = rep(glmpca_res$glmpca_family$nb_theta, ncol(mean_mat)),
                            quantile_shoulder = 0.5,
                            max_num = 50000,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Mouse Brain (Smartseq)\nGLM-PCA:",
                            include_percentage_in_main = T,
                            verbose = T)
graphics.off()
