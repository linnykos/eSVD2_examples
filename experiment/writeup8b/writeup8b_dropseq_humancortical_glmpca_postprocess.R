rm(list=ls())

library(Seurat)

load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_glmpca_nb.RData")
mat2 <- as.matrix(Matrix::t(mat))
source("../writeup8b/plotting.R")

plotting_func_glmpca(seurat_obj = cortical,
                     mat = mat2,
                     glmpca_res,
                     group.by = "celltype",
                     max_value = quantile(mat2[mat2 > 2], probs = 0.999),
                     factor_title = "Human Cortical (Dropseq)\nGLM-PCA, Factor",
                     factor_filename = "../../../../out/fig/writeup8b/dropseq_humancortical_glmpca_factor.png",
                     heatmap_title = "Heatmap for Human Cortical (Dropseq)\nGLM-PCA",
                     heatmap_filename = "../../../../out/fig/writeup8b/dropseq_humancortical_glmpca_heatmap.png",
                     scatter_title = "Human Cortical (Dropseq)\nGLM-PCA:",
                     scatter_filename = "../../../../out/fig/writeup8b/dropseq_humancortical_glmpca_scatter.png")
