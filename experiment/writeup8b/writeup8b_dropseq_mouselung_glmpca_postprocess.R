rm(list=ls())

load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")
mat2 <- as.matrix(Matrix::t(mat))
source("../writeup8b/plotting.R")

plotting_func_glmpca(seurat_obj = lung,
                     mat = mat2,
                     glmpca_res,
                     group.by = "celltype",
                     max_value = quantile(mat2[mat2 > 2], probs = 0.999),
                     factor_title = "Mouse Lung (Dropseq)\nGLM-PCA, Factor",
                     factor_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_glmpca_factor.png",
                     heatmap_title = "Heatmap for Mouse Lung (Dropseq)\nGLM-PCA",
                     heatmap_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_glmpca_heatmap.png",
                     scatter_title = "Mouse Lung (Dropseq)\nGLM-PCA:",
                     scatter_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_glmpca_scatter.png")
