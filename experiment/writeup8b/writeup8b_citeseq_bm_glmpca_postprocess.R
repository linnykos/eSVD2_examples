rm(list=ls())

load("../../../../out/writeup6/writeup6_citeseq_bm_glmpca_nb.RData")
mat2 <- as.matrix(Matrix::t(mat))
source("../writeup8b/plotting.R")

plotting_func_glmpca(seurat_obj = bm,
                     mat = mat2,
                     glmpca_res,
                     group.by = "celltype.l2",
                     max_value = quantile(mat2[mat2 > 2], probs = 0.999),
                     factor_title = "Bone marrow (CITE-seq)\nGLM-PCA, Factor",
                     factor_filename = "../../../../out/fig/writeup8b/citeseq_bm_glmpca_factor.png",
                     heatmap_title = "Heatmap for Bone marrow (CITE-seq)\nGLM-PCA",
                     heatmap_filename = "../../../../out/fig/writeup8b/citeseq_bm_glmpca_heatmap.png",
                     scatter_title = "Bone marrow (CITE-seq)\nGLM-PCA:",
                     scatter_filename = "../../../../out/fig/writeup8b/citeseq_bm_glmpca_scatter.png")
