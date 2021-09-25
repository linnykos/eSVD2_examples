rm(list=ls())

load("../../../../out/writeup6/writeup6_10x_mouseretinal_glmpca_nb.RData")
mat2 <- as.matrix(Matrix::t(mat))

plotting_func_glmpca(seurat_obj = retinal,
                     mat = mat2,
                     glmpca_res,
                     group.by = "Cluster",
                     max_value = quantile(mat2[mat2 > 2], probs = 0.999),
                     factor_title = "Mouse Retinal (10x)\nGLM-PCA, Factor",
                     factor_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_glmpca_factor.png",
                     heatmap_title = "Heatmap for Mouse Retinal (10x)\nGLM-PCA",
                     heatmap_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_glmpca_heatmap.png",
                     scatter_title = "Mouse Retinal (10x)\nGLM-PCA:",
                     scatter_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_glmpca_scatter.png")
