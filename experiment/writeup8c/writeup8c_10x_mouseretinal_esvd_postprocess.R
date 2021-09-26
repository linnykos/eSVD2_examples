rm(list=ls())
load("../../../../out/writeup8c/writeup8c_10x_mouseretinal_esvd.RData")
source("../writeup8b/plotting.R")

plotting_func(seurat_obj = retinal,
              mat,
              esvd_res,
              group.by = "Cluster",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8c/10x_mouseretinal_esvd_factor.png",
              heatmap_title = "Heatmap for Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8c/10x_mouseretinal_esvd_heatmap.png",
              scatter_title = "Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8c/10x_mouseretinal_esvd_scatter.png",
              width = 8)

plotting_func(seurat_obj = retinal,
              mat,
              esvd_res2,
              group.by = "Cluster",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8c/10x_mouseretinal_esvd_factor_2.png",
              heatmap_title = "Heatmap for Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8c/10x_mouseretinal_esvd_heatmap_2.png",
              scatter_title = "Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8c/10x_mouseretinal_esvd_scatter_2.png",
              width = 8)
