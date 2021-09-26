rm(list=ls())
load("../../../../out/writeup8c/writeup8c_10x_mouseretinal_esvd.RData")
source("../writeup8b/plotting.R")

plotting_func(seurat_obj = lung,
              mat,
              esvd_res,
              group.by = "celltype",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8c/dropseq_mouselung_esvd_factor.png",
              heatmap_title = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8c/dropseq_mouselung_esvd_heatmap.png",
              scatter_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8c/dropseq_mouselung_esvd_scatter.png")

plotting_func(seurat_obj = lung,
              mat,
              esvd_res2,
              group.by = "celltype",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8c/dropseq_mouselung_esvd_factor_2.png",
              heatmap_title = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8c/dropseq_mouselung_esvd_heatmap_2.png",
              scatter_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8c/dropseq_mouselung_esvd_scatter_2.png")
