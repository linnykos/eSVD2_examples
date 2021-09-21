rm(list=ls())
load("../../../../out/writeup8b/writeup8b_dropseq_mouselung_esvd2.RData")
library("plotting.R")

plotting_func(seurat_obj = lung,
              esvd_res,
              covariates = covariates,
              group.by = "celltype",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_esvd2_factor.png",
              heatmap_title = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_esvd2_heatmap.png",
              scatter_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_esvd2_scatter.png")


plotting_func(seurat_obj = lung,
              esvd_res2,
              covariates = covariates,
              group.by = "celltype",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_esvd2_factor_2.png",
              heatmap_title = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_esvd2_heatmap_2.png",
              scatter_title = "Mouse Lung (Dropseq)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8b/dropseq_mouselung_esvd2_scatter_2.png")
