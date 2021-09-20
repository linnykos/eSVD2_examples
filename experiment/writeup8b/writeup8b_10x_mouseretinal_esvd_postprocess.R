rm(list=ls())

library(Seurat)
load("../../../../out/writeup8b/writeup8b_10x_mouseretinal_esvd.RData")
source("plotting.R")

quantile(esvd_res$nuisance_param_vec)
quantile(esvd_res$b_mat[,1])
quantile(esvd_res$b_mat[,2])
quantile(esvd_res3$nuisance_param_vec)
quantile(esvd_res3$b_mat[,1])

plotting_func(seurat_obj = retinal,
              esvd_res,
              covariates = covariates,
              group.by = "Cluster",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_factor1.png",
              heatmap_title = "Heatmap for Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_heatmap1.png",
              scatter_title = "Mouse Retinal (10x)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_scatter1.png")


plotting_func(seurat_obj = retinal,
              esvd_res2,
              covariates = covariates,
              group.by = "Cluster",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Retinal (10x)\neSVD (NB2, gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_factor2.png",
              heatmap_title = "Heatmap for Mouse Retinal (10x)\neSVD (NB2, gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_heatmap2.png",
              scatter_title = "Mouse Retinal (10x)\neSVD (NB2, gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_scatter2.png")

plotting_func(seurat_obj = retinal,
              esvd_res3,
              covariates = covariates[,1,drop = F],
              group.by = "Cluster",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Retinal (10x)\neSVD (NB2, GLM-PCA-like), Factor",
              factor_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_factor3.png",
              heatmap_title = "Heatmap for Mouse Retinal (10x)\neSVD (NB2, GLM-PCA-like)",
              heatmap_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_heatmap3.png",
              scatter_title = "Mouse Retinal (10x)\neSVD (NB2, GLM-PCA-like):",
              scatter_filename = "../../../../out/fig/writeup8b/10x_mouseretinal_esvd_scatter3.png")

