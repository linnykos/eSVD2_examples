rm(list=ls())
load("../../../../out/writeup8c/writeup8c_smartseq_mousebrain_esvd.RData")
source("../writeup8b/plotting.R")

plotting_func(seurat_obj = brain,
              mat,
              esvd_res,
              group.by = "subclass",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Brain (Smart-seq)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8c/smartseq_mousebrain_esvd_factor.png",
              heatmap_title = "Heatmap for Mouse Brain (Smart-seq)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8c/smartseq_mousebrain_esvd_heatmap.png",
              scatter_title = "Mouse Brain (Smart-seq)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8c/smartseq_mousebrain_esvd_scatter.png",
              width = 6)

plotting_func(seurat_obj = brain,
              mat,
              esvd_res2,
              group.by = "subclass",
              max_value = quantile(mat[mat > 2], probs = 0.999),
              factor_title = "Mouse Brain (Smart-seq)\neSVD (NB2, reparam+gene-specific theta), Factor",
              factor_filename = "../../../../out/fig/writeup8c/smartseq_mousebrain_esvd_factor_2.png",
              heatmap_title = "Heatmap for Mouse Brain (Smart-seq)\neSVD (NB2, reparam+gene-specific theta)",
              heatmap_filename = "../../../../out/fig/writeup8c/smartseq_mousebrain_esvd_heatmap_2.png",
              scatter_title = "Mouse Brain (Smart-seq)\neSVD (NB2, reparam+gene-specific theta):",
              scatter_filename = "../../../../out/fig/writeup8c/smartseq_mousebrain_esvd_scatter_2.png",
              width = 6)
