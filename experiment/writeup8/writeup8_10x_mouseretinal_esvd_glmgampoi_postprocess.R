rm(list=ls())
load("../../../../out/writeup8/writeup8_10x_mouseretinal_esvd_nb_glmgampoi.RData")

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(retinal@meta.data)

retinal[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

options(ggrepel.max.overlaps = Inf)
plot1 <- Seurat::DimPlot(retinal, reduction = "esvdfactorumap", group.by = "Cluster", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse retinal (10x)\neSVD (via GLMGamPoi), Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/10x_mouseretinal_esvd_factor_nb_glmgampoi_umap.png",
                plot1, device = "png", width = 8, height = 5, units = "in")
