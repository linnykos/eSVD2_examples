rm(list=ls())
load("../../../../out/writeup8/writeup8_smartseq_mousebrain_esvd_nb_glmgampoi.RData")

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(brain@meta.data)

brain[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(brain, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse brain (Smartseq)\neSVD (via GLMGamPoi), Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/smartseq_mousebrain_esvd_factor_nb_glmgampoi_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")
