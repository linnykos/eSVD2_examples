rm(list=ls())
load("../../../../out/writeup8/writeup8_citeseq_bm_esvd_nb_glmgampoi.RData")

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

options(ggrepel.max.overlaps = Inf)
plot1 <- Seurat::DimPlot(bm, reduction = "esvdfactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\neSVD (via GLMGamPoi), Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/10x_mousebm_esvd_factor_nb_glmgampoi_umap.png",
                plot1, device = "png", width = 7, height = 5, units = "in")
