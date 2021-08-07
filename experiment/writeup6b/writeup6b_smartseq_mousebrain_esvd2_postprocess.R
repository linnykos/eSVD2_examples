rm(list=ls())
load("../../../../out/writeup6b/writeup6b_smartseq_mousebrain_esvd2.RData")

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
brain[["pred"]] <- glm_assay

Seurat::DefaultAssay(brain) <- "pred"
brain <- Seurat::NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
brain[["pred"]]@var.features <- rownames(brain)
brain <-  Seurat::ScaleData(brain)
brain <- Seurat::RunPCA(brain, features = Seurat::VariableFeatures(brain), verbose = F,
                           reduction.name = "esvdpca")

set.seed(10)
brain <- Seurat::RunUMAP(brain, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(brain, reduction = "esvdpca", group.by = "subclass", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Brain (Smartseq)\neSVD (Alt), Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/smartseq_mousebrain_esvd2_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(brain@meta.data)

brain[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(brain, reduction = "esvdfactorumap", group.by = "subclass", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Brain (Smartseq)\neSVD (Alt), Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/smartseq_mousebrain_esvd2_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
