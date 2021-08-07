rm(list=ls())
load("../../../../out/writeup6b/writeup6b_dropseq_mouselung_esvd2.RData")

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
lung[["pred"]] <- glm_assay

Seurat::DefaultAssay(lung) <- "pred"
lung <- Seurat::NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung[["pred"]]@var.features <- rownames(lung)
lung <-  Seurat::ScaleData(lung)
lung <- Seurat::RunPCA(lung, features = Seurat::VariableFeatures(lung), verbose = F,
                           reduction.name = "esvdpca")

set.seed(10)
lung <- Seurat::RunUMAP(lung, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdpca", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (Alt), Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_mouselung_esvd2_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (Alt), Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_mouselung_esvd2_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
