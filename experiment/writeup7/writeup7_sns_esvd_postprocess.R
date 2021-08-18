rm(list=ls())

library(Seurat)

load("../../../../out/writeup7/writeup7_sns_esvd.RData")

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
sns[["pred"]] <- glm_assay

Seurat::DefaultAssay(sns) <- "pred"
sns <- Seurat::NormalizeData(sns, normalization.method = "LogNormalize", scale.factor = 10000)
sns[["pred"]]@var.features <- rownames(sns)
sns <-  Seurat::ScaleData(sns)
sns <- Seurat::RunPCA(sns, features = Seurat::VariableFeatures(sns), verbose = F,
                           reduction.name = "esvdpca")

set.seed(10)
sns <- Seurat::RunUMAP(sns, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(sns, reduction = "esvdpca", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS)\neSVD, Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_esvd_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

######################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(sns@meta.data)

sns[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(sns, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS)\neSVD, Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_esvd_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
