rm(list=ls())

library(Seurat)
load("../../../../out/writeup8/writeup8_10x_mousepancreas_esvd_nb.RData")
quantile(nuisance_vec)

nat_mat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(covariates, esvd_res2$b_mat)
pred_mat <- eSVD2:::.mult_mat_vec(exp(nat_mat)/(1-exp(nat_mat)), nuisance_vec)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
pancreas[["pred"]] <- glm_assay

Seurat::DefaultAssay(pancreas) <- "pred"
pancreas <- Seurat::NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas[["pred"]]@var.features <- rownames(pancreas)
pancreas <-  Seurat::ScaleData(pancreas)
pancreas <- Seurat::RunPCA(pancreas, features = Seurat::VariableFeatures(pancreas), verbose = F,
                           reduction.name = "esvdpca")

set.seed(10)
pancreas <- Seurat::RunUMAP(pancreas, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(pancreas, reduction = "esvdpca", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Pancreas (10x)\neSVD (Alt), Full, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/10x_mousepancreas_esvd_full_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res2$x_mat))@cell.embeddings
rownames(tmp) <- rownames(pancreas@meta.data)

pancreas[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(pancreas, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Pancreas (10x)\neSVD (Alt), Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/10x_mousepancreas_esvd_factor_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
