rm(list=ls())

library(Seurat)

load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_esvd2.RData")

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
cortical[["pred"]] <- glm_assay

Seurat::DefaultAssay(cortical) <- "pred"
cortical <- Seurat::NormalizeData(cortical, normalization.method = "LogNormalize", scale.factor = 10000)
cortical[["pred"]]@var.features <- rownames(cortical)
cortical <-  Seurat::ScaleData(cortical)
cortical <- Seurat::RunPCA(cortical, features = Seurat::VariableFeatures(cortical), verbose = F,
                           reduction.name = "esvdpca")

set.seed(10)
cortical <- Seurat::RunUMAP(cortical, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(cortical, reduction = "esvdpca", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\neSVD (Alt), Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_esvd2_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

######################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(cortical@meta.data)

cortical[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(cortical, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\neSVD (Alt), Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_esvd2_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
