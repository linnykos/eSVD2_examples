rm(list=ls())
load("../../../../out/writeup6b/writeup6b_citeseq_bm_esvd2.RData")

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
bm[["pred"]] <- glm_assay

Seurat::DefaultAssay(bm) <- "pred"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm[["pred"]]@var.features <- rownames(bm)
bm <-  Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, features = Seurat::VariableFeatures(bm), verbose = F,
                           reduction.name = "esvdpca")

set.seed(10)
bm <- Seurat::RunUMAP(bm, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(bm, reduction = "esvdpca", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\neSVD (Alt), Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/citeseq_bm_esvd2_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "esvdfactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\neSVD (Alt), Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/citeseq_bm_esvd2_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
