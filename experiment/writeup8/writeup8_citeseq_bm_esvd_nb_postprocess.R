rm(list=ls())

library(Seurat)
load("../../../../out/writeup8/writeup8_citeseq_bm_esvd_nb.RData")
quantile(nuisance_vec)

nat_mat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(covariates, esvd_res2$b_mat)
pred_mat <- eSVD2:::.mult_mat_vec(exp(nat_mat)/(1-exp(nat_mat)), nuisance_vec)
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
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\neSVD, Full, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/citeseq_bm_esvd_full_nb_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res2$x_mat))@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "esvdfactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\neSVD, Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/citeseq_bm_esvd_factor_nb_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")
