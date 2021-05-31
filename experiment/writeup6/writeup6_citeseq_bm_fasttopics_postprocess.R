rm(list=ls())
load("../../../../out/writeup6/writeup6_citeseq_bm_fasttopics.RData")

.mult_vec_mat <- function(vec, mat){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == nrow(mat))
  vec * mat
}

pred_mat_raw <- topic_res$L %*% t(topic_res$F)
pred_mat <- .mult_vec_mat(rowSums(mat), topic_res$L %*% t(topic_res$F))

ft_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
bm[["ft"]] <- ft_assay

Seurat::DefaultAssay(bm) <- "ft"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm[["ft"]]@var.features <- rownames(bm)
bm <-  Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, features = Seurat::VariableFeatures(bm), verbose = F,
                     reduction.name = "ftpca")

set.seed(10)
bm <- Seurat::RunUMAP(bm, reduction = "ftpca", dims = 1:30, reduction.name = "ftumap")

plot1 <- Seurat::DimPlot(bm, reduction = "ftumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nFastTopics, Full")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_fasttopics_full_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(topic_res$L))@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["ftfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "ftfactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nFastTopics, Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_fasttopics_factor_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

####################
#####################

membership_vec <- as.factor(bm@meta.data$celltype.l2)
celltype_vec <- names(which(table(membership_vec) > nrow(mat)/50))
pval_mat <- pvalue_pairwise(pred_mat_raw, membership_vec, celltype_vec)

