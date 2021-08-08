rm(list=ls())
load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_fasttopics.RData")

.mult_vec_mat <- function(vec, mat){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == nrow(mat))
  vec * mat
}

pred_mat_raw <- topic_res$L %*% t(topic_res$F)
pred_mat <- .mult_vec_mat(rowSums(as.matrix(mat)), topic_res$L %*% t(topic_res$F))

ft_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
cortical[["ft"]] <- ft_assay

Seurat::DefaultAssay(cortical) <- "ft"
cortical <- Seurat::NormalizeData(cortical, normalization.method = "LogNormalize", scale.factor = 10000)
cortical[["ft"]]@var.features <- rownames(cortical)
cortical <-  Seurat::ScaleData(cortical)
cortical <- Seurat::RunPCA(cortical, features = Seurat::VariableFeatures(cortical), verbose = F,
                     reduction.name = "ftpca")

set.seed(10)
cortical <- Seurat::RunUMAP(cortical, reduction = "ftpca", dims = 1:30, reduction.name = "ftumap")

plot1 <- Seurat::DimPlot(cortical, reduction = "ftumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\nFastTopics, Full")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_fasttopics_full_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(topic_res$L))@cell.embeddings
rownames(tmp) <- rownames(cortical@meta.data)

cortical[["ftfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(cortical, reduction = "ftfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\nFastTopics, Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_fasttopics_factor_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
