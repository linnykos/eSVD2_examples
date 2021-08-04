rm(list=ls())

library(Seurat)

load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_glmpca_nb.RData")
load("../../../../data/dropseq_humancortical/dropseq_humancortical_formatted.RData")

glm_assay <- Seurat::CreateAssayObject(counts = pred_mat)
cortical[["glm"]] <- glm_assay

Seurat::DefaultAssay(cortical) <- "glm"
cortical <- Seurat::NormalizeData(cortical, normalization.method = "LogNormalize", scale.factor = 10000)
cortical[["glm"]]@var.features <- rownames(cortical)
cortical <-  Seurat::ScaleData(cortical)
cortical <- Seurat::RunPCA(cortical, features = Seurat::VariableFeatures(cortical), verbose = F,
                          reduction.name = "glmpca")

set.seed(10)
cortical <- Seurat::RunUMAP(cortical, reduction = "glmpca", dims = 1:30, reduction.name = "glmumap")

plot1 <- Seurat::DimPlot(cortical, reduction = "glmumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\nGLM-PCA, Full, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_glmpca_full_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(retinal@meta.data)

cortical[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(cortical, reduction = "glmfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\nGLM-PCA, Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_glmpca_factor_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

