rm(list=ls())

load("../../../../out/writeup6/writeup6_10x_mousepancreas_glmpca_nb.RData")

glm_assay <- Seurat::CreateAssayObject(counts = pred_mat)
pancreas[["glm"]] <- glm_assay

Seurat::DefaultAssay(pancreas) <- "glm"
pancreas <- Seurat::NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas[["glm"]]@var.features <- rownames(pancreas)
pancreas <-  Seurat::ScaleData(pancreas)
pancreas <- Seurat::RunPCA(pancreas, features = Seurat::VariableFeatures(pancreas), verbose = F,
                           reduction.name = "glmpca")

set.seed(10)
pancreas <- Seurat::RunUMAP(pancreas, reduction = "glmpca", dims = 1:30, reduction.name = "glmumap")

plot1 <- Seurat::DimPlot(pancreas, reduction = "glmumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Pancreas (10x)\nGLM-PCA, Full, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/10x_mousepancreas_glmpca_full_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(pancreas@meta.data)

pancreas[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(pancreas, reduction = "glmfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Pancreas (10x)\nGLM-PCA, Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/10x_mousepancreas_glmpca_factor_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
