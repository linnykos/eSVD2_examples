rm(list=ls())

load("../../../../out/writeup6/writeup6_smartseq_mousebrain_glmpca_nb.RData")

glm_assay <- Seurat::CreateAssayObject(counts = pred_mat)
brain[["glm"]] <- glm_assay

Seurat::DefaultAssay(brain) <- "glm"
brain <- Seurat::NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
brain[["glm"]]@var.features <- rownames(brain)
brain <-  Seurat::ScaleData(brain)
brain <- Seurat::RunPCA(brain, features = Seurat::VariableFeatures(brain), verbose = F,
                       reduction.name = "glmpca")

set.seed(10)
brain <- Seurat::RunUMAP(brain, reduction = "glmpca", dims = 1:30, reduction.name = "glmumap")

plot1 <- Seurat::DimPlot(brain, reduction = "glmumap", group.by = "subclass", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Brain (Smartseq)\nGLM-PCA, Full, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/smartseq_mousebrain_glmpca_full_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(brain@meta.data)

brain[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(brain, reduction = "glmfactorumap", group.by = "subclass", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Brain (Smartseq)\nGLM-PCA, Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/smartseq_mousebrain_glmpca_factor_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
