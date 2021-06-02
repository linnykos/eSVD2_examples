rm(list=ls())

load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")

glm_assay <- Seurat::CreateAssayObject(counts = pred_mat)
lung[["glm"]] <- glm_assay

Seurat::DefaultAssay(lung) <- "glm"
lung <- Seurat::NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung[["glm"]]@var.features <- rownames(lung)
lung <-  Seurat::ScaleData(lung)
lung <- Seurat::RunPCA(lung, features = Seurat::VariableFeatures(lung), verbose = F,
                          reduction.name = "glmpca")

set.seed(10)
lung <- Seurat::RunUMAP(lung, reduction = "glmpca", dims = 1:30, reduction.name = "glmumap")

plot1 <- Seurat::DimPlot(lung, reduction = "glmumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\nGLM-PCA, Full, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/dropseq_mouselung_glmpca_full_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "glmfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\nGLM-PCA, Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/dropseq_mouselung_glmpca_factor_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

####################################

source("plotting_func.R")

celltype <- as.factor(lung@meta.data["celltype"][,1])
png("../../../../out/fig/writeup6/dropseq_mouselung_glmpca_nb_score.png", height = 1500,
    width = 3000, res = 300, units = "px")
par(mfrow = c(1,2))
plot_scores_heatmap(glmpca_res$factors, membership_vec = celltype,
                    bool_center = T, bool_scale = F,
                    bool_log = T, scaling_power = 2,
                    main = "Scores (GLM-PCA, NB, Unscaled)")
plot_scores_heatmap(glmpca_res$factors, membership_vec = celltype,
                    bool_center = T, bool_scale = T,
                    bool_log = T, scaling_power = 2,
                    main = "Scores (GLM-PCA, NB, Scaled)")
graphics.off()


