rm(list=ls())

load("../../../../out/writeup6/writeup6_10x_mouseretinal_glmpca_nb.RData")

glm_assay <- Seurat::CreateAssayObject(counts = pred_mat)
retinal[["glm"]] <- glm_assay

Seurat::DefaultAssay(retinal) <- "glm"
retinal <- Seurat::NormalizeData(retinal, normalization.method = "LogNormalize", scale.factor = 10000)
retinal[["glm"]]@var.features <- rownames(retinal)
retinal <-  Seurat::ScaleData(retinal)
retinal <- Seurat::RunPCA(retinal, features = Seurat::VariableFeatures(retinal), verbose = F,
                           reduction.name = "glmpca")

set.seed(10)
retinal <- Seurat::RunUMAP(retinal, reduction = "glmpca", dims = 1:30, reduction.name = "glmumap")

plot1 <- Seurat::DimPlot(retinal, reduction = "glmumap", group.by = "Cluster", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Retinal (10x)\nGLM-PCA, Full, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/10x_mouseretinal_glmpca_full_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(retinal@meta.data)

retinal[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(retinal, reduction = "glmfactorumap", group.by = "Cluster", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Retinal (10x)\nGLM-PCA, Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/10x_mouseretinal_glmpca_factor_nb_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")


####################################

source("plotting_func.R")

celltype <- as.factor(retinal@meta.data["Cluster"][,1])
png("../../../../out/fig/writeup6/10x_mouseretinal_glmpca_nb_score.png", height = 1500,
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

