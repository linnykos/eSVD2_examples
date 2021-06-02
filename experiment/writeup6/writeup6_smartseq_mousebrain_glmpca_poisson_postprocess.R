rm(list=ls())

load("../../../../out/writeup6/writeup6_smartseq_mousebrain_glmpca_poisson.RData")

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
plot1 <- plot1 + ggplot2::ggtitle("Mouse Brain (Smartseq)\nGLM-PCA, Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/smartseq_mousebrain_glmpca_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(brain@meta.data)

brain[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(brain, reduction = "glmfactorumap", group.by = "subclass", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Brain (Smartseq)\nGLM-PCA, Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/smartseq_mousebrain_glmpca_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

####################################

source("plotting_func.R")

celltype <- as.factor(brain@meta.data["subclass"][,1])
png("../../../../out/fig/writeup6/smartseq_mousebrain_glmpca_poisson_score.png", height = 1500,
    width = 3000, res = 300, units = "px")
par(mfrow = c(1,2))
plot_scores_heatmap(glmpca_res$factors, membership_vec = celltype,
                    bool_center = T, bool_scale = F,
                    bool_log = T, scaling_power = 2,
                    main = "Scores (GLM-PCA, Poisson, Unscaled)")
plot_scores_heatmap(glmpca_res$factors, membership_vec = celltype,
                    bool_center = T, bool_scale = T,
                    bool_log = T, scaling_power = 2,
                    main = "Scores (GLM-PCA, Poisson, Scaled)")
graphics.off()




