rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

DefaultAssay(sns) <- "RNA"
sns <- Seurat::ScaleData(sns)
sns <- Seurat::RunPCA(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns,
                       reduction = 'pca',
                       dims = 1:30)

indiv_tab <- table(sns$individual, sns$diagnosis)
case_indiv <- rownames(indiv_tab)[which(indiv_tab[,"ASD"] != 0)]
control_indiv <- rownames(indiv_tab)[which(indiv_tab[,"Control"] != 0)]

case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(length(case_indiv))
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(length(control_indiv))
col_vec <- c(case_color_palette, control_color_palette)
names(col_vec) <- c(case_indiv, control_indiv)

plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                         group.by = "individual",
                         cols = col_vec)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_umap_individual_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

######################################

plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                         group.by = "diagnosis")
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_umap_diagnosis_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

########################

col_vec <- c("olivedrab4", "gray")
names(col_vec) <- c("F", "M")

plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                         group.by = "sex",
                         cols = col_vec) + Seurat::NoAxes()
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_umap_sex_cleaned.png"),
                plot1, device = "png", width = 2, height = 2, units = "in")

########################

plot1 <- Seurat::FeaturePlot(sns,
                             features = "age",
                             reduction = "umap")
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_umap_age_cleaned.png"),
                plot1, device = "png", width = 2, height = 2, units = "in")

########################

col_vec <- c("olivedrab4", "gray")
names(col_vec) <- c("F", "M")

plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                         group.by = "sex",
                         cols = col_vec)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_umap_sex_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

########################

col_vec <- c("olivedrab4", "gray")
names(col_vec) <- c("ACC", "PFC")

plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                         group.by = "region",
                         cols = col_vec)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_umap_region_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

########################

plot1 <- Seurat::FeaturePlot(sns,
                             features = "post.mortem.hours",
                             reduction = "umap")
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_umap_postmortem_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")



