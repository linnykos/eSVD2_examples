rm(list=ls())
library(Seurat)
library(eSVD2)
library(dimRed)

load("../../../../out/main/sns_layer23_esvd.RData")

# now compute the embeddings
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::ScaleData(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunPCA(sns, verbose = F)

# https://rdrr.io/cran/dimRed/man/Isomap-class.html
set.seed(10)
pca_mat <- scale(sns[["pca"]]@cell.embeddings[,1:30])
isomap_original <- dimRed::embed(pca_mat, "Isomap", knn = 30)
isomap_original_mat <- isomap_original@data@data
rownames(isomap_original_mat) <- colnames(sns)
colnames(isomap_original_mat) <- paste0("Isomap_", 1:2)

#####################

x <- scale(eSVD_obj$fit_Second$x_mat)
set.seed(10)
isomap_esvd <- dimRed::embed(x, "Isomap", knn = 30)
isomap_esvd_mat <- isomap_esvd@data@data
rownames(isomap_esvd_mat) <- colnames(sns)
colnames(isomap_esvd_mat) <- paste0("eSVDIsomap_", 1:2)

#####################

sns[["isomap"]] <- Seurat::CreateDimReducObject(isomap_original_mat)
sns[["esvd"]] <- Seurat::CreateDimReducObject(isomap_esvd_mat)

######################

tab_mat <- table(sns$individual, sns$diagnosis)
case_indiv <- rownames(tab_mat)[which(tab_mat[,"ASD"] > 0)]
num_cases <- length(case_indiv)
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(num_cases)
names(case_color_palette) <- case_indiv

control_indiv <- rownames(tab_mat)[which(tab_mat[,"Control"] > 0)]
num_controls <- length(control_indiv)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(num_controls)
names(control_color_palette) <- control_indiv
col_palette <- c(case_color_palette, control_color_palette)

gender_col_palette <- c(rgb(235, 134, 47, maxColorValue = 255),
                        rgb(184, 54, 220, maxColorValue = 255))
names(gender_col_palette) <- c("F", "M")

region_col_palette <- c(rgb(217, 198, 159, maxColorValue = 255),
                        rgb(43, 85, 99, maxColorValue = 255))
names(region_col_palette) <- c("ACC", "PFC")

#########################

plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                         group.by = "individual",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup13b/sns_isomap_layer23.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                         group.by = "sex",
                         cols = gender_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_isomap_layer23_by-gender.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                         group.by = "region",
                         cols = region_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_isomap_layer23_by-region.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

##################################

plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                         group.by = "individual",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup13b/sns_esvd-isomap_layer23.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                         group.by = "sex",
                         cols = gender_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_esvd-isomap_layer23_by-gender.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                         group.by = "region",
                         cols = region_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_esvd-isomap_layer23_by-region.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

