rm(list=ls())
library(Seurat)
load("../../../out/main/adams_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(adams))
keep_vec[which(adams$Manuscript_Identity == "T")] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
adams

adams[["pca"]] <- NULL
adams[["umap"]] <- NULL

set.seed(10)
adams <- Seurat::NormalizeData(adams,
                               normalization.method = "LogNormalize", scale.factor = 10000)
adams <- Seurat::FindVariableFeatures(adams,
                                      selection.method = "vst", nfeatures = 5000)
adams <- Seurat::ScaleData(adams)

set.seed(10)
adams <- Seurat::RunPCA(adams, verbose = F)
set.seed(10)
adams <- Seurat::RunUMAP(adams, dims = 1:50)

save(adams, date_of_run, session_info,
     file = "../../../out/main/adams_T_preprocessed.RData")

#############

tab_mat <- table(adams$Subject_Identity, adams$Disease_Identity)
case_indiv <- rownames(tab_mat)[which(tab_mat[,"IPF"] > 0)]
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

plot1 <- Seurat::DimPlot(adams, reduction = "umap",
                         group.by = "Subject_Identity",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/adams_T_umap_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)


