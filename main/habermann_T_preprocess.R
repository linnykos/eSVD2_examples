rm(list=ls())
library(Seurat)

load("../../../out/main/habermann_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$celltype == "T Cells")] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)

Seurat::DefaultAssay(habermann) <- "RNA"
habermann[["SCT"]] <- NULL
habermann[["umap.rna"]] <- NULL
set.seed(10)
habermann <- Seurat::SCTransform(habermann, variable.features.n = 3000)

set.seed(10); habermann <- Seurat::RunPCA(habermann, verbose = F)
set.seed(10)
habermann <- Seurat::RunUMAP(habermann, dims = 1:50)

save(habermann, date_of_run, session_info,
     file = "../../../out/main/habermann_T_preprocessed.RData")

#############

tab_mat <- table(habermann$Sample_Name, habermann$Diagnosis)
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

plot1 <- Seurat::DimPlot(habermann, reduction = "umap",
                         group.by = "Sample_Name",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_umap_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)




