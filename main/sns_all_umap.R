rm(list=ls())
library(Seurat)
library(eSVD2)

file_vec <- c("../../../out/main/sns_astpp_esvd.RData",
              "../../../out/main/sns_endothelial_esvd.RData",
              "../../../out/main/sns_insst_esvd.RData",
              "../../../out/main/sns_invip_esvd.RData",
              "../../../out/main/sns_layer4_esvd.RData",
              "../../../out/main/sns_layer23_esvd.RData",
              "../../../out/main/sns_layer56_esvd.RData",
              "../../../out/main/sns_layer56cc_esvd.RData",
              "../../../out/main/sns_microglia_esvd.RData",
              "../../../out/main/sns_oligo_esvd.RData",
              "../../../out/main/sns_opc_esvd.RData")
names(file_vec) <- c("astpp", "endothelial", "insst", "invip", "layer4", "layer23",
                     "layer56", "layer56cc", "microglia", "oligo", "opc")

for(kk in 1:length(file_vec)){
  file <- file_vec[kk]
  print(file)
  load(file)
  celltype <- names(file_vec)[kk]

  # now compute the embeddings
  Seurat::DefaultAssay(sns) <- "RNA"
  sns <- Seurat::NormalizeData(sns)
  sns <- Seurat::ScaleData(sns, verbose = F)
  set.seed(10)
  sns <- Seurat::RunPCA(sns, verbose = F)
  sns <- Seurat::RunUMAP(sns, dims = 1:30,
                         verbose = F)

  score_mat <- scale(eSVD_obj$fit_Second$x_mat)
  sns[["esvdumap"]] <- Seurat::RunUMAP(score_mat,
                                       verbose = F)

  ##########################

  # now plot inflamed
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


  # now plot original
  plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                           group.by = "individual",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-umap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                           group.by = "sex",
                           cols = gender_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-umap_cleaned_by-gender.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                           group.by = "Seqbatch")
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-umap_cleaned_by-seqbatch.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                           group.by = "region",
                           cols = region_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-umap_cleaned_by-region.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::FeaturePlot(sns, reduction = "umap",
                               features = "age")
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-umap_cleaned_by-age.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  ########################
  # now plot esvd
  plot1 <- Seurat::DimPlot(sns, reduction = "esvdumap",
                           group.by = "individual",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-umap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::DimPlot(sns, reduction = "esvdumap",
                           group.by = "sex",
                           cols = gender_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-umap_cleaned_by-gender.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::DimPlot(sns, reduction = "esvdumap",
                           group.by = "Seqbatch")
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-umap_cleaned_by-seqbatch.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::DimPlot(sns, reduction = "esvdumap",
                           group.by = "region",
                           cols = region_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-umap_cleaned_by-region.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  plot1 <- Seurat::FeaturePlot(sns, reduction = "esvdumap",
                               features = "age")
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-umap_cleaned_by-age.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

}

print("Done! :)")


