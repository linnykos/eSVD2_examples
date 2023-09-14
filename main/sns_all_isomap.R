rm(list=ls())
library(Seurat)
library(eSVD2)
library(dimRed)

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
  plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                           group.by = "individual",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-isomap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  df <- data.frame(pred =  factor(sns$sex)); df <- cbind(df, sns[["pca"]]@cell.embeddings[,1:30])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null
  r2

  plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                           group.by = "sex",
                           cols = gender_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-isomap_by-gender.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)
  
  plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                           group.by = "sex",
                           cols = gender_col_palette,
                           pt.size = 1.25)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-isomap_cleaned_by-gender.png"),
                  plot1, device = "png", width = 2, height = 2, units = "in",
                  dpi = 500)

  df <- data.frame(pred =  factor(sns$region)); df <- cbind(df, sns[["pca"]]@cell.embeddings[,1:30])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null
  r2

  plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                           group.by = "region",
                           cols = region_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-isomap_by-region.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)
  
  plot1 <- Seurat::DimPlot(sns, reduction = "isomap",
                           group.by = "region",
                           cols = region_col_palette,
                           pt.size = 1.25)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_original-isomap_cleaned_by-region.png"),
                  plot1, device = "png", width = 2, height = 2, units = "in",
                  dpi = 500)

  ########################
  # now plot esvd
  plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                           group.by = "individual",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-isomap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  df <- data.frame(pred =  factor(sns$sex)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:30])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null
  r2

  plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                           group.by = "sex",
                           cols = gender_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-isomap_by-gender.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)
  
  plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                           group.by = "sex",
                           cols = gender_col_palette,
                           pt.size = 1.25)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-isomap_cleaned_by-gender.png"),
                  plot1, device = "png", width = 2, height = 2, units = "in",
                  dpi = 500)

  df <- data.frame(pred =  factor(sns$region)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:30])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null
  r2

  plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                           group.by = "region",
                           cols = region_col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-isomap_by-region.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)
  
  plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                           group.by = "region",
                           cols = region_col_palette,
                           pt.size = 1.25)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_", celltype, "_esvd-isomap_cleaned_by-region.png"),
                  plot1, device = "png", width = 2, height = 2, units = "in",
                  dpi = 500)

}

print("Done! :)")


