rm(list=ls())
library(Seurat)
library(eSVD2)
library(dimRed)

celltype_names <- c("cyclingta", "entprog", "ta1", "ta2")

gender_col_palette <- c(rgb(235, 134, 47, maxColorValue = 255),
                        rgb(184, 54, 220, maxColorValue = 255))
names(gender_col_palette) <- c("Female", "Male")

smoking_col_palette <- c(rgb(214, 55, 55, maxColorValue = 255),
                         rgb(214, 55, 55, maxColorValue = 255),
                         rgb(146, 146, 146, maxColorValue = 255))
names(smoking_col_palette) <- c("Current", "Former", "Never")

status <- "inflamed"

for(celltype in celltype_names){
  print(celltype)

  load(paste0("../../../out/main/regevEpi_", celltype, "-", status, "_esvd.RData"))

  # now compute the embeddings, one for each of the original data
  set.seed(10)
  regevEpi <- Seurat::RunPCA(regevEpi, verbose = F)
  set.seed(10)
  pca_mat <- scale(regevEpi[["pca"]]@cell.embeddings[,1:15])
  isomap_inf <- dimRed::embed(pca_mat, "Isomap", knn = 30)
  isomap_inf_mat <- isomap_inf@data@data
  rownames(isomap_inf_mat) <- colnames(regevEpi)
  colnames(isomap_inf_mat) <- paste0("Isomap_", 1:2)

  ##########################

  # now compute the embeddings, one for each of the esvd fit
  x <- scale(eSVD_obj$fit_Second$x_mat)
  set.seed(10)
  isomap_eSVD_obj <- dimRed::embed(x, "Isomap", knn = 30)
  isomap_eSVD_obj_mat <- isomap_eSVD_obj@data@data
  rownames(isomap_eSVD_obj_mat) <- colnames(regevEpi)
  colnames(isomap_eSVD_obj_mat) <- paste0("eSVDIsomap_", 1:2)

  #####################

  regevEpi[["isomap"]] <- Seurat::CreateDimReducObject(isomap_inf_mat)
  regevEpi[["esvd"]] <- Seurat::CreateDimReducObject(isomap_eSVD_obj_mat)

  ##########################

  # now plot inflamed
  tab_mat <- table(regevEpi$Subject, regevEpi$Subject_Disease)
  case_indiv <- rownames(tab_mat)[which(tab_mat[,"Colitis"] > 0)]
  num_cases <- length(case_indiv)
  case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                      rgb(244, 84, 84, maxColorValue = 255)))(num_cases)
  names(case_color_palette) <- case_indiv

  control_indiv <- rownames(tab_mat)[which(tab_mat[,"HC"] > 0)]
  num_controls <- length(control_indiv)
  control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                         rgb(27, 198, 245, maxColorValue = 255)))(num_controls)
  names(control_color_palette) <- control_indiv
  col_palette <- c(case_color_palette, control_color_palette)

  ###################

  plot1 <- Seurat::DimPlot(regevEpi, reduction = "isomap",
                           group.by = "Subject",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_original-isomap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  df <- data.frame(pred =  factor(regevEpi$Subject_Gender)); df <- cbind(df, regevEpi[["pca"]]@cell.embeddings[,1:15])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null

  plot1 <- Seurat::DimPlot(regevEpi, reduction = "isomap",
                           group.by = "Subject_Gender",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2, 2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_original-isomap_cleaned_by-gender.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  df <- data.frame(pred =  factor(regevEpi$Subject_Smoking)); df <- cbind(df, regevEpi[["pca"]]@cell.embeddings[,1:15])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null

  plot1 <- Seurat::DimPlot(regevEpi, reduction = "isomap",
                           group.by = "Subject_Smoking",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2, 2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_original-isomap_cleaned_by-smoking.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  ##################################
  ##################################
  ##################################

  plot1 <- Seurat::DimPlot(regevEpi, reduction = "esvd",
                           group.by = "Subject",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_esvd-isomap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  df <- data.frame(pred =  factor(regevEpi$Subject_Gender)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:15])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null

  plot1 <- Seurat::DimPlot(regevEpi, reduction = "esvd",
                           group.by = "Subject_Gender",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2, 2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_esvd-isomap_cleaned_by-gender.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  df <- data.frame(pred =  factor(regevEpi$Subject_Smoking)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:15])
  fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
  null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
  deviance_current <- stats::deviance(fitted_model)
  deviance_null <- stats::deviance(null_model)
  r2 <- 1 - deviance_current/deviance_null

  plot1 <- Seurat::DimPlot(regevEpi, reduction = "esvd",
                           group.by = "Subject_Smoking",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2, 2)))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_esvd-isomap_cleaned_by-smoking.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

}

print("Done! :)")


