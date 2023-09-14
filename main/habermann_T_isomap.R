rm(list=ls())
library(Seurat)
library(eSVD2)
library(dimRed)

load("../../../out/main/habermann_T_esvd.RData")

# now compute the embeddings
Seurat::DefaultAssay(habermann) <- "RNA"
habermann <- Seurat::NormalizeData(habermann)
habermann <- Seurat::ScaleData(habermann, verbose = F)
set.seed(10)
habermann <- Seurat::RunPCA(habermann, verbose = F)

# https://rdrr.io/cran/dimRed/man/Isomap-class.html
set.seed(10)
pca_mat <- scale(habermann[["pca"]]@cell.embeddings[,1:15])
isomap_original <- dimRed::embed(pca_mat, "Isomap", knn = 30)
isomap_original_mat <- isomap_original@data@data
rownames(isomap_original_mat) <- colnames(habermann)
colnames(isomap_original_mat) <- paste0("Isomap_", 1:2)

#####################

x <- scale(eSVD_obj$fit_Second$x_mat)
set.seed(10)
isomap_esvd <- dimRed::embed(x, "Isomap", knn = 30)
isomap_esvd_mat <- isomap_esvd@data@data
rownames(isomap_esvd_mat) <- colnames(habermann)
colnames(isomap_esvd_mat) <- paste0("eSVDIsomap_", 1:2)

#####################

habermann[["isomap"]] <- Seurat::CreateDimReducObject(isomap_original_mat)
habermann[["esvd"]] <- Seurat::CreateDimReducObject(isomap_esvd_mat)

##########################

# now plot inflamed
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

gender_col_palette <- c(rgb(235, 134, 47, maxColorValue = 255),
                        rgb(184, 54, 220, maxColorValue = 255))
names(gender_col_palette) <- c("F", "M")

smoking_col_palette <- c(rgb(214, 55, 55, maxColorValue = 255),
                        rgb(146, 146, 146, maxColorValue = 255))
names(smoking_col_palette) <- c("Y", "N")

# now plot original
plot1 <- Seurat::DimPlot(habermann, reduction = "isomap",
                         group.by = "Sample_Name",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_original-isomap_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

df <- data.frame(pred =  factor(habermann$Gender)); df <- cbind(df, habermann[["pca"]]@cell.embeddings[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(habermann, reduction = "isomap",
                         group.by = "Gender",
                         cols = gender_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2, 2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_original-isomap_cleaned_by-gender.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

df <- data.frame(pred =  factor(habermann$Tobacco)); df <- cbind(df, habermann[["pca"]]@cell.embeddings[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(habermann, reduction = "isomap",
                         group.by = "Tobacco",
                         cols = smoking_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_original-isomap_cleaned_by-smoking.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

########################
# now plot esvd
plot1 <- Seurat::DimPlot(habermann, reduction = "esvd",
                         group.by = "Sample_Name",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_esvd-isomap_cleaned.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

df <- data.frame(pred =  factor(habermann$Gender)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(habermann, reduction = "esvd",
                         group.by = "Gender",
                         cols = gender_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_esvd-isomap_cleaned_by-gender.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

df <- data.frame(pred =  factor(habermann$Tobacco)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(habermann, reduction = "esvd",
                         group.by = "Tobacco",
                         cols = smoking_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_esvd-isomap_cleaned_by-smoking.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)


print("Done! :)")


