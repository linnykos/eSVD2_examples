rm(list=ls())
library(Seurat)
library(eSVD2)
library(dimRed)

load("../../../out/main/adams_T_esvd.RData")

# now compute the embeddings
Seurat::DefaultAssay(adams) <- "RNA"
adams <- Seurat::NormalizeData(adams)
adams <- Seurat::ScaleData(adams, verbose = F)
set.seed(10)
adams <- Seurat::RunPCA(adams, verbose = F)

# https://rdrr.io/cran/dimRed/man/Isomap-class.html
set.seed(10)
pca_mat <- scale(adams[["pca"]]@cell.embeddings[,1:15])
isomap_original <- dimRed::embed(pca_mat, "Isomap", knn = 30)
isomap_original_mat <- isomap_original@data@data
rownames(isomap_original_mat) <- colnames(adams)
colnames(isomap_original_mat) <- paste0("Isomap_", 1:2)

#####################

x <- scale(eSVD_obj$fit_Second$x_mat)
set.seed(10)
isomap_esvd <- dimRed::embed(x, "Isomap", knn = 30)
isomap_esvd_mat <- isomap_esvd@data@data
rownames(isomap_esvd_mat) <- colnames(adams)
colnames(isomap_esvd_mat) <- paste0("eSVDIsomap_", 1:2)

#####################

adams[["isomap"]] <- Seurat::CreateDimReducObject(isomap_original_mat)
adams[["esvd"]] <- Seurat::CreateDimReducObject(isomap_esvd_mat)

##########################

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")

# now plot inflamed
tab_mat <- table(adams$Subject_Identity, adams$Disease_Identity)
case_indiv <- rownames(tab_mat)[which(tab_mat[,"IPF"] > 0)]
num_cases <- length(case_indiv)
case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(num_cases)
names(case_color_palette) <- case_indiv

control_indiv <- rownames(tab_mat)[which(tab_mat[,"Control"] > 0)]
num_controls <- length(control_indiv)
control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(num_controls)
names(control_color_palette) <- control_indiv
col_palette <- c(case_color_palette, control_color_palette)

gender_col_palette <- c(rgb(235, 134, 47, maxColorValue = 255),
                        rgb(184, 54, 220, maxColorValue = 255))
names(gender_col_palette) <- c("Female", "Male")

smoking_col_palette <- c(rgb(214, 55, 55, maxColorValue = 255),
                        rgb(146, 146, 146, maxColorValue = 255))
names(smoking_col_palette) <- c("Yes", "No")

# now plot original
isomap_mat <- adams[["isomap"]]@cell.embeddings
n <- nrow(isomap_mat)
color_vec <- rep(NA, n)
for(i in 1:length(case_indiv)){
  color_vec[adams$Subject_Identity == case_indiv[i]] <- case_color_palette[i]
}
for(i in 1:length(control_indiv)){
  color_vec[adams$Subject_Identity == control_indiv[i]] <- control_color_palette[i]
}
png("../../../out/fig/main/adams_T_original-isomap_cleaned.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(isomap_mat[,1], isomap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
title(xlab="Isomap 1", ylab="Isomap 2")
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

df <- data.frame(pred =  factor(adams$Gender)); df <- cbind(df, adams[["pca"]]@cell.embeddings[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(adams, reduction = "isomap",
                         group.by = "Gender",
                         cols = gender_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2, 2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/adams_T_original-isomap_by-gender.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

color_vec <- rep(NA, n)
color_vec[adams$Gender == "Female"] <- rgb(235, 134, 47, maxColorValue = 255)
color_vec[adams$Gender == "Male"] <- rgb(184, 54, 220, maxColorValue = 255)
png("../../../out/fig/main/adams_T_original-isomap_cleaned_by-gender.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(isomap_mat[,1], isomap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16, cex = 0.5)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

df <- data.frame(pred =  factor(adams$Tobacco)); df <- cbind(df, adams[["pca"]]@cell.embeddings[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(adams, reduction = "isomap",
                         group.by = "Tobacco",
                         cols = smoking_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/adams_T_original-isomap_by-smoking.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

color_vec <- rep(NA, n)
color_vec[adams$Tobacco == "Yes"] <- rgb(214, 55, 55, maxColorValue = 255)
color_vec[adams$Tobacco == "No"] <- rgb(146, 146, 146, maxColorValue = 255)
png("../../../out/fig/main/adams_T_original-isomap_cleaned_by-smoking.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(isomap_mat[,1], isomap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16, cex = 0.5)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

########################
# now plot esvd
isomap_mat <- adams[["esvd"]]@cell.embeddings
n <- nrow(isomap_mat)
color_vec <- rep(NA, n)
for(i in 1:length(case_indiv)){
  color_vec[adams$Subject_Identity == case_indiv[i]] <- case_color_palette[i]
}
for(i in 1:length(control_indiv)){
  color_vec[adams$Subject_Identity == control_indiv[i]] <- control_color_palette[i]
}
png("../../../out/fig/main/adams_T_esvd-isomap_cleaned.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(isomap_mat[,1], isomap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
title(xlab="Isomap 1", ylab="Isomap 2")
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

df <- data.frame(pred =  factor(adams$Gender)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(adams, reduction = "esvd",
                         group.by = "Gender",
                         cols = gender_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/adams_T_esvd-isomap_by-gender.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

color_vec <- rep(NA, n)
color_vec[adams$Gender == "Female"] <- rgb(235, 134, 47, maxColorValue = 255)
color_vec[adams$Gender == "Male"] <- rgb(184, 54, 220, maxColorValue = 255)
png("../../../out/fig/main/adams_T_esvd-isomap_cleaned_by-gender.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(isomap_mat[,1], isomap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16, cex = 0.5)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

df <- data.frame(pred =  factor(adams$Tobacco)); df <- cbind(df, eSVD_obj$fit_Second$x_mat[,1:15])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(adams, reduction = "esvd",
                         group.by = "Tobacco",
                         cols = smoking_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("R2: ", round(r2,2)))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/fig/main/adams_T_esvd-isomap_by-smoking.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

color_vec <- rep(NA, n)
color_vec[adams$Tobacco == "Yes"] <- rgb(214, 55, 55, maxColorValue = 255)
color_vec[adams$Tobacco == "No"] <- rgb(146, 146, 146, maxColorValue = 255)
png("../../../out/fig/main/adams_T_esvd-isomap_cleaned_by-smoking.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(isomap_mat[,1], isomap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16, cex = 0.5)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

print("Done! :)")


