rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/sns_layer23_esvd.RData")
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]

# now compute the embeddings
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::ScaleData(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunPCA(sns, verbose = F)

######################

x <- eSVD_obj$fit_Second$x_mat
colnames(x) <- paste0("eSVD_", 1:ncol(x))
sns[["esvd"]] <- Seurat::CreateDimReducObject(x)

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

plot1 <- Seurat::DimPlot(sns, reduction = "pca",
                         group.by = "individual",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup13b/sns_pca_layer23.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

# https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html
df <- data.frame(pred =  factor(sns$sex)); df <- cbind(df, sns[["pca"]]@cell.embeddings[,1:30])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(sns, reduction = "pca",
                         group.by = "sex",
                         cols = gender_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_pca_layer23_by-gender.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

# https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html
df <- data.frame(pred =  factor(sns$region)); df <- cbind(df, sns[["pca"]]@cell.embeddings[,1:30])
fitted_model <- stats::glm(pred ~ ., data = df, family = binomial)
null_model <- stats::glm(pred ~ 1, data = df, family = binomial)
deviance_current <- stats::deviance(fitted_model)
deviance_null <- stats::deviance(null_model)
r2 <- 1 - deviance_current/deviance_null
r2

plot1 <- Seurat::DimPlot(sns, reduction = "pca",
                         group.by = "region",
                         cols = region_col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_pca_layer23_by-region.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

##################################

plot1 <- Seurat::DimPlot(sns, reduction = "esvd",
                         group.by = "individual",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup13b/sns_esvd_layer23.png"),
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
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_esvd_layer23_by-gender.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)

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
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/Writeup13b/sns_esvd_layer23_by-region.png",
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)
