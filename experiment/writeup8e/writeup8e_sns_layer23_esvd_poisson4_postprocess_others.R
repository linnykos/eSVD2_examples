
# try the two naive methods -- MAST and also SCTransform
library(Seurat); library(MAST)
sns <- Seurat::CreateSeuratObject(counts = t(mat), meta.data = metadata)
Seurat::Idents(sns) <- "diagnosis"
set.seed(10)
mast_res <- Seurat::FindMarkers(sns, ident.1 = "ASD", ident.2 = "Control", test.use = "MAST",
                                logfc.threshold = 0, min.pct = 0)

case_idx <- which(esvd_res_full$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res_full$covariates[,"diagnosis_ASD"] == 0)
x_vec <- sapply(1:ncol(mat), function(j){
  log2(mean(expm1(mat[case_idx,j]))) - log2(mean(expm1(mat[control_idx,j])))
})

hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]
hk_idx <- which(rownames(mast_res) %in% hk_genes)
de_idx <- which(rownames(mast_res) %in% de_genes)
col_vec <- rep(rgb(0.5,0.5,0.5,0.5), nrow(mast_res))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
idx <- sample(c(de_idx, hk_idx))
png("../../out/fig/writeup8e/sns_layer23_volcano_mast1.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-150,150), ylim = c(0, max(-log10(mast_res$p_val))), bty = "n",
     main = "Volcano plot for Layer 2/3\n(MAST, unadjusted)",
     xlab = "Log2 fold change", ylab = "-Log10(P value)")
points(x = mast_res$avg_log2FC,
       y = -log10(mast_res$p_val),
       pch = 16, col = col_vec)
points(x = mast_res$avg_log2FC[idx],
       y = -log10(mast_res$p_val[idx]),
       pch = 16, col = "white", cex = 1.5)
points(x = mast_res$avg_log2FC[idx],
       y = -log10(mast_res$p_val[idx]),
       pch = 16, col = col_vec[idx])
legend("topright", c("Published DE gene", "Housekeeping gene", "Other"),
       fill = c(2,3,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()

##########################

sns <- Seurat::CreateSeuratObject(counts = t(mat), meta.data = metadata)
Seurat::Idents(sns) <- "diagnosis"
set.seed(10)
de_res <- Seurat::FindMarkers(sns, ident.1 = "ASD", ident.2 = "Control",
                              test.use = "negbinom",
                              logfc.threshold = 0, min.pct = 0)
hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]
hk_idx <- which(rownames(de_res) %in% hk_genes)
de_idx <- which(rownames(de_res) %in% de_genes)
col_vec <- rep(rgb(0.5,0.5,0.5,0.5), nrow(de_res))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
idx <- sample(c(de_idx, hk_idx))
png("../../out/fig/writeup8e/sns_layer23_volcano_negbinom1.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = range(de_res$avg_log2FC), ylim = c(0, max(-log10(de_res$p_val))), bty = "n",
     main = "Volcano plot for Layer 2/3\n(Neg-binom, unadjusted)",
     xlab = "Log2 fold change", ylab = "-Log10(P value)")
points(x = de_res$avg_log2FC,
       y = -log10(de_res$p_val),
       pch = 16, col = col_vec)
points(x = de_res$avg_log2FC[idx],
       y = -log10(de_res$p_val[idx]),
       pch = 16, col = "white", cex = 1.5)
points(x = de_res$avg_log2FC[idx],
       y = -log10(de_res$p_val[idx]),
       pch = 16, col = col_vec[idx])
legend("topright", c("Published DE gene", "Housekeeping gene", "Other"),
       fill = c(2,3,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()

##########################

sns2 <- Seurat::CreateSeuratObject(counts = t(mat), meta.data = as.data.frame(esvd_res_full$covariates))
sns2 <- Seurat::NormalizeData(sns2)
vars.to.regress <- colnames(esvd_res_full$covariates)
vars.to.regress <- vars.to.regress[!vars.to.regress %in% c("Intercept", "Log_UMI", "diagnosis_ASD", "nFeature_RNA")]
sns2 <- Seurat::ScaleData(sns2, vars.to.regress = vars.to.regress)
sns2$diagnosis <- metadata[,"diagnosis"]
Seurat::Idents(sns2) <- "diagnosis"
de_res <- Seurat::FindMarkers(sns2, ident.1 = "ASD", ident.2 = "Control",
                              slot = "scale.data",
                              logfc.threshold = 0, min.pct = 0)
p_vec <- de_res$p_val
p_vec[p_vec == 0] <- min(p_vec[p_vec != 0])/2
hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]
hk_idx <- which(rownames(de_res) %in% hk_genes)
de_idx <- which(rownames(de_res) %in% de_genes)
col_vec <- rep(rgb(0.5,0.5,0.5,0.5), nrow(de_res))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
idx <- sample(c(de_idx, hk_idx))
png("../../out/fig/writeup8e/sns_layer23_volcano_lognorm2.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = range(de_res$avg_diff), ylim = c(0, max(-log10(de_res$p_val))), bty = "n",
     main = "Volcano plot for Layer 2/3\n(Log-norm, adjusted)",
     xlab = "Mean change", ylab = "-Log10(P value)")
points(x = de_res$avg_diff,
       y = -log10(p_vec),
       pch = 16, col = col_vec)
points(x = de_res$avg_diff[idx],
       y = -log10(p_vec[idx]),
       pch = 16, col = "white", cex = 1.5)
points(x = de_res$avg_diff[idx],
       y = -log10(p_vec[idx]),
       pch = 16, col = col_vec[idx])
legend("topright", c("Published DE gene", "Housekeeping gene", "Other"),
       fill = c(2,3,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()

sns2 <- Seurat::RunPCA(sns2, features = rownames(sns2), verbose = F)
set.seed(10)
sns2 <- Seurat::RunUMAP(sns2, dims = 1:50)
sns2$diagnosis <- metadata$diagnosis
sns2$sex <- metadata$sex
sns2$individual <- metadata$individual
sns2$region <- metadata$region
sns2$Capbatch <- metadata$Capbatch
sns2$Seqbatch <- metadata$Seqbatch

covariates <- c("diagnosis", "sex", "individual", "region", "Capbatch", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns2, reduction = "umap",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via Log-normalization\n(Adjusted): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../out/fig/writeup8e/sns_layer23_lognorm2_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns2,
                               features = covariate,
                               reduction = "umap")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via Log-normalization\n(Adjusted): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../out/fig/writeup8e/sns_layer23_lognorm2_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

##########################

sns2 <- Seurat::CreateSeuratObject(counts = t(mat), meta.data = as.data.frame(esvd_res_full$covariates))
vars.to.regress <- colnames(esvd_res_full$covariates)
vars.to.regress <- vars.to.regress[!vars.to.regress %in% c("Intercept", "Log_UMI", "diagnosis_ASD", "nFeature_RNA")]
sns2 <- Seurat::SCTransform(sns2, vars.to.regress = vars.to.regress, verbose = T)
sns2$diagnosis <- metadata[,"diagnosis"]
Seurat::Idents(sns2) <- "diagnosis"
Seurat::DefaultAssay(sns2) <- "SCT"
de_res <- Seurat::FindMarkers(sns2, ident.1 = "ASD", ident.2 = "Control",
                              slot = "scale.data",
                              logfc.threshold = 0, min.pct = 0)
p_vec <- de_res$p_val
p_vec[p_vec == 0] <- min(p_vec[p_vec != 0])/2
hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]
hk_idx <- which(rownames(de_res) %in% hk_genes)
de_idx <- which(rownames(de_res) %in% de_genes)
col_vec <- rep(rgb(0.5,0.5,0.5,0.5), nrow(de_res))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
idx <- sample(c(de_idx, hk_idx))
png("../../out/fig/writeup8e/sns_layer23_volcano_sctransform2.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = range(de_res$avg_diff), ylim = c(0, max(-log10(p_vec))), bty = "n",
     main = "Volcano plot for Layer 2/3\n(SCTransform, adjusted)",
     xlab = "Mean change", ylab = "-Log10(P value)")
points(x = de_res$avg_diff,
       y = -log10(p_vec),
       pch = 16, col = col_vec)
points(x = de_res$avg_diff[idx],
       y = -log10(p_vec[idx]),
       pch = 16, col = "white", cex = 1.5)
points(x = de_res$avg_diff[idx],
       y = -log10(p_vec[idx]),
       pch = 16, col = col_vec[idx])
legend("topright", c("Published DE gene", "Housekeeping gene", "Other"),
       fill = c(2,3,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()


###################################333333333

# vizualizing batch effects
sns <- Seurat::CreateSeuratObject(counts = t(mat), meta.data = metadata)
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::ScaleData(sns)
sns <- Seurat::RunPCA(sns, features = rownames(sns), verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns, dims = 1:50)

covariates <- c("diagnosis", "sex", "individual", "region", "Capbatch", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) (No adjustment):\n", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../out/fig/writeup8e/sns_layer23_noadjust_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "umap")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via (No adjustment):\n", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../out/fig/writeup8e/sns_layer23_noadjust_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}
