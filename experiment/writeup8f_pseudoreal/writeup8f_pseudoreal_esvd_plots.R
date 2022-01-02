rm(list=ls())
load("../../../../out/writeup8f/writeup8f_sns_pseudoreal_esvd_poisson.RData")
library(Seurat)

# vizualizing batch effects
sns <- Seurat::CreateSeuratObject(counts = t(mat), meta.data = data.frame(esvd_res_full$covariates))
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::ScaleData(sns)
sns <- Seurat::RunPCA(sns, features = rownames(sns), verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns, dims = 1:50)

sns$diagnosis <- metadata$diagnosis
sns$sex <- metadata$sex
sns$individual <- metadata$individual
sns$region <- metadata$region
sns$Seqbatch <- metadata$Seqbatch

covariates <- c("diagnosis", "sex", "individual", "region", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "umap",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real (No adjustment):\n", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8f/sns_pseudoreal_noadjust_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "umap")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real (No adjustment):\n", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8f/sns_pseudoreal_noadjust_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

################################

set.seed(10)
tmp <- Seurat::RunUMAP(esvd_res_full$x_mat)
sns[["esvdumap"]] <- Seurat::CreateDimReducObject(embeddings = tmp@cell.embeddings)

covariates <- c("diagnosis", "sex", "individual", "region", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "esvdumap",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real (eSVD):\n", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8f/sns_pseudoreal_esvd_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "esvdumap")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real (eSVD):\n", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8f/sns_pseudoreal_esvd_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

###############################

png("../../../../out/fig/writeup8f/sns_pseudoreal_histograms.png",
    height = 1200, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
hist(esvd_res_full$b_mat[,"Intercept"],
     xlab = "Intercept", ylab = "Frequency", main = "Intercept",
     breaks = 50)
hist(esvd_res_full$b_mat[,"Log_UMI"],
     xlab = "Log_UMI", ylab = "Frequency", main = "Log_UMI",
     breaks = 50)
hist(esvd_res_full$b_mat[,"diagnosis_ASD"],
     xlab = "diagnosis_ASD", ylab = "Frequency", main = "diagnosis_ASD",
     breaks = 50)
rug(esvd_res_full$b_mat[true_objects$up_idx,"diagnosis_ASD"], col = 2, lwd = 2)
rug(esvd_res_full$b_mat[true_objects$down_idx,"diagnosis_ASD"], col = 3, lwd = 2)
legend("topright", c("Upexpressed DE", "Downexpressed DE"),
       fill = c(2,3))
graphics.off()
