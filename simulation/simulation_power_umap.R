rm(list=ls())
library(Seurat)
library(SummarizedExperiment)

#########################################

load("../eSVD2_examples/simulation/simulation_power.RData")

tab <- table(seurat_obj$individual, seurat_obj$cc)
case_subj <- rownames(tab)[which(tab[,"1"] != 0)]
control_subj <- rownames(tab)[which(tab[,"1"] == 0)]

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(length(case_subj))
control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(length(control_subj))

umap_mat <- seurat_obj[["isomap"]]@cell.embeddings
n <- nrow(umap_mat)
color_vec <- rep(NA, n)
for(i in 1:length(case_subj)){
  color_vec[seurat_obj$individual == case_subj[i]] <- case_color_palette[i]
}
for(i in 1:length(control_subj)){
  color_vec[seurat_obj$individual == control_subj[i]] <- control_color_palette[i]
}

set.seed(10)
shuff_idx <- sample(1:nrow(umap_mat))

png("../../out/fig/simulation/simulation_power_isomap-observed.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[shuff_idx,1], umap_mat[shuff_idx,2], col = color_vec[shuff_idx],
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
title(xlab="UMAP 1", ylab="UMAP 2")
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

############################

tmp_mat <- nat_mat
tmp_mat2 <- tmp_mat

for(j in 1:ncol(tmp_mat)){
  df <- data.frame(cbind(tmp_mat[,j], covariates))
  colnames(df)[1] <- "y"
  lm_res <- stats::lm(y ~ . - 1, data = df)

  tmp_mat2[,j] <- stats::residuals(lm_res) + stats::coef(lm_res)["cc"]*covariates[,"cc"]
}

mean_mat <- exp(tmp_mat2)
mean_mat <- pmin(mean_mat, 10)
mean_mat <- scale(mean_mat)

de_idx <- which(true_fdr_vec <= 0.05)
mean_mat <- mean_mat[,de_idx]
col_idx <- which(apply(mean_mat, 2, function(x){any(is.na(x))}))
if(length(col_idx) > 0) mean_mat <- mean_mat[,-col_idx]

set.seed(10)
svd_res <- irlba::irlba(mean_mat, nv = 5)
pca_mat <- svd_res$u %*% diag(svd_res$d)

isomap_res <- dimRed::embed(pca_mat, "Isomap", knn = 10)
isomap_2 <- isomap_res@data@data
rownames(isomap_2) <- colnames(seurat_obj)
colnames(isomap_2) <- paste0("Isomap_", 1:2)

set.seed(10)
shuff_idx <- sample(1:nrow(isomap_2))

png("../../out/fig/simulation/simulation_power_isomap-latent.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(isomap_2[shuff_idx,1], isomap_2[shuff_idx,2],
     xlim = c(min(isomap_2[,1]), 6),
     col = color_vec[shuff_idx],
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
title(xlab="UMAP 1", ylab="UMAP 2")
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

