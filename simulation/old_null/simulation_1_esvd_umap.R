rm(list=ls())
library(Seurat)
library(SummarizedExperiment)

#########################################

load("../eSVD2_examples/simulation/simulation_1.RData")

tab <- table(seurat_obj$individual, seurat_obj$cc)
case_subj <- rownames(tab)[which(tab[,"1"] != 0)]
control_subj <- rownames(tab)[which(tab[,"1"] == 0)]

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(length(case_subj))
control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(length(control_subj))

umap_mat <- seurat_obj[["umap"]]@cell.embeddings
n <- nrow(umap_mat)
color_vec <- rep(NA, n)
for(i in 1:length(case_subj)){
  color_vec[seurat_obj$individual == case_subj[i]] <- case_color_palette[i]
}
for(i in 1:length(control_subj)){
  color_vec[seurat_obj$individual == control_subj[i]] <- control_color_palette[i]
}
png("../../out/fig/simulation/simulation_1_umap-observed.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
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

de_idx <- which(abs(true_teststat_vec) >= 2)

svd_res <- irlba::irlba(mean_mat[,de_idx], nv = 11)
pca_mat <- svd_res$u %*% diag(svd_res$d)
umap_mat2 <- Seurat::RunUMAP(pca_mat)
umap_mat2 <- umap_mat2@cell.embeddings
png("../../out/fig/simulation/simulation_1_umap-latent.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat2[,1], umap_mat2[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
title(xlab="UMAP 1", ylab="UMAP 2")
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

