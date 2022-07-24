rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta1-inflamed_esvd.RData")

tab <- table(regevEpi$Subject, regevEpi$Subject_Disease)
case_subj <- rownames(tab)[which(tab[,"Colitis"] != 0)]
control_subj <- rownames(tab)[which(tab[,"Colitis"] == 0)]

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(length(case_subj))
control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(length(control_subj))

#################################

regevEpi[["pca"]] <- NULL
regevEpi[["umap"]] <- NULL

set.seed(10)
regevEpi <- Seurat::RunPCA(regevEpi, verbose = F)
set.seed(10)
regevEpi <- Seurat::RunUMAP(regevEpi, dims = 1:50)

umap_mat <- regevEpi[["umap"]]@cell.embeddings
n <- nrow(umap_mat)
color_vec <- rep(NA, n)
for(i in 1:length(case_subj)){
  color_vec[regevEpi$Subject == case_subj[i]] <- case_color_palette[i]
}
for(i in 1:length(control_subj)){
  color_vec[regevEpi$Subject == control_subj[i]] <- control_color_palette[i]
}
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_prior_subject.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
title(xlab="UMAP 1", ylab="UMAP 2")
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

color_vec <- rep(NA, n)
color_vec[regevEpi$Subject_Smoking == "Former"] <- rgb(235, 134, 47, maxColorValue = 255)
color_vec[regevEpi$Subject_Smoking == "Never"] <- rgb(184, 54, 220, maxColorValue = 255)
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_prior_smoking.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

color_vec <- rep(NA, n)
color_vec[regevEpi$Subject_Gender == "Female"] <- rgb(235, 134, 47, maxColorValue = 255)
color_vec[regevEpi$Subject_Gender == "Male"] <- rgb(184, 54, 220, maxColorValue = 255)
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_prior_gender.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, name = "Greens")[2:9])(100)
mt_range <- range(regevEpi$percent.mt)
mt_breaks <- seq(mt_range[1], mt_range[2], length.out = 100)
color_vec <- sapply(1:n, function(i){
  color_palette[which.min(abs(regevEpi$percent.mt[i] - mt_breaks))]
})
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_prior_percentmt.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, name = "Greens")[2:9])(100)
lib_range <- range(eSVD_obj$covariates[,"Log_UMI"])
lib_breaks <- seq(lib_range[1], lib_range[2], length.out = 100)
color_vec <- sapply(1:n, function(i){
  color_palette[which.min(abs(eSVD_obj$covariates[i,"Log_UMI"] - lib_breaks))]
})
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_prior_logumi.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

##############################################

x_mat <- eSVD_obj$fit_Second$x_mat
# cc_cov_l2 <- sqrt(sum(eSVD_obj$covariates[,"Subject_Disease_Colitis"]))
# cc_z_l2 <- sqrt(sum(eSVD_obj$fit_Second$z_mat[,"Subject_Disease_Colitis"]))
# mid_val <- sqrt(cc_cov_l2*cc_z_l2)
# x_mat <- cbind(x_mat, eSVD_obj$covariates[,"Subject_Disease_Colitis"]*mid_val/cc_cov_l2)
# apply(x_mat, 2, function(x){sqrt(sum(x^2))})

set.seed(10)
umap_res <- Seurat::RunUMAP(x_mat)
umap_mat <- umap_res@cell.embeddings
n <- nrow(umap_mat)
color_vec <- rep(NA, n)
for(i in 1:length(case_subj)){
  color_vec[regevEpi$Subject == case_subj[i]] <- case_color_palette[i]
}
for(i in 1:length(control_subj)){
  color_vec[regevEpi$Subject == control_subj[i]] <- control_color_palette[i]
}
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_after_subject.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
title(xlab="UMAP 1", ylab="UMAP 2")
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, name = "Greens")[2:9])(100)
mt_range <- range(regevEpi$percent.mt)
mt_breaks <- seq(mt_range[1], mt_range[2], length.out = 100)
color_vec <- sapply(1:n, function(i){
  color_palette[which.min(abs(regevEpi$percent.mt[i] - mt_breaks))]
})
idx <- order(regevEpi$percent.mt, decreasing = F)
set.seed(10)
width <- 3
for(i in 1:n){
  idx[max(1, i-round(n/width)):min(n,i+round(n/width))] <- sample(idx[max(1, i-round(n/width)):min(n,i+round(n/width))])
}
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_after_percentmt.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[idx,1], umap_mat[idx,2], col = color_vec[idx],
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, name = "Greens")[2:9])(100)
lib_range <- range(eSVD_obj$covariates[,"Log_UMI"])
lib_breaks <- seq(lib_range[1], lib_range[2], length.out = 100)
color_vec <- sapply(1:n, function(i){
  color_palette[which.min(abs(eSVD_obj$covariates[i,"Log_UMI"] - lib_breaks))]
})
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_after_logumi.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

color_vec <- rep(NA, n)
color_vec[regevEpi$Subject_Gender == "Female"] <- rgb(235, 134, 47, maxColorValue = 255)
color_vec[regevEpi$Subject_Gender == "Male"] <- rgb(184, 54, 220, maxColorValue = 255)
png("../../../out/fig/main/regevEpi_ta1-inflamed_umap_after_gender.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

###########################################3






