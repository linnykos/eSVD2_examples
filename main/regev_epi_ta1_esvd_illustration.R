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

png("../../../out/fig/main/regevEpi_ta1-inflamed_colorpalette.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
plot(1:(length(case_color_palette) + length(control_color_palette)),
     pch = 16, cex = 5, col = c(case_color_palette, control_color_palette))
graphics.off()

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

###########################################

mentioned_genes1 <- c("ASCL2", "LGR5", "LEFTY1", "CDC25C",
                      "SLC26A2", "FABP1", "CA1", "EPCAM",
                      "RBP2", "SLC26A3", "CCL20", "TNFRSF11A",
                      "BEST4", "TFF1", "MUC2", "HCK", "TRPM5",
                      "SCGN", "CHGA", "CCL7", "CCL13")
mentioned_genes2 <- c("PLA2G2A", "MUC12", "REG4", "S100P",
                      "ADIRF", "DUOX2", "MT1G", "B3GNT7",
                      "FAM3B", "ACADS", "CCL20", "FABP1", "SOD3",
                      "FFAR4", "ALDH2", "MUC4", "HOXD13", "NPSR1",
                      "SLC39A8", "URAD", "TFF1", "PRAC1", "FOSB",
                      "P13", "MIA", "MMP7", "AQP8", "SAA1",
                      "FAM213A", "SERPINB5", "CLCA1", "SLC6A14",
                      "CRABP2", "SAA2", "GABRP", "CDH3",
                      "PRSS21", "IL1RN", "FGFR2", "SAA1")
mentioned_genes1 <- unique(mentioned_genes1[mentioned_genes1 %in% colnames(eSVD_obj$dat)])
mentioned_genes2 <- unique(mentioned_genes2[mentioned_genes2 %in% colnames(eSVD_obj$dat)])

eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      alpha_max = 1e10,
                                      bool_adjust_covariates = F,
                                      bool_covariates_as_library = T)

sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Health"))
inf_de_genes <- sheet2[sheet2$ident == "TA 1","gene"]

genes <- c("CXCL3",  "RPL15")
y_vec <- eSVD_obj$fit_Second$z_mat[,"Log_UMI"]
x_vec <- log10(eSVD_obj$fit_Second$nuisance_vec)
p <- length(y_vec)
color_vec <- rep("gray", times = p)
color_vec[names(y_vec) %in% inf_de_genes] <- 2
idx <- c(which(color_vec == "gray"), which(color_vec == "2"))
png("../../../out/fig/main/regevEpi_ta1-inflamed_scatterplot_logumi-nuisance.png",
    height = 2500, width = 2000,
    units = "px", res = 500)
par(mar = c(4,4,0.1,0.1))
plot(x_vec[idx], y_vec[idx],
     col = color_vec[idx], pch = 16,
     cex = cex_vec[idx],
     xaxt = "n", yaxt = "n", bty = "n", type = "n",
     xlab = "", ylab = "")
lines(x = c(-1e4,1e4), y = rep(1,2), col = "red", lwd = 2, lty = 2)
lines(x = rep(0,2), y = c(-1e4,1e4), col = "red", lwd = 2, lty = 2)
points(x_vec[idx], y_vec[idx],
       col = color_vec[idx], pch = 16)
points(x_vec[genes], y_vec[genes],
       col = "white", pch = 16, cex = 4)
points(x_vec[genes], y_vec[genes],
       col = 2, pch = 16, cex = 3)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

den <- stats::density(y_vec)
png("../../../out/fig/main/regevEpi_ta1-inflamed_density_logumi.png",
    height = 500, width = 2500,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = range(y_vec),
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(1,2), y = c(-1e4,1e4), col = "red", lwd = 2, lty = 2)
graphics.off()


den <- stats::density(x_vec)
png("../../../out/fig/main/regevEpi_ta1-inflamed_density_nuisance.png",
    height = 500, width = 2000,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = range(x_vec),
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = "red", lwd = 2, lty = 2)
graphics.off()

########################################################

zz = cbind(round(y_vec[inf_de_genes],2), round(x_vec[inf_de_genes],2))
zz[order(zz[,2], decreasing = T),]
zz = cbind(round(y_vec[mentioned_genes1],2), round(x_vec[mentioned_genes1],2))
zz[rownames(zz) %in% inf_de_genes,,drop = F]
zz = cbind(round(y_vec[mentioned_genes2],2), round(x_vec[mentioned_genes2],2))
zz[rownames(zz) %in% inf_de_genes,,drop = F]

case_control_variable <- eSVD2:::.get_object(eSVD_obj = eSVD_obj, what_obj = "init_case_control_variable", which_fit = "param")
library_size_variables <- "Log_UMI"
library_size_variables <- unique(c(library_size_variables, setdiff(colnames(eSVD_obj$covariates), c("Intercept", case_control_variable))))
library_size_variables <- unique(c("Intercept", library_size_variables))
library_idx <- which(colnames(eSVD_obj$covariates) %in% library_size_variables)
library_mat <- exp(tcrossprod(
  eSVD_obj$covariates[,library_idx], eSVD_obj$fit_Second$z_mat[,library_idx]
))

nat_mat1 <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,-library_idx], eSVD_obj$fit_Second$z_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)

# cor_vec <- sapply(inf_de_genes, function(gene){
#   stats::cor(eSVD_obj$fit_Second$posterior_mean_mat[,gene], mean_mat_nolib[,gene])
# })
# cor_vec2 <- sapply(inf_de_genes, function(gene){
#   vec <- as.numeric(eSVD_obj$dat[,gene])/library_mat[,gene]
#   idx <- which(vec != 0)
#   stats::cor(eSVD_obj$fit_Second$posterior_mean_mat[idx,gene], vec[idx])
# })
# round(cor_vec,2)[order(cor_vec, decreasing = T)]
# round(cor_vec2,2)[order(cor_vec2, decreasing = T)]
# zz = cbind(round(y_vec[inf_de_genes],2), round(x_vec[inf_de_genes],2))
# zz[order(cor_vec, decreasing = T)[1:20],]
# zz[order(cor_vec, decreasing = F)[1:20],]

genes <- c("CXCL3",  "RPL15")

for(gene in genes){
  obs_vec <- as.numeric(eSVD_obj$dat[,gene])/library_mat[,gene]
  fit_vec <- mean_mat_nolib[,gene]
  posterior_vec <- eSVD_obj$fit_Second$posterior_mean_mat[,gene]

  png(paste0("../../../out/fig/main/regevEpi_ta1-inflamed_gene-", gene, ".png"),
      height = 1250, width = 1250,
      units = "px", res = 500)
  par(mar = c(3,3,0.5,0.5))
  set.seed(10)
  xupper <- min(c(quantile(obs_vec, probs = 0.99), quantile(fit_vec, probs = 0.99)))
  xlim <- c(-.05*xupper, xupper)
  idx <- intersect(which(obs_vec <= xupper), which(fit_vec <= xupper))
  plot(x = obs_vec[idx], y = fit_vec[idx], col = rgb(0.8, 0.8, 0.8, 0.4),
       xlim = xlim, ylim = xlim,
       xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "",
       lwd = 2, cex = 2, asp = T)
  lines(x = c(-1e4,1e4), y = c(-1e4,1e4), col = 2, lwd = 2, lty = 2)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, col = "white", cex = 3)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, col = rgb(0,0,0,0.2), cex = 2)
  axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  graphics.off()

  png(paste0("../../../out/fig/main/regevEpi_ta1-inflamed_gene-", gene, "_fulldetail.png"),
      height = 3000, width = 3000,
      units = "px", res = 300)
  set.seed(10)
  xupper <- min(c(quantile(obs_vec, probs = 0.99), quantile(fit_vec, probs = 0.99)))
  xlim <- c(-.05*xupper, xupper)
  idx <- intersect(which(obs_vec <= xupper), which(fit_vec <= xupper))
  idx2 <- which(obs_vec != 0)
  plot(x = obs_vec[idx], y = fit_vec[idx], col = "gray",
       xlim = xlim, ylim = xlim,
       xaxt = "n", yaxt = "n", bty = "n",
       main = paste0(gene, ": Nuisance: ", round(eSVD_obj$fit_Second$nuisance_vec[gene], 2), ",\nLog_UMI: ",
                     round(eSVD_obj$fit_Second$z_mat[gene,"Log_UMI"], 2), ", Cor: ",
                     round(stats::cor(obs_vec[idx2], posterior_vec[idx2]), 2)),
       lwd = 2, cex = 2, asp = T)
  lines(x = c(-1e4,1e4), y = c(-1e4,1e4), col = 2, lwd = 2, lty = 2)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, col = "white", cex = 3)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, cex = 2)
  axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  graphics.off()
}

# #################
#
# genes <- colnames(mean_mat_nolib)
# cor_vec <- sapply(genes, function(gene){
#   stats::cor(eSVD_obj$fit_Second$posterior_mean_mat[,gene], mean_mat_nolib[,gene])
# })
# cor_vec2 <- sapply(genes, function(gene){
#   vec <- as.numeric(eSVD_obj$dat[,gene])/library_mat[,gene]
#   idx <- which(vec != 0)
#   stats::cor(eSVD_obj$fit_Second$posterior_mean_mat[idx,gene], vec[idx])
# })
# round(cor_vec,2)[order(cor_vec, decreasing = T)][1:100]
# round(cor_vec2,2)[order(cor_vec2, decreasing = T)][1:100]
# zz <- intersect(order(cor_vec, decreasing = F)[1:500], order(cor_vec2, decreasing = T)[1:500])
# tmp_genes <- genes[zz]
# sapply(tmp_genes, function(gene){length(which(eSVD_obj$dat[,gene]!=0))})
# any(tmp_genes %in% inf_de_genes)
# zz <- intersect(order(cor_vec, decreasing = T)[1:200], order(cor_vec2, decreasing = F)[1:200])
# tmp_genes <- genes[zz]
# any(tmp_genes %in% inf_de_genes)
#
# genes <- c("RBCK1", "PRAC1")
# for(gene in genes){
#   obs_vec <- as.numeric(eSVD_obj$dat[,gene])/library_mat[,gene]
#   fit_vec <- mean_mat_nolib[,gene]
#   posterior_vec <- eSVD_obj$fit_Second$posterior_mean_mat[,gene]
#
#   png(paste0("../../../out/fig/main/regevEpi_ta1-inflamed_gene-", gene, ".png"),
#       height = 1500, width = 1500,
#       units = "px", res = 500)
#   par(mar = c(3,3,0.5,0.5))
#   set.seed(10)
#   xupper <- min(c(quantile(obs_vec, probs = 0.99), quantile(fit_vec, probs = 0.99)))
#   xlim <- c(-.05*xupper, xupper)
#   idx <- intersect(which(obs_vec <= xupper), which(fit_vec <= xupper))
#   plot(x = obs_vec[idx], y = fit_vec[idx], col = "gray",
#        xlim = xlim, ylim = xlim,
#        xaxt = "n", yaxt = "n", bty = "n",
#        lwd = 2, cex = 2, asp = T)
#   lines(x = c(-1e4,1e4), y = c(-1e4,1e4), col = 2, lwd = 2, lty = 2)
#   points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, col = "white", cex = 3)
#   points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, cex = 2)
#   axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
#   axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
#   graphics.off()
# }

#############################

nuisance_vec <- eSVD_obj$fit_Second$nuisance_vec
genes <- names(nuisance_vec)[which(nuisance_vec <= 0.01)]
set.seed(10); genes <- sample(genes, size = 15)

for(gene in genes){
  obs_vec <- as.numeric(eSVD_obj$dat[,gene])/library_mat[,gene]
  fit_vec <- mean_mat_nolib[,gene]
  posterior_vec <- eSVD_obj$fit_Second$posterior_mean_mat[,gene]

  png(paste0("../../../out/fig/main/regevEpi_ta1-inflamed_gene-", gene, ".png"),
      height = 1500, width = 1500,
      units = "px", res = 500)
  par(mar = c(3,3,0.5,0.5))
  set.seed(10)
  xupper <- min(c(quantile(obs_vec, probs = 0.99), quantile(fit_vec, probs = 0.99)))
  xlim <- c(-.05*xupper, xupper)
  idx <- intersect(which(obs_vec <= xupper), which(fit_vec <= xupper))
  plot(x = obs_vec[idx], y = fit_vec[idx], col = "gray",
       xlim = xlim, ylim = xlim,
       xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "",
       lwd = 2, cex = 2, asp = T)
  lines(x = c(-1e4,1e4), y = c(-1e4,1e4), col = 2, lwd = 2, lty = 2)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, col = "white", cex = 3)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, cex = 2)
  axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  graphics.off()

  png(paste0("../../../out/fig/main/regevEpi_ta1-inflamed_gene-", gene, "_fulldetail.png"),
      height = 3000, width = 3000,
      units = "px", res = 300)
  set.seed(10)
  xupper <- min(c(quantile(obs_vec, probs = 0.99), quantile(fit_vec, probs = 0.99)))
  xlim <- c(-.05*xupper, xupper)
  idx <- intersect(which(obs_vec <= xupper), which(fit_vec <= xupper))
  idx2 <- which(obs_vec != 0)
  plot(x = obs_vec[idx], y = fit_vec[idx], col = "gray",
       xlim = xlim, ylim = xlim,
       xaxt = "n", yaxt = "n", bty = "n",
       main = paste0(gene, ": Nuisance: ", round(eSVD_obj$fit_Second$nuisance_vec[gene], 2), ",\nLog_UMI: ",
                     round(eSVD_obj$fit_Second$z_mat[gene,"Log_UMI"], 2), ", Cor: ",
                     round(stats::cor(obs_vec[idx2], posterior_vec[idx2]), 2)),
       lwd = 2, cex = 2, asp = T)
  lines(x = c(-1e4,1e4), y = c(-1e4,1e4), col = 2, lwd = 2, lty = 2)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, col = "white", cex = 3)
  points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, cex = 2)
  axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  graphics.off()
}

