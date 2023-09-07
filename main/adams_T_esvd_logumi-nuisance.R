rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)

load("../../../out/main/adams_T_esvd.RData")

tab <- table(adams$Subject_Identity, adams$Disease_Identity)
case_subj <- rownames(tab)[which(tab[,"Control"] != 0)]
control_subj <- rownames(tab)[which(tab[,"Control"] == 0)]

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(length(case_subj))
control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(length(control_subj))

png("../../../out/fig/main/adams_colorpalette.png",
    height = 3000, width = 3000,
    units = "px", res = 500)
plot(1:(length(case_color_palette) + length(control_color_palette)),
     pch = 16, cex = 5, col = c(case_color_palette, control_color_palette))
graphics.off()

########################################

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]

genes <- c("RPS24",  "SMAD3")
y_vec <- eSVD_obj$fit_Second$z_mat[,"Log_UMI"]
x_vec <- -log10(eSVD_obj$fit_Second$nuisance_vec)
p <- length(y_vec)
cex_vec <- rep(1, p)
color_vec <- rep("gray", times = p)
color_vec[names(y_vec) %in% adams_df_genes] <- 2
idx <- c(which(color_vec == "gray"), which(color_vec == "2"))
png("../../../out/fig/main/adams_T_scatterplot_logumi-nuisance.png",
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
png("../../../out/fig/main/adams_T_density_logumi.png",
    height = 500, width = 2500,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = range(y_vec),
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(1,2), y = c(-1e4,1e4), col = "red", lwd = 3, lty = 2)
graphics.off()


den <- stats::density(x_vec)
png("../../../out/fig/main/adams_T_density_nuisance.png",
    height = 500, width = 2000,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = range(x_vec),
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = "red", lwd = 3, lty = 2)
graphics.off()

#################################

adams_df_genes <- intersect(adams_df_genes, colnames(eSVD_obj$dat))

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

genes <- c("RPS24",  "SMAD3")
for(gene in genes){
  obs_vec <- as.numeric(eSVD_obj$dat[,gene])/library_mat[,gene]
  fit_vec <- mean_mat_nolib[,gene]
  posterior_vec <- eSVD_obj$fit_Second$posterior_mean_mat[,gene]

  png(paste0("../../../out/fig/main/adams_T_gene-", gene, ".png"),
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

  png(paste0("../../../out/fig/main/adams_T_gene-", gene, "_fulldetail.png"),
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

