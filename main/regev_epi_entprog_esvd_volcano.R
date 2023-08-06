rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_entprog-inflamed_esvd.RData")
eSVD_obj_inflamed <- eSVD_obj
load("../../../out/main/regevEpi_entprog-noninflamed_esvd.RData")
eSVD_obj_noninflamed <- eSVD_obj

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Non-Inflamed vs. He"))
sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Health"))
sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Non-In"))
noninf_de_genes <- sheet1[sheet1$ident == "Enterocyte Progenitors","gene"]
inf_de_genes <- sheet2[sheet2$ident == "Enterocyte Progenitors","gene"]
other_de_genes <- sheet3[sheet3$ident == "Enterocyte Progenitors","gene"]
other_de_genes <- setdiff(other_de_genes, c(noninf_de_genes, inf_de_genes))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
cycling_genes <- setdiff(cycling_genes, c(other_de_genes, noninf_de_genes, inf_de_genes))
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
inf_de_idx <- which(gene_names %in% inf_de_genes)
other_idx <- which(gene_names %in% other_de_genes)
noninf_de_idx <- which(gene_names %in% noninf_de_genes)
hk_idx <- which(gene_names %in% hk_genes)

other_idx <- setdiff(other_idx, c(inf_de_idx, noninf_de_idx))
cycling_idx <- setdiff(cycling_idx, c(inf_de_idx, noninf_de_idx, other_idx))
hk_idx <- setdiff(hk_idx, c(inf_de_idx, noninf_de_idx, other_idx, cycling_idx))

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

#################################

obj <- eSVD_obj_inflamed
x_vec <- log2(obj$case_mean) - log2(obj$control_mean)
xlim <- range(x_vec)
xlim <- c(-max(abs(xlim)), max(abs(xlim)))
y_vec <- obj$pvalue_list$log10pvalue
ylim <- range(y_vec)
idx_inflamed <- order(y_vec, decreasing = T)[1:length(inf_de_idx)]

png("../../../out/fig/main/regevEpi_entprog-inflamed_volcano.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in pretty(seq(0,max(ylim),by = 10))){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}
points(x = x_vec, y = y_vec,
       col = rgb(0.6, 0.6, 0.6, 0.3), pch = 16)
points(x = x_vec[idx_inflamed], y = y_vec[idx_inflamed],
       col = orange_col, pch = 16, cex = 1.5)
points(x = x_vec[setdiff(inf_de_idx, idx_inflamed)], y = y_vec[setdiff(inf_de_idx, idx_inflamed)],
       col = purple_col, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = green_col_trans, pch = 16, cex = 1)
points(x = x_vec[intersect(idx_inflamed, inf_de_idx)], y = y_vec[intersect(idx_inflamed, inf_de_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx_inflamed, inf_de_idx)], y = y_vec[intersect(idx_inflamed, inf_de_idx)],
       col = purple_col, pch = 1, cex = 2, lwd = 2)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

###########

obj <- eSVD_obj_noninflamed
x_vec <- log2(obj$case_mean) - log2(obj$control_mean)
xlim <- range(x_vec)
xlim <- c(-max(abs(xlim)), max(abs(xlim)))
y_vec <- obj$pvalue_list$log10pvalue
ylim <- range(y_vec)
idx_noninflamed <- order(y_vec, decreasing = T)[1:length(noninf_de_idx)]

png("../../../out/fig/main/regevEpi_entprog-noninflamed_volcano.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in pretty(seq(0,max(ylim),by = 10))){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}
points(x = x_vec, y = y_vec,
       col = rgb(0.6, 0.6, 0.6, 0.3), pch = 16)
points(x = x_vec[idx_noninflamed], y = y_vec[idx_noninflamed],
       col = orange_col, pch = 16, cex = 1.5)
points(x = x_vec[setdiff(noninf_de_idx, idx_noninflamed)], y = y_vec[setdiff(noninf_de_idx, idx_noninflamed)],
       col = purple_col, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = green_col_trans, pch = 16, cex = 1)
points(x = x_vec[intersect(idx_noninflamed, noninf_de_idx)], y = y_vec[intersect(idx_noninflamed, noninf_de_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx_noninflamed, noninf_de_idx)], y = y_vec[intersect(idx_noninflamed, noninf_de_idx)],
       col = purple_col, pch = 1, cex = 2, lwd = 2)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()
