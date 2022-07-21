rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/Writeup11h/Writeup11h_regev_esvd_ta1_inflamed.RData")
eSVD_obj_inflamed <- eSVD_obj
load("../../../../out/Writeup11h/Writeup11h_regev_esvd_ta1_noninflamed.RData")
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
noninf_de_genes <- sheet1[sheet1$ident == "TA 1","gene"]
inf_de_genes <- sheet2[sheet2$ident == "TA 1","gene"]
other_de_genes <- sheet3[sheet3$ident == "TA 1","gene"]
other_de_genes <- setdiff(other_de_genes, c(noninf_de_genes, inf_de_genes))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
cycling_genes <- setdiff(cycling_genes, c(other_de_genes, noninf_de_genes, inf_de_genes))

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
inf_de_idx <- which(gene_names %in% inf_de_genes)
other_idx <- which(gene_names %in% other_de_genes)
noninf_de_idx <- which(gene_names %in% noninf_de_genes)

#########################################

png(paste0("../../../../out/fig/Writeup11h/Writeup11h_regev_ta1_inflamed_diagnostic_gene.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj_inflamed,
                  what_1 = "nuisance",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[inf_de_idx],
                                   gene_names[unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx))]),
                  color_palette = c(3,2,4))

eSVD2:::gene_plot(eSVD_obj_inflamed,
                  what_1 = "nuisance",
                  what_2 = "sparsity",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[inf_de_idx],
                                   gene_names[unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx))]),
                  color_palette = c(3,2,4))

eSVD2:::gene_plot(eSVD_obj_inflamed,
                  what_1 = "Subject_Disease_Colitis",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[inf_de_idx],
                                   gene_names[unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx))]),
                  color_palette = c(3,2,4))

eSVD2:::gene_plot(eSVD_obj_inflamed,
                  what_1 = "sparsity",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[inf_de_idx],
                                   gene_names[unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx))]),
                  color_palette = c(3,2,4))
graphics.off()


png(paste0("../../../../out/fig/Writeup11h/Writeup11h_regev_ta1_noninflamed_diagnostic_gene.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj_noninflamed,
                  what_1 = "nuisance",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[noninf_de_idx],
                                   gene_names[unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx))]),
                  color_palette = c(3,2,4))

eSVD2:::gene_plot(eSVD_obj_noninflamed,
                  what_1 = "nuisance",
                  what_2 = "sparsity",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[noninf_de_idx],
                                   gene_names[unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx))]),
                  color_palette = c(3,2,4))

eSVD2:::gene_plot(eSVD_obj_noninflamed,
                  what_1 = "Subject_Disease_Colitis",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[noninf_de_idx],
                                   gene_names[unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx))]),
                  color_palette = c(3,2,4))

eSVD2:::gene_plot(eSVD_obj_noninflamed,
                  what_1 = "sparsity",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[noninf_de_idx],
                                   gene_names[unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx))]),
                  color_palette = c(3,2,4))

graphics.off()

#####################################3

max_val <- max(abs(eSVD_obj_inflamed$teststat_vec))
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
png(paste0("../../../../out/fig/Writeup11h/Writeup11h_regev_ta1_inflamed_diagnostic_gene_histogram.png"),
    height = 2000, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2))
eSVD2:::gene_plot(eSVD_obj_inflamed,
                  what_1 = "teststat",
                  breaks = break_vec,
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[inf_de_idx],
                                   gene_names[unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx))]),
                  color_palette = c(3,2,4))

uniq_col_vec <- c(3,2,4)
idx_list <- list(cycling_idx, inf_de_idx, unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx)))
names(idx_list) <- c("Cell-cycle", "Inflamed", "Other")
for(kk in 1:length(idx_list)){
  idx <- idx_list[[kk]]
  hist(eSVD_obj_inflamed$teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = names(idx_list)[kk],
       xlab = "teststat", ylab = "Frequency", freq = T)
  lines(rep(median(eSVD_obj_inflamed$teststat_vec),2), c(0, 1e5), lwd = 5, lty = 2, col = 2)
  rug(eSVD_obj_inflamed$teststat_vec[idx], col = uniq_col_vec[kk], lwd = 2)
}
graphics.off()

max_val <- max(abs(eSVD_obj_noninflamed$teststat_vec))
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
png(paste0("../../../../out/fig/Writeup11h/Writeup11h_regev_ta1_noninflamed_diagnostic_gene_histogram.png"),
    height = 2000, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2))
eSVD2:::gene_plot(eSVD_obj_noninflamed,
                  what_1 = "teststat",
                  breaks = break_vec,
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[noninf_de_idx],
                                   gene_names[unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx))]),
                  color_palette = c(3,2,4))

uniq_col_vec <- c(3,2,4)
idx_list <- list(cycling_idx, noninf_de_idx, unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx)))
names(idx_list) <- c("Cell-cycle", "Non-inflamed", "Other")
for(kk in 1:length(idx_list)){
  idx <- idx_list[[kk]]
  hist(eSVD_obj_noninflamed$teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = names(idx_list)[kk],
       xlab = "teststat", ylab = "Frequency", freq = T)
  lines(rep(median(eSVD_obj_noninflamed$teststat_vec),2), c(0, 1e5), lwd = 5, lty = 2, col = 2)
  rug(eSVD_obj_noninflamed$teststat_vec[idx], col = uniq_col_vec[kk], lwd = 2)
}
graphics.off()

########################

xlim <- quantile(eSVD_obj_noninflamed$teststat_vec, probs = c(0.001, 0.999))
ylim <- quantile(eSVD_obj_inflamed$teststat_vec, probs = c(0.001, 0.999))
mean_x <- mean(eSVD_obj_noninflamed$teststat_vec)
mean_y <- mean(eSVD_obj_inflamed$teststat_vec)

png(paste0("../../../../out/fig/Writeup11h/Writeup11h_regev_ta1_agreement_teststatistic.png"),
    height = 1500, width = 3500,
    units = "px", res = 300)
par(mfrow = c(1,3))
plot(x = eSVD_obj_noninflamed$teststat_vec,
     y = eSVD_obj_inflamed$teststat_vec,
     xlab = "Noninflamed", ylab = "Inflamed", pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
     main = "TA 1", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_noninflamed$teststat_vec[gene_names[cycling_idx]],
     y = eSVD_obj_inflamed$teststat_vec[gene_names[cycling_idx]],
     xlab = "Noninflamed", ylab = "Inflamed", pch = 16, col = 3,
     main = "Cycling genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_noninflamed$teststat_vec[gene_names[unique(c(noninf_de_idx, inf_de_idx, other_idx))]],
     y = eSVD_obj_inflamed$teststat_vec[gene_names[unique(c(noninf_de_idx, inf_de_idx, other_idx))]],
     xlab = "Noninflamed", ylab = "Inflamed", pch = 16, col = 2,
     main = "DE genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

graphics.off()


