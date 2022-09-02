rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/habermann_T_esvd.RData")
date_of_run

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("B", "Macrophage", "Macrophage Alveolar", "NK"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("Macrophages_disease_vs_control_.csv", "Monocytes_disease_vs_control_.csv",
              "B_Cells_disease_vs_control_.csv", "NK_Cells_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
adams_df_genes <- setdiff(adams_df_genes, habermann_df_genes)
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
de_genes_others <- unique(c(adams_df_genes_others, habermann_df_genes_others))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
de_genes_others <- setdiff(de_genes_others, de_genes)
cycling_genes <- setdiff(cycling_genes, c(de_genes_others, de_genes))
hk_genes <- setdiff(hk_genes, c(cycling_genes, de_genes_others, de_genes))

gene_names <- names(eSVD_obj$teststat_vec)
adam_idx <- which(gene_names %in% adams_df_genes)
habermann_idx <- which(gene_names %in% habermann_df_genes)
de_other_idx <- which(gene_names %in% de_genes_others)
cycling_idx <- which(gene_names %in% cycling_genes)
hk_idx <- which(gene_names %in% hk_genes)

#########################################

png(paste0("../../../out/fig/main/habermann_T_diagnostic_gene.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "teststat",
                  gene_list = list(gene_names[adam_idx],
                                   gene_names[habermann_idx],
                                   gene_names[de_other_idx],
                                   gene_names[hk_idx],
                                   gene_names[cycling_idx]),
                  color_palette = c(1,2,4,3,5))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "Log_UMI",
                  gene_list = list(gene_names[adam_idx],
                                   gene_names[habermann_idx],
                                   gene_names[de_other_idx],
                                   gene_names[hk_idx],
                                   gene_names[cycling_idx]),
                  color_palette = c(1,2,4,3,5))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "Diagnosis_IPF",
                  what_2 = "teststat",
                  gene_list = list(gene_names[adam_idx],
                                   gene_names[habermann_idx],
                                   gene_names[de_other_idx],
                                   gene_names[hk_idx],
                                   gene_names[cycling_idx]),
                  color_palette = c(1,2,4,3,5))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "sparsity",
                  what_2 = "teststat",
                  gene_list = list(gene_names[adam_idx],
                                   gene_names[habermann_idx],
                                   gene_names[de_other_idx],
                                   gene_names[hk_idx],
                                   gene_names[cycling_idx]),
                  color_palette = c(1,2,4,3,5))

graphics.off()


max_val <- max(abs(eSVD_obj$teststat_vec))
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
png(paste0("../../../out/fig/main/habermann_T_diagnostic_gene_histogram.png"),
    height = 2000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(2,3))
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "teststat",
                  breaks = break_vec,
                  gene_list = list(gene_names[adam_idx],
                                   gene_names[habermann_idx],
                                   gene_names[de_other_idx],
                                   gene_names[hk_idx],
                                   gene_names[cycling_idx]),
                  color_palette = c(1,2,4,3,5))

uniq_col_vec <- c(1,2,4,3,5)
idx_list <- list(adam_idx, habermann_idx, de_other_idx, hk_idx, cycling_idx)
names(idx_list) <- c("Adams", "Habermann", "Other DE", "HK", "Cell-cycle")
for(kk in 1:length(idx_list)){
  idx <- idx_list[[kk]]
  hist(eSVD_obj$teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = names(idx_list)[kk],
       xlab = "teststat", ylab = "Frequency", freq = T)
  lines(rep(median(eSVD_obj$teststat_vec),2), c(0, 1e5), lwd = 5, lty = 2, col = 2)
  rug(eSVD_obj$teststat_vec[idx], col = uniq_col_vec[kk], lwd = 2)
}
graphics.off()


##########################################

# make volcano plot
df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
teststat_vec <- eSVD_obj$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
fdr_vec <- locfdr_res$fdr
names(fdr_vec) <- names(gaussian_teststat)
null_mean <- locfdr_res$fp0["mlest", "delta"]
null_sd <- locfdr_res$fp0["mlest", "sigma"]
# null_mean <- mean(gaussian_teststat)
# tmp_range <- stats::quantile(gaussian_teststat, probs = c(0.01, 0.99))
# null_sd <- stats::sd(intersect(gaussian_teststat >= tmp_range[1], gaussian_teststat <= tmp_range[2]))
# null_sd <- stats::sd(gaussian_teststat)
logpvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec <- -(logpvalue_vec/log(10) + log10(2))
idx <- order(logpvalue_vec, decreasing = T)[1:length(habermann_idx)]
length(idx)
length(intersect(idx, c(habermann_idx)))
length(intersect(idx, c(adam_idx, habermann_idx)))

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(habermann_idx)
n <- length(gaussian_teststat) - m
k <- length(idx)
x <- length(intersect(idx, c(habermann_idx)))
fisher <- stats::dhyper(x = x, m = m, n = n, k = k, log = F)
fisher

# tab <- table(habermann$Sample_Name, habermann$Diagnosis)
# indiv_cases <- rownames(tab)[which(tab[,"IPF"] != 0)]
# indiv_controls <- rownames(tab)[which(tab[,"IPF"] == 0)]
# indiv_vec <- factor(as.character(habermann$Sample_Name))
# p <- length(gaussian_teststat)
# posterior_mat <- eSVD_obj$fit_Second$posterior_mean_mat
# x_vec <- sapply(1:p, function(j){
#   log2(mean(posterior_mat[which(indiv_vec %in% indiv_cases),j])) - log2(mean(posterior_mat[which(indiv_vec %in% indiv_controls),j]))
# })

x_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
x_vec <- pmax(pmin(x_vec, 1), -1)
xlim <- range(x_vec)
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- abs(eSVD_obj$teststat_vec)
ylim <- range(y_vec)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

png("../../../out/fig/main/habermann_T_volcano.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in seq(0,4.5,by = 0.25)){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}
points(x = x_vec, y = y_vec,
       col = rgb(0.6, 0.6, 0.6, 0.3), pch = 16)
points(x = x_vec[idx], y = y_vec[idx],
       col = orange_col, pch = 16, cex = 1.5)
points(x = x_vec[setdiff(habermann_idx, idx)], y = y_vec[setdiff(habermann_idx, idx)],
       col = purple_col, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = green_col_trans, pch = 16, cex = 1)
points(x = x_vec[intersect(idx, habermann_idx)], y = y_vec[intersect(idx, habermann_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx, habermann_idx)], y = y_vec[intersect(idx, habermann_idx)],
       col = purple_col, pch = 1, cex = 2, lwd = 2)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

png("../../../out/fig/main/habermann_T_volcano-stats.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
plot(x = 1:10, y = 1:10, type = "n",
     main = paste0("Fisher = ", fisher, "\nTotal: ",
                   length(habermann_idx), ", inter: ", length(intersect(idx, habermann_idx))))
graphics.off()

