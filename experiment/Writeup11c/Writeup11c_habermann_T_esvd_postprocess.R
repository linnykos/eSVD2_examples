rm(list=ls())
load("../../../../out/Writeup11c/Writeup11c_habermann_T_esvd.RData")

library(Seurat)
library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

eSVD_obj$param$init_pval_thres
quantile(eSVD_obj$fit_Second$nuisance_vec)

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
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
other_genes <- unique(c(adams_df_genes_others, habermann_df_genes_others))

hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(names(eSVD_obj$teststat_vec) %in% c(hk_genes, cycling_genes))
de_idx <- which(names(eSVD_obj$teststat_vec) %in% de_genes)
other_idx <- which(names(eSVD_obj$teststat_vec) %in% other_genes)

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj$teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

teststat_vec <- pmax(pmin(eSVD_obj$teststat_vec, 30), -30)
max_val <- max(abs(eSVD_obj$teststat_vec))
png("../../../../out/fig/Writeup11c/Writeup11c_habermann_T_esvd_teststat_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
hist(teststat_vec, breaks = break_vec,
     xlim = c(-max_val, max_val),
     main = "Histogram of test statistic",
     xlab = "Z-score", ylab = "Frequency", freq = T)
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(teststat_vec[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Other interest gene", "Housekeeping gene"),
       fill = c(2,4,3), cex = 0.6)
graphics.off()

png("../../../../out/fig/Writeup11c/Writeup11c_habermann_T_esvd_teststat_histogram_separate.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
main_vec <- c("Published DE", "Other interest\n(SFARI, DE other region)", "Housekeeping+Cell-cycle")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency", freq = T)
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()

###################

eSVD_obj_habermann <- eSVD_obj
load("../../../../out/Writeup11c/Writeup11c_adams_T_esvd.RData")
eSVD_obj_adams <- eSVD_obj

gene_list <- list(
  de_genes,
  other_genes,
  c(hk_genes, cycling_genes)
)

col_template_vec <- c(2,4,3)
idx_list <- lapply(gene_list, function(gene_vec){
  which(names(eSVD_obj_habermann$teststat_vec) %in% gene_vec)
})
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj_habermann$teststat_vec))
for(i in 1:length(idx_list)){
  col_vec[idx_list[[i]]] <- col_template_vec[i]
}

xlim <- quantile(eSVD_obj_habermann$teststat_vec, probs = c(0.01, 0.99))*1.1
ylim <- quantile(eSVD_obj_adams$teststat_vec, probs = c(0.01, 0.99))*1.1

png("../../../../out/fig/Writeup11c/Writeup11c_habermann_adams_teststatistic.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(eSVD_obj_habermann$teststat_vec, eSVD_obj_adams$teststat_vec,
     xlab = "Habermann", ylab = "Adams", pch = 16, col = col_vec,
     xlim = xlim, ylim = ylim)
graphics.off()

png("../../../../out/fig/Writeup11c/Writeup11c_habermann_adams_teststatistic_separate.png",
    height = 1200, width = 3600, res = 300, units = "px")
par(mfrow = c(1,3))
for(i in 1:length(idx_list)){
  plot(eSVD_obj_habermann$teststat_vec[idx_list[[i]]],
       eSVD_obj_adams$teststat_vec[idx_list[[i]]],
       xlab = "Habermann", ylab = "Adams", pch = 16, col = col_template_vec[i],
       xlim = xlim, ylim = ylim)
  lines(c(-15,15), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-15,15), lwd = 2, lty = 2)
}
graphics.off()

