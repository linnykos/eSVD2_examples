rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/Writeup11f/Writeup11f_sns_invip_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "IN-VIP"),]
de_gene_specific <- tmp[,"Gene name"]
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list[-1], function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
de_genes <- de_genes[!de_genes %in% de_gene_specific]
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(names(eSVD_obj$teststat_vec) %in% c(hk_genes, cycling_genes))
de_idx <- which(names(eSVD_obj$teststat_vec) %in% de_gene_specific)
other_idx <- which(names(eSVD_obj$teststat_vec) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj$teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

teststat_vec <- pmax(pmin(eSVD_obj$teststat_vec, 30), -30)
max_val <- max(abs(eSVD_obj$teststat_vec))
png("../../../../out/fig/Writeup11f/Writeup11f_sns_invip_esvd_teststat_histogram.png",
    height = 1200, width = 1200,
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

png("../../../../out/fig/Writeup11f/Writeup11f_sns_invip_esvd_teststat_histogram_separate.png",
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
