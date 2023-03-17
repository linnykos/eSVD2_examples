rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)

load("../../../../out/Writeup13/habermann_T_esvd.RData")
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
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
de_genes_others <- setdiff(de_genes_others, de_genes)
cycling_genes <- setdiff(cycling_genes, c(de_genes_others, de_genes))
hk_genes <- setdiff(hk_genes, c(cycling_genes, de_genes_others, de_genes))

gene_names <- names(eSVD_obj$teststat_vec)
adam_idx <- which(gene_names %in% adams_df_genes)
habermann_idx <- which(gene_names %in% habermann_df_genes)
de_other_idx <- which(gene_names %in% de_genes_others)
cycling_idx <- which(gene_names %in% cycling_genes)
hk_idx <- which(gene_names %in% hk_genes)

##########################################

# make volcano plot
df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
teststat_vec <- eSVD_obj$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

.compute_multtest <- function(gaussian_teststat){
  tryCatch(expr = {
    locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
    null_mean <- locfdr_res$fp0["mlest", "delta"]
    null_sd <- locfdr_res$fp0["mlest", "sigma"]
    return(list(null_mean = null_mean, null_sd = null_sd))
  }, error = function(e){
    null_mean <- mean(gaussian_teststat)
    null_sd <- sd(gaussian_teststat)
    return(list(null_mean = null_mean, null_sd = null_sd))
  })
}
tmp <- .compute_multtest(gaussian_teststat)
null_mean <- tmp$null_mean; null_sd <- tmp$null_sd
logpvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec <- -(logpvalue_vec/log(10) + log10(2))
idx <- order(logpvalue_vec, decreasing = T)[1:length(adams_df_genes)]
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

x_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control_mean)
xlim <- range(x_vec)
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- abs(eSVD_obj$teststat_vec)
ylim <- range(y_vec)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

png("../../../../out/fig/Writeup13/habermann_T_esvd_volcano.png",
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

png("../../../../out/fig/Writeup13/habermann_T_esvd_volcano-stats.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
plot(x = 1:10, y = 1:10, type = "n",
     main = paste0("Fisher = ", fisher, "\nTotal: ",
                   length(habermann_idx), ", inter: ", length(intersect(idx, habermann_idx))))
graphics.off()

