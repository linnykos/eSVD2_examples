rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/Writeup13/regevEpi_ta1-inflamed_esvd.RData")
regevEpi_inflamed <- regevEpi
eSVD_obj_inflamed <- eSVD_obj
load("../../../../out/Writeup13/regevEpi_ta1-noninflamed_esvd.RData")
regevEpi_noninflamed <- regevEpi
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
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
inf_de_idx <- which(gene_names %in% inf_de_genes)
other_idx <- which(gene_names %in% other_de_genes)
noninf_de_idx <- which(gene_names %in% noninf_de_genes)
hk_idx <- which(gene_names %in% hk_genes)

other_idx <- setdiff(other_idx, c(inf_de_idx, noninf_de_idx))
cycling_idx <- setdiff(cycling_idx, c(inf_de_idx, noninf_de_idx, other_idx))
hk_idx <- setdiff(hk_idx, c(inf_de_idx, noninf_de_idx, other_idx, cycling_idx))

##########################################

# make volcano plots

# first for inflamed
df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_inflamed)
teststat_vec <- eSVD_obj_inflamed$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_inflamed <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_inflamed <- locfdr::locfdr(gaussian_teststat_inflamed, plot = 0)
fdr_vec <- locfdr_res_inflamed$fdr
names(fdr_vec) <- names(gaussian_teststat_inflamed)
null_mean <- locfdr_res_inflamed$fp0["mlest", "delta"]
null_sd <- locfdr_res_inflamed$fp0["mlest", "sigma"]
logpvalue_vec_inflamed <- sapply(gaussian_teststat_inflamed, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec_inflamed <- -(logpvalue_vec_inflamed/log(10) + log10(2))
idx_inflamed <- order(logpvalue_vec_inflamed, decreasing = T)[1:length(inf_de_idx)]
length(idx_inflamed)
length(intersect(idx_inflamed, c(inf_de_idx)))
length(intersect(idx_inflamed, c(inf_de_idx, noninf_de_idx)))

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(inf_de_idx)
n <- length(gaussian_teststat_inflamed) - m
k <- length(idx_inflamed)
x <- length(intersect(idx_inflamed, c(inf_de_idx)))
fisher <- stats::dhyper(x = x, m = m, n = n, k = k, log = F)
fisher

x_vec <- log2(eSVD_obj_inflamed$case_mean) - log2(eSVD_obj_inflamed$control_mean)
xlim <- range(x_vec)
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- abs(eSVD_obj_inflamed$teststat_vec)
ylim <- range(y_vec)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

png("../../../../out/fig/Writeup13/regevEpi_ta1-inflamed_volcano.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in seq(0,10,by = 1)){
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

png("../../../../out/fig/Writeup13/regevEpi_ta1-inflamed_volcano-stats.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
plot(x = 1:10, y = 1:10, type = "n",
     main = paste0("Fisher = ", fisher, "\nTotal: ",
                   length(inf_de_idx), ", inter: ", length(intersect(idx_inflamed, inf_de_idx))))
graphics.off()

###########

# next for non-inflamed
df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_noninflamed)
teststat_vec <- eSVD_obj_noninflamed$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_noninflamed <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_noninflamed <- locfdr::locfdr(gaussian_teststat_noninflamed, plot = 0)
fdr_vec <- locfdr_res_noninflamed$fdr
names(fdr_vec) <- names(gaussian_teststat_noninflamed)
null_mean <- locfdr_res_noninflamed$fp0["mlest", "delta"]
null_sd <- locfdr_res_noninflamed$fp0["mlest", "sigma"]
logpvalue_vec_noninflamed <- sapply(gaussian_teststat_noninflamed, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec_noninflamed <- -(logpvalue_vec_noninflamed/log(10) + log10(2))
idx_noninflamed <- order(logpvalue_vec_noninflamed, decreasing = T)[1:length(noninf_de_idx)]
length(idx_noninflamed)
length(intersect(idx_noninflamed, c(noninf_de_idx)))
length(intersect(idx_noninflamed, c(inf_de_idx, noninf_de_idx)))

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(noninf_de_idx)
n <- length(gaussian_teststat_inflamed) - m
k <- length(idx_noninflamed)
x <- length(intersect(idx_noninflamed, c(noninf_de_idx)))
fisher <- stats::dhyper(x = x, m = m, n = n, k = k, log = F)
fisher

x_vec <- log2(eSVD_obj_noninflamed$case_mean) - log2(eSVD_obj_noninflamed$control_mean)
xlim <- range(x_vec) # quantile(x_vec, probs = c(0.01,0.99))
# xlim <- c(-1,1)*max(abs(xlim))
xlim <- c(-3.1,3.1)
y_vec <- abs(eSVD_obj_noninflamed$teststat_vec)
ylim <- range(y_vec)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

png("../../../../out/fig/Writeup13/regevEpi_ta1-noninflamed_volcano.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in seq(0,10,by = 1)){
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

png("../../../../out/fig/Writeup13/regevEpi_ta1-noninflamed_volcano-stats.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
plot(x = 1:10, y = 1:10, type = "n",
     main = paste0("Fisher = ", fisher, "\nTotal: ",
                   length(noninf_de_idx), ", inter: ", length(intersect(idx_noninflamed, noninf_de_idx))))
graphics.off()


##############################################

png("../../../../out/fig/Writeup13/regevEpi_ta1-agreement_de-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
y1 <- gaussian_teststat_inflamed[unique(c(inf_de_idx, noninf_de_idx))]
y2 <- gaussian_teststat_noninflamed[unique(c(inf_de_idx, noninf_de_idx))]
xbnds <- range(gaussian_teststat_inflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
ybnds <- range(gaussian_teststat_noninflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
bin <- hexbin::hexbin(y1, y2, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(y1, y2, method = "spearman"), 2)))
graphics.off()

# hk_idx2 <- setdiff(hk_idx, c(idx_inflamed, idx_noninflamed, noninf_de_genes, inf_de_genes, other_de_genes))
png("../../../../out/fig/Writeup13/regevEpi_ta1-agreement_hk-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
y1 <- gaussian_teststat_inflamed[hk_idx]
y2 <- gaussian_teststat_noninflamed[hk_idx]
xbnds <- range(gaussian_teststat_inflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
ybnds <- range(gaussian_teststat_noninflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
bin <- hexbin::hexbin(y1, y2, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(y1, y2, method = "spearman"), 2)))
graphics.off()

