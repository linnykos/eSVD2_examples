rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta1-inflamed_esvd.RData")
regevEpi_inflamed <- regevEpi
eSVD_obj_inflamed <- eSVD_obj
load("../../../out/main/regevEpi_ta1-noninflamed_esvd.RData")
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

###################

max_val <- max(abs(eSVD_obj_inflamed$teststat_vec))
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
png(paste0("../../../out/fig/main/regev_ta1_inflamed_gene_histogram.png"),
    height = 2000, width = 3500,
    units = "px", res = 300)
par(mfrow = c(2,3))
eSVD2:::gene_plot(eSVD_obj_inflamed,
                  what_1 = "teststat",
                  breaks = break_vec,
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[hk_idx],
                                   gene_names[inf_de_idx],
                                   gene_names[unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx))]),
                  color_palette = c(5,3,2,4))

uniq_col_vec <- c(5,3,2,4)
idx_list <- list(cycling_idx, hk_idx, inf_de_idx, unique(setdiff(c(noninf_de_idx, other_idx), inf_de_idx)))
names(idx_list) <- c("Cell-cycle", "Housekeeping", "Inflamed", "Other")
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
png(paste0("../../../out/fig/main/regev_ta1_noninflamed_gene_histogram.png"),
    height = 2000, width = 3500,
    units = "px", res = 300)
par(mfrow = c(2,3))
eSVD2:::gene_plot(eSVD_obj_noninflamed,
                  what_1 = "teststat",
                  breaks = break_vec,
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[hk_idx],
                                   gene_names[noninf_de_idx],
                                   gene_names[unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx))]),
                  color_palette = c(5,3,2,4))

uniq_col_vec <- c(5,3,2,4)
idx_list <- list(cycling_idx, hk_idx, noninf_de_idx, unique(setdiff(c(inf_de_idx, other_idx), noninf_de_idx)))
names(idx_list) <- c("Cell-cycle", "Housekeeping",  "Non-inflamed", "Other")
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

png(paste0("../../../out/fig/main/regev_ta1_agreement.png"),
    height = 1200, width = 4000,
    units = "px", res = 300)
par(mfrow = c(1,4))
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
     xlab = "Noninflamed", ylab = "Inflamed", pch = 16, col = 5,
     main = "Cycling genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_noninflamed$teststat_vec[gene_names[hk_idx]],
     y = eSVD_obj_inflamed$teststat_vec[gene_names[hk_idx]],
     xlab = "Noninflamed", ylab = "Inflamed", pch = 16, col = 4,
     main = "Housekeeping genes", asp = T,
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


##########################################

# make volcano plots

# first for inflamed
eSVD_obj_inflamed$fit_Second$posterior_mean_mat <- NULL
eSVD_obj_inflamed$fit_Second$posterior_var_mat <- NULL
eSVD_obj_inflamed$teststat_vec <- NULL
eSVD_obj_inflamed <- eSVD2:::compute_posterior(input_obj = eSVD_obj_inflamed,
                                               bool_adjust_covariates = F,
                                               alpha_max = NULL,
                                               bool_covariates_as_library = T,
                                               library_min = 1e-4)
metadata <- regevEpi_inflamed@meta.data
metadata[,"Sample"] <- as.factor(metadata[,"Sample"])
eSVD_obj_inflamed <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj_inflamed,
                                                    covariate_individual = "Sample",
                                                    metadata = metadata,
                                                    verbose = 1)

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_inflamed,
                             metadata = regevEpi_inflamed@meta.data,
                             covariate_individual = "Sample")
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

x_vec <- log2(eSVD_obj_inflamed$case_mean) - log2(eSVD_obj_inflamed$control_mean)
xlim <- range(x_vec)
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- abs(eSVD_obj_inflamed$teststat_vec)
ylim <- range(y_vec)

png("../../../out/fig/main/regevEpi_ta1-inflamed_volcano.png",
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
       col = 2, pch = 16, cex = 1.5)
points(x = x_vec[setdiff(inf_de_idx, idx_inflamed)], y = y_vec[setdiff(inf_de_idx, idx_inflamed)],
       col = 3, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = rgb(34, 151, 230, 255*.35, maxColorValue = 255), pch = 16, cex = 1)
points(x = x_vec[intersect(idx_inflamed, inf_de_idx)], y = y_vec[intersect(idx_inflamed, inf_de_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx_inflamed, inf_de_idx)], y = y_vec[intersect(idx_inflamed, inf_de_idx)],
       col = 3, pch = 1, cex = 2, lwd = 2)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

png("../../../out/fig/main/regevEpi_ta1-inflamed_volcano-stats.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
plot(x = 1:10, y = 1:10, type = "n",
     main = paste0("Fisher = ", fisher, "\nTotal: ",
                   length(inf_de_idx), ", inter: ", length(intersect(idx_inflamed, inf_de_idx))))
graphics.off()

###########

# next for non-inflamed
eSVD_obj_noninflamed$fit_Second$posterior_mean_mat <- NULL
eSVD_obj_noninflamed$fit_Second$posterior_var_mat <- NULL
eSVD_obj_noninflamed$teststat_vec <- NULL
eSVD_obj_noninflamed <- eSVD2:::compute_posterior(input_obj = eSVD_obj_noninflamed,
                                                  bool_adjust_covariates = F,
                                                  alpha_max = NULL,
                                                  bool_covariates_as_library = T,
                                                  library_min = 1e-4)
metadata <- regevEpi@meta.data
metadata[,"Sample"] <- as.factor(metadata[,"Sample"])
eSVD_obj_noninflamed <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj_noninflamed,
                                                       covariate_individual = "Sample",
                                                       metadata = metadata,
                                                       verbose = 1)

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_noninflamed,
                             metadata = regevEpi_noninflamed@meta.data,
                             covariate_individual = "Sample")
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

x_vec <- log2(eSVD_obj_noninflamed$case_mean) - log2(eSVD_obj_noninflamed$control_mean)
xlim <- range(x_vec) # quantile(x_vec, probs = c(0.01,0.99))
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- abs(eSVD_obj_noninflamed$teststat_vec)
ylim <- range(y_vec)

png("../../../out/fig/main/regevEpi_ta1-noninflamed_volcano.png",
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
       col = 2, pch = 16, cex = 1.5)
points(x = x_vec[setdiff(noninf_de_idx, idx_noninflamed)], y = y_vec[setdiff(noninf_de_idx, idx_noninflamed)],
       col = 3, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = rgb(34, 151, 230, 255*.35, maxColorValue = 255), pch = 16, cex = 1)
points(x = x_vec[intersect(idx_noninflamed, noninf_de_idx)], y = y_vec[intersect(idx_noninflamed, noninf_de_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx_noninflamed, noninf_de_idx)], y = y_vec[intersect(idx_noninflamed, noninf_de_idx)],
       col = 3, pch = 1, cex = 2, lwd = 2)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

png("../../../out/fig/main/regevEpi_ta1-noninflamed_volcano-stats.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
plot(x = 1:10, y = 1:10, type = "n",
     main = paste0("Fisher = ", fisher, "\nTotal: ",
                   length(noninf_de_idx), ", inter: ", length(intersect(idx_noninflamed, noninf_de_idx))))
graphics.off()


###############################

# png("../../../out/fig/main/regevEpi_ta1-agreement_lfdr-genes.png",
#     height = 2000, width = 2000,
#     units = "px", res = 500)
# par(mar = c(3,3,0.1,0.1))
# y1 <- gaussian_teststat_inflamed[unique(c(idx_inflamed, idx_noninflamed))]
# names(y1) <- NULL
# y2 <- gaussian_teststat_noninflamed[unique(c(idx_inflamed, idx_noninflamed))]
# names(y2) <- NULL
# bin <- hexbin::hexbin(y1, y2, xbins=20)
# my_colors <- colorRampPalette(viridis::viridis(11))
# hexbin::plot(bin, main="", colramp=my_colors , legend=F,
#               xlab = "", ylab = "")
# graphics.off()

png("../../../out/fig/main/regevEpi_ta1-agreement_de-genes.png",
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
png("../../../out/fig/main/regevEpi_ta1-agreement_hk-genes.png",
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

##############################################

tab <- table(regevEpi_inflamed$Sample, regevEpi_inflamed$Subject_Disease)
sample_names <- rownames(tab)
sample_names <- sample_names[which(rowSums(tab) >= 150)]
subject_names <- sapply(sample_names, function(x){
  strsplit(x, split = "\\.")[[1]][1]
})
valid_subject_names <- names(table(subject_names)[which(table(subject_names) == 2)])
column_pairs <- lapply(valid_subject_names, function(subject_name){
  grep(paste0("Sample_", subject_name, "\\."), colnames(eSVD_obj_inflamed$covariates))
})
names(column_pairs) <- valid_subject_names
column_pairs <- column_pairs[which(sapply(column_pairs, length) == 2)]
mat_inflamed <- do.call(rbind, lapply(column_pairs, function(vec){
  zz <- cbind(eSVD_obj_inflamed$fit_Second$z_mat[,vec])
  idx1 <- intersect(which(zz[,1] >= quantile(zz[,1], probs = 0.01)),
                    which(zz[,1] <= quantile(zz[,1], probs = 0.99)))
  idx2 <- intersect(which(zz[,2] >= quantile(zz[,2], probs = 0.01)),
                    which(zz[,2] <= quantile(zz[,2], probs = 0.99)))
  zz <- zz[intersect(idx1, idx2),]
  scale(zz)
  print(cor(zz[,1], zz[,2]))
  zz
}))
cor(mat_inflamed[,1], mat_inflamed[,2])

png("../../../out/fig/main/regevEpi_ta1-sample-agreement.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
par(mar = c(3,3,0.5,0.5))
plot(mat_inflamed[,1], mat_inflamed[,2],
     type = "n", xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n", asp = T)
lines(c(-1e4,1e4), c(-1e4,1e4), col = 2, lwd = 2, lty = 2)
points(mat_inflamed[,1], mat_inflamed[,2], col = "gray", pch = 16)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
