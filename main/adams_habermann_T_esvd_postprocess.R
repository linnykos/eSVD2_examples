rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)

load("../../../out/main/habermann_T_esvd.RData")
eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj$teststat_vec <- NULL
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      alpha_max = NULL,
                                      bool_covariates_as_library = T,
                                      library_min = NULL)
metadata <- habermann@meta.data
metadata[,"Sample_Name"] <- as.factor(metadata[,"Sample_Name"])
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Sample_Name",
                                           metadata = metadata,
                                           verbose = 1)
eSVD_obj_habermann <- eSVD_obj

load("../../../out/main/adams_T_esvd.RData")
eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj$teststat_vec <- NULL
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      alpha_max = NULL,
                                      bool_covariates_as_library = T,
                                      library_min = NULL)
metadata <- adams@meta.data
metadata[,"Subject_Identity"] <- as.factor(metadata[,"Subject_Identity"])
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Subject_Identity",
                                           metadata = metadata,
                                           verbose = 1)
eSVD_obj_adams <- eSVD_obj

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
habermann_df_genes <- setdiff(habermann_df_genes, adams_df_genes)
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

xlim <- quantile(eSVD_obj_adams$teststat_vec, probs = c(0.001, 0.999))
ylim <- quantile(eSVD_obj_habermann$teststat_vec, probs = c(0.001, 0.999))
mean_x <- mean(eSVD_obj_adams$teststat_vec)
mean_y <- mean(eSVD_obj_habermann$teststat_vec)

png(paste0("../../../out/fig/main/adams_habermann_T_teststatistic.png"),
    height = 3000, width = 4000,
    units = "px", res = 300)
par(mfrow = c(2,3))
plot(x = eSVD_obj_adams$teststat_vec,
     y = eSVD_obj_habermann$teststat_vec,
     xlab = "Adams", ylab = "Habermann", pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
     main = "All genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[adam_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[adam_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 1,
     main = "Adams genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[habermann_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[habermann_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 2,
     main = "Habermann genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[de_other_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[de_other_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 4,
     main = "DE other genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[hk_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[hk_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 3,
     main = "HK genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[cycling_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[cycling_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 5,
     main = "Cycling genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
graphics.off()

#########################################

xlim <- quantile(eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"], probs = c(0.001, 0.999))
ylim <- quantile(eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"], probs = c(0.001, 0.999))
mean_x <- mean(eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"])
mean_y <- mean(eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"])

png(paste0("../../../out/fig/main/adams_habermann_T_coefficient.png"),
    height = 3000, width = 4000,
    units = "px", res = 300)
par(mfrow = c(2,3))
plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
     main = "All genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[adam_idx]],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[adam_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 1,
     main = "Adams genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[habermann_idx]],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[habermann_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 2,
     main = "Habermann genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[de_other_idx]],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[de_other_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 4,
     main = "DE other genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[hk_idx]],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[hk_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 3,
     main = "HK genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[cycling_idx]],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[cycling_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 5,
     main = "Cycling genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

graphics.off()


##########################################

##########################################

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_adams,
                             metadata = adams@meta.data,
                             covariate_individual = "Subject_Identity")
teststat_vec <- eSVD_obj_adams$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_adams <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_adams <- locfdr::locfdr(gaussian_teststat_adams, plot = 0)
fdr_vec <- locfdr_res_adams$fdr
names(fdr_vec) <- names(gaussian_teststat_adams)
null_mean <- locfdr_res_adams$fp0["mlest", "delta"]
null_sd <- locfdr_res_adams$fp0["mlest", "sigma"]
logpvalue_vec <- sapply(gaussian_teststat_adams, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec <- -(logpvalue_vec/log10(exp(1)) + log10(2))
idx_adams <- which(fdr_vec < 0.005)
idx_adams2 <- order(logpvalue_vec, decreasing = T)[1:length(unique(c(adam_idx, habermann_idx)))]
length(idx_adams)
length(idx_adams2)

##

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_habermann,
                             metadata = habermann@meta.data,
                             covariate_individual = "Sample_Name")
teststat_vec <- eSVD_obj_habermann$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_habermann <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_habermann <- locfdr::locfdr(gaussian_teststat_habermann, plot = 0)
fdr_vec <- locfdr_res_habermann$fdr
names(fdr_vec) <- names(gaussian_teststat_habermann)
null_mean <- locfdr_res_habermann$fp0["mlest", "delta"]
null_sd <- locfdr_res_habermann$fp0["mlest", "sigma"]
logpvalue_vec <- sapply(gaussian_teststat_habermann, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec <- -(logpvalue_vec/log10(exp(1)) + log10(2))
idx_habermann <- which(fdr_vec < 0.005)
idx_habermann2 <- order(logpvalue_vec, decreasing = T)[1:length(unique(c(adam_idx, habermann_idx)))]
length(idx_habermann)
length(idx_habermann2)

###############################

# png("../../../out/fig/main/adams_habermann_T-agreement_lfdr-genes.png",
#     height = 2000, width = 2000,
#     units = "px", res = 500)
# y1 <- gaussian_teststat_adams[unique(c(idx_adams, idx_habermann))]
# # y1 <- gaussian_teststat_adams[unique(c(idx_adams2, idx_habermann2))]
# # y1 <- eSVD_obj_adams$fit_Second$z_mat[unique(c(idx_adams, idx_habermann)),"Disease_Identity_IPF"]
# y2 <- gaussian_teststat_habermann[unique(c(idx_adams, idx_habermann))]
# # y2 <- gaussian_teststat_habermann[unique(c(idx_adams2, idx_habermann2))]
# # y2 <- eSVD_obj_habermann$fit_Second$z_mat[unique(c(idx_adams, idx_habermann)),"Diagnosis_IPF"]
# bin <- hexbin::hexbin(y1, y2, xbins=20)
# my_colors <- colorRampPalette(viridis::viridis(11))
# hexbin::plot(bin, colramp=my_colors , legend=F,
#              xlab = "", ylab = "",
#              main = paste0("Cor: ", round(stats::cor(y1, y2), 2)))
# graphics.off()

png("../../../out/fig/main/adams_habermann_T-agreement_de-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
y1 <- gaussian_teststat_adams[unique(c(adam_idx, habermann_idx))]
y2 <- gaussian_teststat_habermann[unique(c(adam_idx, habermann_idx))]
xbnds <- range(gaussian_teststat_adams[c(adam_idx, habermann_idx, hk_idx)])
ybnds <- range(gaussian_teststat_habermann[c(adam_idx, habermann_idx, hk_idx)])
bin <- hexbin::hexbin(y1, y2, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(y1, y2, method = "spearman"), 2)))
graphics.off()

# hk_idx2 <- setdiff(hk_idx, c(idx_inflamed, idx_noninflamed, noninf_de_genes, inf_de_genes, other_de_genes))
png("../../../out/fig/main/adams_habermann_T-agreement_hk-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
y1 <- gaussian_teststat_adams[hk_idx]
y2 <- gaussian_teststat_habermann[hk_idx]
xbnds <- range(gaussian_teststat_adams[c(adam_idx, habermann_idx, hk_idx)])
ybnds <- range(gaussian_teststat_habermann[c(adam_idx, habermann_idx, hk_idx)])
bin <- hexbin::hexbin(y1, y2, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(y1, y2, method = "spearman"), 2)))
graphics.off()


