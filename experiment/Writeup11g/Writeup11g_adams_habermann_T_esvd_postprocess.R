rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/Writeup11g/Writeup11g_habermann_T_esvd.RData")
eSVD_obj$fit_Second$nuisance_vec <- pmax(eSVD_obj$fit_Second$nuisance_vec, 1e-4)
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      bool_covariates_as_library = T)
metadata <- habermann@meta.data
metadata[,"Sample_Name"] <- as.factor(metadata[,"Sample_Name"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Sample_Name",
                                           metadata = metadata,
                                           verbose = 1)
time_end5 <- Sys.time()
eSVD_obj_habermann <- eSVD_obj

load("../../../../out/Writeup11g/Writeup11g_adams_T_esvd.RData")
eSVD_obj$fit_Second$nuisance_vec <- pmax(eSVD_obj$fit_Second$nuisance_vec, 1e-4)
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      bool_covariates_as_library = T)
metadata <- adams@meta.data
metadata[,"Subject_Identity"] <- as.factor(metadata[,"Subject_Identity"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Subject_Identity",
                                           metadata = metadata,
                                           verbose = 1)
time_end5 <- Sys.time()
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
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
de_idx <- which(gene_names %in% de_genes)

xlim <- quantile(eSVD_obj_adams$teststat_vec, probs = c(0.001, 0.999))
ylim <- quantile(eSVD_obj_habermann$teststat_vec, probs = c(0.001, 0.999))
mean_x <- mean(eSVD_obj_adams$teststat_vec)
mean_y <- mean(eSVD_obj_habermann$teststat_vec)

png(paste0("../../../../out/fig/Writeup11g/Writeup11g_adams_habermann_T_teststatistic.png"),
    height = 1500, width = 3500,
    units = "px", res = 300)
par(mfrow = c(1,3))
plot(x = eSVD_obj_adams$teststat_vec,
     y = eSVD_obj_habermann$teststat_vec,
     xlab = "Adams", ylab = "Habermann", pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
     main = "All genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[cycling_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[cycling_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 3,
     main = "Cycling genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[de_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[de_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 2,
     main = "DE genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

graphics.off()


cor(eSVD_obj_adams$teststat_vec[gene_names[de_idx]],
    eSVD_obj_habermann$teststat_vec[gene_names[de_idx]])


#########################################

xlim <- quantile(eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"], probs = c(0.001, 0.999))
ylim <- quantile(eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"], probs = c(0.001, 0.999))
mean_x <- mean(eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"])
mean_y <- mean(eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"])

png(paste0("../../../../out/fig/Writeup11g/Writeup11g_adams_habermann_T_coefficient.png"),
    height = 1500, width = 3500,
    units = "px", res = 300)
par(mfrow = c(1,3))
plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
     main = "All genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[cycling_idx]],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[cycling_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 3,
     main = "Cycling genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[de_idx]],
     y = eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[de_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 2,
     main = "DE genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
lines(c(-100,100), rep(mean_y,2), col = 1, lwd = 2, lty = 2)
lines(rep(mean_x,2), c(-100,100), col = 1, lwd = 2, lty = 2)

graphics.off()

cor(eSVD_obj_adams$fit_Second$z_mat[,"Disease_Identity_IPF"][gene_names[de_idx]],
    eSVD_obj_habermann$fit_Second$z_mat[,"Diagnosis_IPF"][gene_names[de_idx]])
