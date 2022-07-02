rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/Writeup11e/Writeup11e_habermann_T_esvd9.RData")

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
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
de_idx <- which(gene_names %in% de_genes)

tab <- table(habermann$Sample_Name, habermann$Diagnosis)
indiv_cases <- rownames(tab)[which(tab[,"IPF"] != 0)]
indiv_controls <- rownames(tab)[which(tab[,"Control"] != 0)]
indiv_vec <- factor(as.character(habermann$Sample_Name))

round(apply(eSVD_obj$fit_Second$z_mat, 2, quantile), 2)

#########################################

png(paste0("../../../../out/fig/Writeup11e/Writeup11e_habermann_T_diagnostic_gene.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "sparsity",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "Diagnosis_IPF",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "sparsity",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

graphics.off()

##########################

eSVD_obj_habermann <- eSVD_obj

load("../../../../out/Writeup11e/Writeup11e_adams_T_esvd9.RData")
eSVD_obj_adams <- eSVD_obj

png(paste0("../../../../out/fig/Writeup11e/Writeup11e_adams_T_diagnostic_gene_histogram.png"),
    height = 1500, width = 2500,
    units = "px", res = 300)
eSVD2:::gene_plot(eSVD_obj_adams,
                  what_1 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))
graphics.off()


png(paste0("../../../../out/fig/Writeup11e/Writeup11e_habermann_T_diagnostic_gene_histogram.png"),
    height = 1500, width = 2500,
    units = "px", res = 300)
eSVD2:::gene_plot(eSVD_obj_habermann,
                  what_1 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))
graphics.off()

xlim <- range(eSVD_obj_adams$teststat_vec)
ylim <- range(eSVD_obj_habermann$teststat_vec)

png(paste0("../../../../out/fig/Writeup11e/Writeup11e_adams_habermann_T_teststatistic.png"),
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

plot(x = eSVD_obj_adams$teststat_vec[gene_names[cycling_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[cycling_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 3,
     main = "Cycling genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[de_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[de_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 2,
     main = "DE genes", asp = T,
     xlim = xlim, ylim = ylim)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
graphics.off()

