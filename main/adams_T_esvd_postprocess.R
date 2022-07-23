rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/adams_T_esvd.RData")

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

#########################################

png(paste0("../../../out/fig/main/adams_T_diagnostic_gene.png"),
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
                  what_1 = "Disease_Identity_IPF",
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
png(paste0("../../../out/fig/main/adams_T_diagnostic_gene_histogram.png"),
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

