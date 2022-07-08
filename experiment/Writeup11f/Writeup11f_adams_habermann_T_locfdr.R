rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/habermann_T_esvd.RData")
eSVD_obj_habermann <- eSVD_obj

load("../../../out/main/adams_T_esvd.RData")
eSVD_obj_adams <- eSVD_obj

teststat_habermann <- compute_locfdr(eSVD_obj_habermann,
                                     metadata = habermann@meta.data,
                                     covariate_individual = "Sample_Name")

teststat_adams <- compute_locfdr(eSVD_obj_adams,
                                 metadata = adams@meta.data,
                                 covariate_individual = "Subject_Identity")

#########

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

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(gene_names))
col_vec[cycling_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(cycling_idx, de_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

png("../../../out/fig/Writeup11f/Writeup11f_adams_habermann_locfdr.png",
    height = 1200, width = 1200,
    units = "px", res = 300)
plot(-log10(teststat_adams$locfdr),
     -log10(teststat_habermann$locfdr),
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.3),
     xlab = "-Log10 localFDR Adams",
     ylab = "-Log10 localFDR Habermann")
points(-log10(teststat_adams$locfdr[shuf_idx]),
     -log10(teststat_habermann$locfdr[shuf_idx]),
     pch = 16, col = col_vec[shuf_idx], cex = 2)
graphics.off()
