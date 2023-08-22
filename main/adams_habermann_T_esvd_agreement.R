rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)

load("../../../out/main/habermann_T_esvd.RData")
eSVD_obj_habermann <- eSVD_obj

load("../../../out/main/adams_T_esvd.RData")
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

###########

logpvalue_vec <- eSVD_obj_adams$pvalue_list$log10pvalue
idx_adams2 <- order(logpvalue_vec, decreasing = T)[1:length(unique(c(adam_idx, habermann_idx)))]

logpvalue_vec <- eSVD_obj_habermann$pvalue_list$log10pvalue
idx_habermann2 <- order(logpvalue_vec, decreasing = T)[1:length(unique(c(adam_idx, habermann_idx)))]

gaussian_teststat_adams <- eSVD_obj_adams$teststat_vec
gaussian_teststat_habermann <- eSVD_obj_habermann$teststat_vec

###############################

png("../../../out/fig/main/adams_habermann_T-agreement_lfdr-genes.png",
    height = 2000, width = 2000,
    units = "px", res = 500)
y1 <- gaussian_teststat_adams[unique(c(idx_adams2, idx_habermann2))]
y2 <- gaussian_teststat_habermann[unique(c(idx_adams2, idx_habermann2))]
bin <- hexbin::hexbin(y1, y2, xbins=20)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(y1, y2, method = "spearman"), 2)))
graphics.off()
stats::cor(y1, y2, method = "spearman")

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
stats::cor(y1, y2, method = "spearman")

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
stats::cor(y1, y2, method = "spearman")


