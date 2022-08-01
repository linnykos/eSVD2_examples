rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)

load("../../../out/main/habermann_T_deseq2.RData")
habermann_deseq2 <- deseq2_res

load("../../../out/main/adams_T_deseq2.RData")
adams_deseq2 <- deseq2_res

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

#############################

adam_idx <- sapply(c(adams_df_genes, habermann_df_genes), function(x){
  zz <- which(rownames(adams_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
adams_de <- adams_deseq2[adam_idx,"log2FoldChange"]
habermann_idx <- sapply(c(adams_df_genes, habermann_df_genes), function(x){
  zz <- which(rownames(habermann_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
habermann_de <- habermann_deseq2[habermann_idx,"log2FoldChange"]
rm_idx <- unique(c(which(is.na(adams_de)), which(is.na(habermann_de))))
adams_de <- adams_de[-rm_idx]; habermann_de <- habermann_de[-rm_idx]

adam_idx <- sapply(hk_genes, function(x){
  zz <- which(rownames(adams_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
adams_hk <- adams_deseq2[adam_idx,"log2FoldChange"]
habermann_idx <- sapply(hk_genes, function(x){
  zz <- which(rownames(habermann_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
habermann_hk <- habermann_deseq2[habermann_idx,"log2FoldChange"]
rm_idx <- unique(c(which(is.na(adams_hk)), which(is.na(habermann_hk))))
adams_hk <- adams_hk[-rm_idx]; habermann_hk <- habermann_hk[-rm_idx]

png("../../../out/fig/main/adams_habermann_T-agreement_deseq2_de-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
xbnds <- range(c(adams_de, adams_hk))
ybnds <- range(c(habermann_de, habermann_hk))
bin <- hexbin::hexbin(adams_de, habermann_de, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(adams_de, habermann_de, method = "spearman"), 2)))
graphics.off()

png("../../../out/fig/main/adams_habermann_T-agreement_deseq2_hk-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
xbnds <- range(c(adams_de, adams_hk))
ybnds <- range(c(habermann_de, habermann_hk))
bin <- hexbin::hexbin(adams_hk, habermann_hk, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(adams_hk, habermann_hk, method = "spearman"), 2)))
graphics.off()

