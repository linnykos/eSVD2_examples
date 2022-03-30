rm(list=ls())
load("../../../../out/Writeup11/Writeup11_adams_epithelial_preprocessed.RData")
load("../../../../out/Writeup11/Writeup11_habermann_epithelial_preprocessed.RData")

library(Seurat)
library(eSVD2)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(adams))
keep_vec[intersect(which(adams$Manuscript_Identity == "Ciliated"),
                   which(adams$Disease_Identity %in% c("Control", "IPF")))] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
tab_vec <- table(adams$Subject_Identity)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(adams))
keep_vec[which(adams$Subject_Identity %in% subj_keep)] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)

adams <- Seurat::NormalizeData(adams,
                               normalization.method = "LogNormalize", scale.factor = 10000)
adams <- Seurat::FindVariableFeatures(adams,
                                      selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "Ciliated")]

keep_vec <- rep(0, ncol(habermann))
keep_vec[intersect(which(habermann$celltype == "Ciliated"),
                   which(habermann$Diagnosis %in% c("Control", "IPF")))] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)
tab_vec <- table(habermann$Sample_Name)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$Sample_Name %in% subj_keep)] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)

Seurat::DefaultAssay(habermann) <- "RNA"
habermann <- Seurat::NormalizeData(habermann,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
habermann <- Seurat::FindVariableFeatures(habermann,
                                          selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/Ciliated_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
length(intersect(habermann_df_genes, adams_df_genes))
length(unique(c(habermann_df_genes, adams_df_genes)))
length(intersect(Seurat::VariableFeatures(habermann), Seurat::VariableFeatures(adams)))
length(unique(c(Seurat::VariableFeatures(habermann), Seurat::VariableFeatures(adams))))

hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

all_genes <- unique(c(Seurat::VariableFeatures(habermann),
                      Seurat::VariableFeatures(adams),
                      habermann_df_genes,
                      adams_df_genes,
                      hk_genes,
                      cycling_genes))
all_available_genes <- intersect(rownames(adams), rownames(habermann))
all_genes <- intersect(all_available_genes, all_genes)

adams[["RNA"]]@var.features <- all_genes
habermann[["RNA"]]@var.features <- all_genes
adams <- Seurat::ScaleData(adams)
habermann <- Seurat::ScaleData(habermann)

##############################

adams_mat <- as.matrix(Matrix::t(adams[["RNA"]]@scale.data))
habermann_mat <- as.matrix(Matrix::t(habermann[["RNA"]]@scale.data))

adams_metadata <- adams@meta.data
adams_case_individuals <- unique(adams_metadata[which(adams_metadata$Disease_Identity == "IPF"),"Subject_Identity"])
adams_control_individuals <- unique(adams_metadata[which(adams_metadata$Disease_Identity == "Control"),"Subject_Identity"])
adams_case_idx <- which(adams_metadata[,"Disease_Identity"] == "IPF")
adams_control_idx <- which(adams_metadata[,"Disease_Identity"] == "Control")

habermann_metadata <- habermann@meta.data
habermann_case_individuals <- unique(habermann_metadata[which(habermann_metadata$Diagnosis == "IPF"),"Sample_Name"])
habermann_control_individuals <- unique(habermann_metadata[which(habermann_metadata$Diagnosis == "Control"),"Sample_Name"])
habermann_case_idx <- which(habermann_metadata[,"Diagnosis"] == "IPF")
habermann_control_idx <- which(habermann_metadata[,"Diagnosis"] == "Control")

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "Ciliated")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("AT1", "AT2", "Basal", "Club", "Goblet", "Mesothelial"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/Ciliated_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("AT1_disease_vs_control_.csv", "AT2_disease_vs_control_.csv",
              "Basal_disease_vs_control_.csv", "Differentiating_Ciliated_disease_vs_control_.csv",
              "KRT5-KRT17+_disease_vs_control_.csv", "MUC5AC+_High_disease_vs_control_.csv",
              "MUC5B+_disease_vs_control_.csv", "Proliferating_Epithelial_Cells_disease_vs_control_.csv",
              "SCGB3A2+_disease_vs_control_.csv", "SCGB3A2+_SCGB1A1+_disease_vs_control_.csv",
              "Transitional_AT2_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
other_genes <- unique(c(adams_df_genes_others, habermann_df_genes_others))

# hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(colnames(adams_mat) %in% cycling_genes)
de_idx <- which(colnames(adams_mat) %in% de_genes)
other_idx <- which(colnames(adams_mat) %in% other_genes)

p <- ncol(adams_mat)
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), p)
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- unique(c(other_idx, hk_idx, de_idx))

##############################

adams_case_control <- sapply(1:p, function(j){
  mean(adams_mat[adams_case_idx,j]) - mean(adams_mat[adams_control_idx,j])
})
habermann_case_control <- sapply(1:p, function(j){
  mean(habermann_mat[habermann_case_idx,j]) - mean(habermann_mat[habermann_control_idx,j])
})

cor_val <- stats::cor(adams_case_control, habermann_case_control,
                      method = "spearman")
png("../../../../out/fig/Writeup11/Writeup11_lung_ciliated_adams_habermann_comparingMeans.png",
    width = 2500, height = 2500, res = 300, units = "px")
plot(NA, xlim = range(habermann_case_control), ylim = range(adams_case_control),
     xlab = "Habermann test statistic", ylab = "Adams test statistic",
     main = paste0("Ciliated cells (Comparing means),\nSpearman correlation: ", round(cor_val,2)))
points(x = habermann_case_control[-shuf_idx], y = adams_case_control[-shuf_idx],
       col = col_vec[-shuf_idx], pch = 16)
points(x = habermann_case_control[shuf_idx], y = adams_case_control[shuf_idx],
       col = col_vec[shuf_idx], pch = 16)
graphics.off()

#############################

adams_subject_case_mat <- t(sapply(adams_case_individuals, function(indiv){
  idx <- which(adams_metadata$Subject_Identity == indiv)
  colMeans(adams_mat[idx,])
}))
adams_subject_control_mat <- t(sapply(adams_control_individuals, function(indiv){
  idx <- which(adams_metadata$Subject_Identity == indiv)
  colMeans(adams_mat[idx,])
}))
adams_case_control <- sapply(1:p, function(j){
  mean(adams_subject_case_mat[,j]) - mean(adams_subject_control_mat[,j])
})

habermann_subject_case_mat <- t(sapply(habermann_case_individuals, function(indiv){
  idx <- which(habermann_metadata$Sample_Name == indiv)
  colMeans(habermann_mat[idx,])
}))
habermann_subject_control_mat <- t(sapply(habermann_control_individuals, function(indiv){
  idx <- which(habermann_metadata$Sample_Name == indiv)
  colMeans(habermann_mat[idx,])
}))
habermann_case_control <- sapply(1:p, function(j){
  mean(habermann_subject_case_mat[,j]) - mean(habermann_subject_control_mat[,j])
})

cor_val <- stats::cor(adams_case_control, habermann_case_control,
                      method = "spearman")
png("../../../../out/fig/Writeup11/Writeup11_lung_ciliated_adams_habermann_comparingIndivMeans.png",
    width = 2500, height = 2500, res = 300, units = "px")
plot(NA, xlim = range(habermann_case_control), ylim = range(adams_case_control),
     xlab = "Habermann test statistic", ylab = "Adams test statistic",
     main = paste0("Ciliated cells (Comparing individual's means),\nSpearman correlation: ", round(cor_val,2)))
points(x = habermann_case_control[-shuf_idx], y = adams_case_control[-shuf_idx],
       col = col_vec[-shuf_idx], pch = 16)
points(x = habermann_case_control[shuf_idx], y = adams_case_control[shuf_idx],
       col = col_vec[shuf_idx], pch = 16)
graphics.off()
