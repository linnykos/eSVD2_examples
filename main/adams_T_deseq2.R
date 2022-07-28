rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

load("../../../out/main/adams_T_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

subj_vec <- adams$Subject_Identity

mat <- as.matrix(adams[["RNA"]]@counts[adams[["RNA"]]@var.features,])
categorical_var <- c("Disease_Identity", "Subject_Identity", "Gender", "Tobacco")
numerical_var <- c("Age", "percent.mt")
metadata <- adams@meta.data[,c(categorical_var, numerical_var)]
for(var in categorical_var){
  metadata[,var] <- as.factor(metadata[,var])
}

uniq_subj_vec <- unique(subj_vec)
num_subj <- length(uniq_subj_vec)
mat_pseudobulk <- matrix(NA, nrow = nrow(mat), ncol = num_subj)
rownames(mat_pseudobulk) <- rownames(mat)
colnames(mat_pseudobulk) <- uniq_subj_vec
metadata_pseudobulk <- as.data.frame(matrix(NA, nrow = num_subj, ncol = ncol(metadata)))
colnames(metadata_pseudobulk) <- colnames(metadata)
rownames(metadata_pseudobulk) <- uniq_subj_vec

for(subj in uniq_subj_vec){
  idx <- which(subj_vec == subj)
  mat_pseudobulk[,subj] <- Matrix::rowSums(mat[,idx])

  for(vr in categorical_var){
    metadata_pseudobulk[subj, vr] <- as.character(unique(metadata[idx,vr]))
  }
  for(vr in numerical_var){
    metadata_pseudobulk[subj, vr] <- mean(metadata[idx,vr])
  }
}

for(vr in categorical_var){
  metadata_pseudobulk[,vr] <- factor(metadata_pseudobulk[,vr])
}
for(vr in numerical_var){
  metadata_pseudobulk[,vr] <- scale(metadata_pseudobulk[,vr])
}

metadata_pseudobulk[,"Disease_Identity"] <- relevel(metadata_pseudobulk[,"Disease_Identity"], ref = "Control")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ Gender + Tobacco + Age + percent.mt + Disease_Identity)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="Disease_Identity_IPF_vs_Control")

save(adams, deseq2_res,
     date_of_run, session_info,
     file = "../../../out/main/adams_T_deseq2.RData")
