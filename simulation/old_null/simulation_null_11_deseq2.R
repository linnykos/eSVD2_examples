rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

load("../eSVD2_examples/simulation/simulation_null_11.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(seurat_obj[["RNA"]]@counts)
categorical_var <- c("Sex", "Individual", "CC")
numerical_var <- c("Age")
metadata <- seurat_obj@meta.data[,c(categorical_var, numerical_var)]
for(var in categorical_var){
  metadata[,var] <- as.factor(metadata[,var])
}

superstring_vec <- sapply(1:nrow(metadata), function(i){
  paste0(sapply(categorical_var, function(var){
    as.character(metadata[i,var])
  }), collapse = "-")
})
unique_superstring <- unique(superstring_vec)
num_uniq <- length(unique_superstring)
mat_pseudobulk <- matrix(NA, nrow = nrow(mat), ncol = num_uniq)
rownames(mat_pseudobulk) <- rownames(mat)
colnames(mat_pseudobulk) <- unique_superstring
metadata_pseudobulk <- as.data.frame(matrix(NA, nrow = num_uniq, ncol = ncol(metadata)))
colnames(metadata_pseudobulk) <- colnames(metadata)
rownames(metadata_pseudobulk) <- unique_superstring

for(superstring in unique_superstring){
  idx <- which(superstring_vec == superstring)
  mat_pseudobulk[,superstring] <- Matrix::rowSums(mat[,idx])

  for(vr in categorical_var){
    metadata_pseudobulk[superstring, vr] <- as.character(unique(metadata[idx,vr]))
  }
  for(vr in numerical_var){
    metadata_pseudobulk[superstring, vr] <- mean(metadata[idx,vr])
  }
}

for(vr in categorical_var){
  metadata_pseudobulk[,vr] <- factor(metadata_pseudobulk[,vr])
}
for(vr in numerical_var){
  metadata_pseudobulk[,vr] <- scale(metadata_pseudobulk[,vr])
}

metadata_pseudobulk[,"CC"] <- relevel(metadata_pseudobulk[,"CC"], ref = "0")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ Sex + Age + CC)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="CC_1_vs_0")

pvalue_vec <- deseq2_res$pvalue
names(pvalue_vec) <- rownames(deseq2_res)
plot(sort(pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(pvalue_vec[-c(1:10)])), asp = T)
lines(c(0,1), c(0,1), col = 2, lty = 2)

save(deseq2_res,
     pvalue_vec,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_null_11_deseq2.RData")



