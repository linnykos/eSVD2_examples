rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

file_vec <- c("../../../out/main/sns_astpp_processed.RData",
              "../../../out/main/sns_endothelial_processed.RData",
              "../../../out/main/sns_insst_processed.RData",
              "../../../out/main/sns_invip_processed.RData",
              "../../../out/main/sns_layer4_processed.RData",
              "../../../out/main/sns_layer23_processed.RData",
              "../../../out/main/sns_layer56_processed.RData",
              "../../../out/main/sns_layer56cc_processed.RData",
              "../../../out/main/sns_microglia_processed.RData",
              "../../../out/main/sns_oligo_processed.RData",
              "../../../out/main/sns_opc_processed.RData")
names(file_vec) <- c("astpp", "endothelial", "insst", "invip", "layer4", "layer23",
                     "layer56", "layer56cc", "microglia", "oligo", "opc")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

for(kk in 1:length(file_vec)){
  file <- file_vec[kk]
  print(file)
  load(file)
  celltype <- names(file_vec)[kk]

  mat <- as.matrix(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])
  categorical_var <- c("diagnosis", "individual", "sex", "region", "Seqbatch", "Capbatch")
  numerical_var <- c("age")
  metadata <- sns@meta.data[,c(categorical_var, numerical_var)]
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

  metadata_pseudobulk[,"diagnosis"] <- relevel(metadata_pseudobulk[,"diagnosis"], ref = "Control")

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                        colData = metadata_pseudobulk,
                                        design = ~ sex + region + Seqbatch + Capbatch + age + diagnosis)

  dds <- DESeq2::DESeq(dds)
  nms <- DESeq2::resultsNames(dds)
  deseq2_pval <- DESeq2::results(dds)$pvalue
  stats::quantile(deseq2_pval, na.rm = T)

  deseq2_res <- DESeq2::results(dds, name="diagnosis_ASD_vs_Control")

  save(deseq2_res,
       date_of_run, session_info,
       file = paste0("../../../out/main/sns_", celltype, "_deseq2.RData"))
}
