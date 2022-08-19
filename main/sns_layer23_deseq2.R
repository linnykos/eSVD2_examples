rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])
categorical_var <- c("diagnosis", "individual", "sex", "region", "Seqbatch")
numerical_var <- c("age", "percent.mt", "RNA.Integrity.Number", "post.mortem.hours")
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
                                      design = ~ sex + region + Seqbatch + age + percent.mt + RNA.Integrity.Number + post.mortem.hours + diagnosis)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="diagnosis_ASD_vs_Control")

save(sns, deseq2_res,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_deseq2.RData")

########################

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L2/3"),]
de_gene_specific <- tmp[,"Gene name"]
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list[-1], function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
de_genes <- de_genes[!de_genes %in% de_gene_specific]
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

fdr_val <- stats::p.adjust(deseq2_res$pvalue, method = "BH")
selected_genes <- rownames(deseq2_res)[which(fdr_val <= 0.05)]
length(selected_genes)
length(intersect(selected_genes, de_gene_specific))
length(intersect(selected_genes, de_genes))
length(intersect(selected_genes, sfari_genes))
length(intersect(selected_genes, hk_genes))
length(intersect(selected_genes, cycling_genes))

