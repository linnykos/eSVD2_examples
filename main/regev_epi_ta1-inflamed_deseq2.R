rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

load("../../../out/main/regevEpi_ta1_preprocessed.RData")
# table(regevEpi$Subject, regevEpi$Sample_Health)
# table(regevEpi$Subject, regevEpi$Subject_Disease)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(regevEpi))
keep_vec[which(regevEpi$Sample_Health == "Non-inflamed")] <- 0
regevEpi$keep <- keep_vec
regevEpi <- subset(regevEpi, keep == 1)

# take only half of the healthy subjects
tab <- table(regevEpi$Subject, regevEpi$Subject_Disease)
healthy_subj <- rownames(tab[tab[,"HC"] != 0,])
set.seed(10)
split1 <- sample(healthy_subj, size = round(length(healthy_subj)/2), replace = F)
split2 <- setdiff(healthy_subj, split1)
keep_vec <- rep(1, ncol(regevEpi))
# Non-inflamed analysis uses split1, Inflamed uses split2
if(any(regevEpi$Sample_Health == "Non-inflamed")){
  keep_vec[which(regevEpi$Subject %in% split2)] <- 0
} else {
  keep_vec[which(regevEpi$Subject %in% split1)] <- 0
}
regevEpi$keep <- keep_vec
regevEpi <- subset(regevEpi, keep == 1)

regevEpi[["percent.mt"]] <- Seurat::PercentageFeatureSet(regevEpi, pattern = "^MT-")

##############################

subj_vec <- regevEpi$Sample

mat <- as.matrix(regevEpi[["RNA"]]@counts[regevEpi[["RNA"]]@var.features,])
categorical_var <- c("Subject_Disease", "Sample", "Subject_Gender", "Subject_Location", "Subject_Smoking")
numerical_var <- c("percent.mt")
metadata <- regevEpi@meta.data[,c(categorical_var, numerical_var)]
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

metadata_pseudobulk[,"Subject_Disease"] <- relevel(metadata_pseudobulk[,"Subject_Disease"], ref = "HC")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ Subject_Gender + Subject_Smoking + Subject_Location + percent.mt + Subject_Disease)

dds <- DESeq2::DESeq(dds)
nms <- DESeq2::resultsNames(dds)
deseq2_pval <- DESeq2::results(dds)$pvalue
stats::quantile(deseq2_pval, na.rm = T)

deseq2_res <- DESeq2::results(dds, name="Subject_Disease_Colitis_vs_HC")

save(regevEpi, deseq2_res,
     date_of_run, session_info,
     file = "../../../out/main/regevEpi_ta1-inflamed_deseq2.RData")

