rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

for(ii in 1:3){
  for(jj in 1:3){
    print(paste0("Gene setting: ", ii, ", Individual setting: ", jj))

    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))

    set.seed(10)

    mat <- as.matrix(SeuratObject::LayerData(seurat_obj,
                                             assay = "RNA",
                                             layer = "counts"))
    categorical_var <- c("gender", "tobacco", "cc")
    numerical_var <- c("age")
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

    metadata_pseudobulk[,"cc"] <- relevel(metadata_pseudobulk[,"cc"], ref = "0")

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                          colData = metadata_pseudobulk,
                                          design = ~ gender + tobacco + age + cc)

    dds <- DESeq2::DESeq(dds)
    nms <- DESeq2::resultsNames(dds)
    deseq2_pval <- DESeq2::results(dds)$pvalue

    deseq2_res <- DESeq2::results(dds, name="cc_1_vs_0")

    date_of_run <- Sys.time()
    session_info <- devtools::session_info()
    save(deseq2_res,
         date_of_run, session_info,
         file = paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power-deseq2_geneSetting",
                       ii, "_individualSetting", jj,
                       ".RData"))
  }
}

print("Done! :)")

