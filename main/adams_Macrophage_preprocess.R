rm(list=ls())
library(Seurat)
load("../../../out/main/adams_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(adams))
keep_vec[intersect(which(adams$Manuscript_Identity == "Macrophage"),
                   which(adams$Disease_Identity %in% c("Control", "IPF")))] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
adams

adams[["pca"]] <- NULL
adams[["umap"]] <- NULL

set.seed(10)
adams <- Seurat::NormalizeData(adams,
                               normalization.method = "LogNormalize", scale.factor = 10000)
adams <- Seurat::FindVariableFeatures(adams,
                                      selection.method = "vst", nfeatures = 5000)
adams <- Seurat::ScaleData(adams)

set.seed(10)
adams <- Seurat::RunPCA(adams, verbose = F)
set.seed(10)
adams <- Seurat::RunUMAP(adams, dims = 1:50)

save(adams, date_of_run, session_info,
     file = "../../../out/main/adams_Macrophage_preprocessed.RData")

