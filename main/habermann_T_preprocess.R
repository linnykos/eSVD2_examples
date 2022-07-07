rm(list=ls())
library(Seurat)

load("../../../out/main/habermann_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$celltype == "T Cells")] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)

Seurat::DefaultAssay(habermann) <- "RNA"
habermann[["SCT"]] <- NULL
habermann[["umap.rna"]] <- NULL
set.seed(10)
habermann <- Seurat::SCTransform(habermann, variable.features.n = 3000)

set.seed(10); habermann <- Seurat::RunPCA(habermann, verbose = F)
set.seed(10)
habermann <- Seurat::RunUMAP(habermann, dims = 1:50)

save(habermann, date_of_run, session_info,
     file = "../../../out/main/habermann_T_preprocessed.RData")
