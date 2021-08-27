rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

library(Seurat)

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_large.RData")

print("Forming covariates")
mat <- sns[["RNA"]]@counts[Seurat::VariableFeatures(sns), which(sns@meta.data$celltype == "L2/3")]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

metadata <- metadata[which(metadata$celltype == "L2/3"),]

save(mat, metadata, file = "../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_36501genes.RData")
