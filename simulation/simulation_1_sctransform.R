rm(list=ls())
library(Seurat)
load("../eSVD2_examples/simulation/simulation_1.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

seurat_obj$cc <- factor(seurat_obj$cc)
seurat_obj$gender <- factor(seurat_obj$gender)
seurat_obj$tobacco <- factor(seurat_obj$tobacco)
set.seed(10)
seurat_obj <- Seurat::SCTransform(seurat_obj, method = "glmGamPoi",
                           residual.features = rownames(seurat_obj[["RNA"]]),
                           vars.to.regress = c("age", "gender", "tobacco"),
                           verbose = T)
Seurat::Idents(seurat_obj) <- "cc"
levels(seurat_obj)

Seurat::DefaultAssay(seurat_obj) <- "SCT"
de_result <- Seurat::FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0",
                                 slot = "scale.data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)

save(seurat_obj, de_result,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_1_sctransform.RData")
