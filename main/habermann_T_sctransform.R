rm(list=ls())

library(Seurat)
load("../../../out/main/habermann_T_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

habermann$Gender <- factor(habermann$Gender)
habermann$Tobacco <- factor(habermann$Tobacco)
set.seed(10)
habermann <- Seurat::SCTransform(habermann, method = "glmGamPoi",
                             vars.to.regress = c("Gender", "percent.mt", "Tobacco", "Age"),
                             verbose = T)
Seurat::Idents(habermann) <- "Diagnosis"
levels(habermann)

Seurat::DefaultAssay(habermann) <- "SCT"
de_result <- Seurat::FindMarkers(habermann, ident.1 = "IPF", ident.2 = "Control",
                                 slot = "scale.data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)

save(habermann, de_result,
     date_of_run, session_info,
     file = "../../../out/main/habermann_T_sctransform.RData")
