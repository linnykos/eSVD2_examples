rm(list=ls())

library(Seurat)
load("../../../out/main/adams_T_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

adams$Gender <- factor(adams$Gender)
adams$Tobacco <- factor(adams$Tobacco)
set.seed(10)
adams <- Seurat::SCTransform(adams, method = "glmGamPoi",
                             vars.to.regress = c("Gender", "percent.mt", "Tobacco", "Age"),
                             verbose = T)
Seurat::Idents(adams) <- "Disease_Identity"
levels(adams)

Seurat::DefaultAssay(adams) <- "SCT"
de_result <- Seurat::FindMarkers(adams, ident.1 = "IPF", ident.2 = "Control",
                                 slot = "scale.data",
                                 features = adams[["RNA"]]@var.features,
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)

save(adams, de_result,
     date_of_run, session_info,
     file = "../../../out/main/adams_T_sctransform.RData")
