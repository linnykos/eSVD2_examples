rm(list=ls())
load("../../../../out/writeup8c/writeup8c_sns_layer23_de_seurat.RData")
library(Seurat)

set.seed(10)
sns_de_mast <- Seurat::FindMarkers(sns,
                              assay = "SCT",
                              slot = "data",
                              ident.1 = "Control",
                              ident.2 = "ASD",
                              group.by = "diagnosis",
                              test.use = "MAST",
                              logfc.threshold = 0,
                              min.pct = 0,
                              min.cells.feature = 0,
                              verbose = T)

save(sns, sns_de, sns_de_mast,
     file = "../../../../out/writeup8c/writeup8c_sns_layer23_de_seurat_mast.RData")

