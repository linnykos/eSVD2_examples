rm(list=ls())
library(Seurat)
library(eSVD2)

load("../eSVD2_examples/simulation/simulation_null_4.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

seurat_obj$CC <- factor(seurat_obj$CC)
seurat_obj$Sex <- factor(seurat_obj$Sex)
set.seed(10)
seurat_obj <- Seurat::SCTransform(seurat_obj, method = "glmGamPoi",
                                  residual.features = rownames(seurat_obj[["RNA"]]),
                                  vars.to.regress = c("Age", "Sex"),
                                  verbose = T)
Seurat::Idents(seurat_obj) <- "CC"
levels(seurat_obj)

Seurat::DefaultAssay(seurat_obj) <- "SCT"
de_result <- Seurat::FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0",
                                 slot = "scale.data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)

sctransform_pval <- de_result[rownames(seurat_obj),"p_val"]
plot(sort(sctransform_pval[-c(1:10)]),
     seq(0,1,length.out = length(sctransform_pval[-c(1:10)])), asp = T)
lines(c(0,1), c(0,1), col = 2, lty = 2)
