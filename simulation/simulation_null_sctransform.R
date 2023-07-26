rm(list=ls())
library(Seurat)
load("../eSVD2_examples/simulation/simulation_null.RData")

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


pvalue_vec <- de_result$p_val
names(pvalue_vec) <- rownames(de_result)
pvalue_vec <- pvalue_vec[paste0("g", 1:length(pvalue_vec))]
plot(sort(pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(pvalue_vec[-c(1:10)])), asp = T)
lines(c(0,1), c(0,1), col = 2, lty = 2)

save(seurat_obj, de_result, pvalue_vec,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_null_sctransform.RData")
