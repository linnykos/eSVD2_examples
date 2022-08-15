rm(list=ls())

library(Seurat)
load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sns$region <- factor(sns$region)
sns$sex <- factor(sns$sex)
sns$Seqbatch <- factor(sns$Seqbatch)
sns$Capbatch <- factor(sns$Capbatch)
set.seed(10)
sns <- Seurat::SCTransform(sns, method = "glmGamPoi",
                             residual.features = sns[["RNA"]]@var.features,
                             vars.to.regress = c("region", "sex", "Seqbatch", "Capbatch", "age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt"),
                             verbose = T)
Seurat::Idents(sns) <- "diagnosis"
levels(sns)

Seurat::DefaultAssay(sns) <- "SCT"
de_result <- Seurat::FindMarkers(sns, ident.1 = "ASD", ident.2 = "Control",
                                 slot = "scale.data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)

save(sns, de_result,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_sctransform.RData")
