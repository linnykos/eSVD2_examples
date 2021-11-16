rm(list=ls())
load("../../out/writeup8d/writeup8d_sns_layer23_esvd_extended.RData")

sns <- subset(sns, features = colnames(mat))
set.seed(10)
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::SCTransform(sns, method = "glmGamPoi",
                           vars.to.regress = "percent.mt",
                           verbose = T)

