rm(list=ls())
library(Seurat)

load("~/Dropbox/Collaboration-and-People/Kathryn Roeder - private/eSVD2/data/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- SeuratObject::LayerData(sns,
                               assay = "RNA",
                               features = Seurat::VariableFeatures(sns),
                               layer = "counts")

categorical_var <- c("diagnosis", "individual", "region", "sex", "Seqbatch", "Capbatch")
numerical_var <- c("age", "percent.mt", "RNA.Integrity.Number", "post.mortem.hours")
metadata <- sns@meta.data[,c(categorical_var, numerical_var)]

save(date_of_run, session_info,
     mat, metadata, categorical_var, numerical_var,
     file = "~/Dropbox/Collaboration-and-People/Kathryn Roeder - private/eSVD2/data/sns_layer23_mast-prepared.RData")
