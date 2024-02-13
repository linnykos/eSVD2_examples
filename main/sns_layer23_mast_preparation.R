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

save(date_of_run, session_info,
     mat,
     file = "~/Dropbox/Collaboration-and-People/Kathryn Roeder - private/eSVD2/data/sns_layer23_mast-prepared.RData")
