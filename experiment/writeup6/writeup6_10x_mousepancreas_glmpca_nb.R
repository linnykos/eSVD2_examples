rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
print("Loading in data")
dat <- anndata::read_h5ad("../../../../data/10x_mousepancreas/GSE132188_adata.h5ad.h5")

pancreas <- Seurat::CreateSeuratObject(counts = Matrix::t(dat$X))
pancreas[["celltype"]] <- dat$obs$clusters
pancreas <- Seurat::NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas <-  Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)
pancreas <-  Seurat::ScaleData(pancreas)
pancreas <- Seurat::RunPCA(pancreas, features = Seurat::VariableFeatures(pancreas),
                           verbose = F)

mat <- pancreas[["RNA"]]@counts[Seurat::VariableFeatures(pancreas),]
set.seed(10)
K <- 30
time_start <- Sys.time()
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "nb",
                             ctl = list(verbose = T), minibatch = "stochastic")
time_end <- Sys.time()
pred_mat <- glmpca:::predict.glmpca(glmpca_res)

save.image("../../../../out/writeup6/writeup6_10x_mousepancreas_glmpca_nb.RData")

