rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
print("Loading in data")
dat <- anndata::read_h5ad("../../../../data/10x_mousepancreas/endocrinogenesis_day15.5.h5ad")
tmp <- Matrix::t(dat$X)
clusters <- dat$obs$clusters
rm(list = "dat")
gc()

print("Starting Seurat")
pancreas <- Seurat::CreateSeuratObject(counts = tmp)
pancreas[["celltype"]] <- clusters
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

