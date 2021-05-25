rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
print("Loading in data")
dat <- anndata::read_h5ad("../../../../data/dropseq_mouselung/lung_regeneration_after_bleo")
tmp <- Matrix::t(dat$X)
colnames(tmp) <- gsub(pattern = "_", replacement = "-", x = colnames(tmp))

print("Starting Seurat")
lung <- Seurat::CreateSeuratObject(counts = tmp)
lung[["celltype"]] <- dat$obs$clusters
lung <- Seurat::NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung <-  Seurat::FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)
lung <-  Seurat::ScaleData(lung)
lung <- Seurat::RunPCA(lung, features = Seurat::VariableFeatures(lung),
                       verbose = F)

mat <- lung[["RNA"]]@counts[Seurat::VariableFeatures(lung),]
set.seed(10)
K <- 30
time_start <- Sys.time()
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "nb",
                             ctl = list(verbose = T), minibatch = "stochastic")
time_end <- Sys.time()
pred_mat <- glmpca:::predict.glmpca(glmpca_res)

save.image("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")
