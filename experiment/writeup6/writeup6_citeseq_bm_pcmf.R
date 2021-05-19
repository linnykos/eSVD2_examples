rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
bm <- SeuratData::LoadData(ds = "bmcite")

Seurat::DefaultAssay(bm) <- "RNA"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm <-  Seurat::FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)

mat <- bm[["RNA"]]@counts[Seurat::VariableFeatures(bm),]
mat <- t(as.matrix(mat))
set.seed(10)
K <- 30
zero_inflation <- TRUE
sparisty <- TRUE

time_start <- Sys.time()
pcmf_res <- pCMF::pCMF(mat, K = K, verbose = T, zero_inflation = zero_inflation,
                       sparsity = sparisty, ncores = 4)
time_end <- Sys.time()

save.image("../../../../out/writeup6/writeup6_citeseq_bm_pcmf.RData")


