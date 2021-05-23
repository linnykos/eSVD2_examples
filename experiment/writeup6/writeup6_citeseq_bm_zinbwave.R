rm(list=ls())

## see examples of how to use zinbwave at:
# https://github.com/drisso/zinb_analysis/blob/master/real_data/espresso_covariates.Rmd
# https://github.com/willtownes/scrna2019/blob/master/algs/zinbwave_script.R

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
bm <- SeuratData::LoadData(ds = "bmcite")

Seurat::DefaultAssay(bm) <- "RNA"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm <-  Seurat::FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)

mat <- bm[["RNA"]]@counts[Seurat::VariableFeatures(bm),]
set.seed(10)
K <- 30

time_start <- Sys.time()
zinb_res <- zinbwave::zinbFit(mat, K = K, BPPARAM=BiocParallel::MulticoreParam(4),
                              verbose = T)
time_end <- Sys.time()

save.image("../../../../out/writeup6/writeup6_citeseq_bm_zinbwave.RData")


