rm(list=ls())
load("../../../../data/dropseq_humancortical/dropseq_humancortical_formatted.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- cortical[["RNA"]]@counts[Seurat::VariableFeatures(cortical),]
set.seed(10)
K <- 30

time_start <- Sys.time()
zinb_res <- zinbwave::zinbFit(mat, K = K, BPPARAM=BiocParallel::MulticoreParam(4),
                              verbose = T)
time_end <- Sys.time()

save.image("../../../../out/writeup6b/writeup6b_dropseq_humancortical_zinbwave.RData")


