rm(list=ls())
load("../../../../data/dropseq_humancortical/dropseq_humancortical_formatted.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- cortical[["RNA"]]@counts[Seurat::VariableFeatures(cortical),]
mat <- t(as.matrix(mat))
set.seed(10)
K <- 30
zero_inflation <- TRUE
sparisty <- TRUE

time_start <- Sys.time()
pcmf_res <- pCMF::pCMF(mat, K = K, verbose = T, zero_inflation = zero_inflation,
                       sparsity = sparisty, ncores = 4)
time_end <- Sys.time()

save.image("../../../../out/writeup6b/writeup6b_dropseq_humancortical_pcmf.RData")


