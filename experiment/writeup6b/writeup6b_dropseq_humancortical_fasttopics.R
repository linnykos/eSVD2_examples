rm(list=ls())
load("../../../../data/dropseq_humancortical/dropseq_humancortical_formatted.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- cortical[["RNA"]]@counts[Seurat::VariableFeatures(cortical),]
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
time_start <- Sys.time()
topic_res <- fastTopics::fit_topic_model(mat, k = K)
time_end <- Sys.time()

save.image("../../../../out/writeup6b/writeup6b_dropseq_humancortical_fasttopics.RData")




