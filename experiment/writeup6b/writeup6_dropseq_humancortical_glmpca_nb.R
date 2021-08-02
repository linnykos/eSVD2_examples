rm(list=ls())
load("../../../../data/dropseq_humancortical/dropseq_humancortical_formatted.RData")

mat <- cortical[["RNA"]]@counts[Seurat::VariableFeatures(cortical),]
set.seed(10)
K <- 30
time_start <- Sys.time()
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "nb",
                             ctl = list(verbose = T), minibatch = "stochastic")
time_end <- Sys.time()
pred_mat <- glmpca:::predict.glmpca(glmpca_res)

save.image("../../../../out/writeup6b/writeup6b_dropseq_humancortical_glmpca_nb.RData")
