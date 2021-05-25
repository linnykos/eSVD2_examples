rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
load("../../../../data/10x_mouseretinal/10x_mouseretinal_formatted.RData")

rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

retinal <- Seurat::CreateSeuratObject(counts = Matrix::t(dat),
                                      meta.data = metadata, min.cells = 5)
largest_batch <- names(which.max(table(retinal@meta.data$BatchID)))
cell_keep <- rownames(retinal@meta.data[retinal@meta.data$BatchID == largest_batch,])
retinal <- retinal[,cell_keep]
retinal <- Seurat::NormalizeData(retinal, normalization.method = "LogNormalize", scale.factor = 10000)
retinal <-  Seurat::FindVariableFeatures(retinal, selection.method = "vst", nfeatures = 2000)
retinal <-  Seurat::ScaleData(retinal)
retinal <- Seurat::RunPCA(retinal, features = Seurat::VariableFeatures(retinal),
                          verbose = F)

mat <- retinal[["RNA"]]@counts[Seurat::VariableFeatures(retinal),]
set.seed(10)
K <- 30
time_start <- Sys.time()
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "nb",
                             ctl = list(verbose = T), minibatch = "stochastic")
time_end <- Sys.time()
pred_mat <- glmpca:::predict.glmpca(glmpca_res)

save.image("../../../../out/writeup6/writeup6_10x_mouseretinal_glmpca_nb.RData")


