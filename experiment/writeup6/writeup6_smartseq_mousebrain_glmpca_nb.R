rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
load("../../../../data/smartseq_mousebrain/smartseq_mousebrain_formatted.RData")


brain <- Seurat::CreateSeuratObject(counts = Matrix::t(dat), meta.data = metadata, min.cells = 10)
low_q_cells <- rownames(brain@meta.data[brain@meta.data$class %in% c('Low Quality', 'No Class'), ])
ok_cells <- rownames(brain@meta.data)[!(rownames(x = brain@meta.data) %in% low_q_cells)]
brain <- brain[, ok_cells]
rm(list = "dat")
gc()

brain <- Seurat::NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
brain <-  Seurat::FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <-  Seurat::ScaleData(brain)
brain <- Seurat::RunPCA(brain, features = Seurat::VariableFeatures(brain),
                        verbose = F)

mat <- brain[["RNA"]]@counts[Seurat::VariableFeatures(brain),]
set.seed(10)
K <- 30
time_start <- Sys.time()
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "nb",
                             ctl = list(verbose = T), minibatch = "stochastic")
time_end <- Sys.time()
pred_mat <- glmpca:::predict.glmpca(glmpca_res)

save.image("../../../../out/writeup6/writeup6_smartseq_mousebrain_glmpca_nb.RData")


