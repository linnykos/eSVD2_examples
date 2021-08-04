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

mat <- retinal[["RNA"]]@counts[Seurat::VariableFeatures(retinal),]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

# run eSVD
print("Estimating Poisson")
set.seed(10)
n <- nrow(mat)
library_size_vec <- rowSums(mat)
covariates <- cbind(matrix(1, nrow = n, ncol = 1), log(library_size_vec))
K <- 30
init <- eSVD2::initialize_esvd(mat, k = K, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = 1,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)

save.image("../../../../out/writeup6b/writeup6b_10x_mouseretinal_esvd.RData")

