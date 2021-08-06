rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
print("Loading in data")
dat <- anndata::read_h5ad("../../../../data/10x_mousepancreas/endocrinogenesis_day15.5.h5ad")
tmp <- Matrix::t(dat$X)
tmp <- Matrix::Matrix(as.matrix(tmp), sparse = T)
clusters <- dat$obs$clusters
rm(list = "dat")
gc()

print("Starting Seurat")
pancreas <- Seurat::CreateSeuratObject(counts = tmp)
pancreas[["celltype"]] <- clusters
pancreas <- Seurat::NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas <-  Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)

# extract data
mat <- pancreas[["RNA"]]@counts[Seurat::VariableFeatures(pancreas),]
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

save.image("../../../../out/writeup6b/writeup6b_10x_mousepancreas_esvd.RData")

###############3

print("Estimating NB")
init_nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
nuisance_vec <- eSVD2::initialize_nuisance_param(mat, init_nat_mat, family = "neg_binom",
                                          library_size_vec = rep(1, n))

init <- eSVD2::initialize_esvd(mat, k = K, family = "neg_binom",
                               nuisance_param_vec = nuisance_vec,
                               library_size_vec = 1,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "neg_binom",
                            nuisance_param_vec = nuisance_vec, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)

save.image("../../../../out/writeup6b/writeup6b_10x_mousepancreas_esvd.RData")

