rm(list=ls())
load("../../../../data/dropseq_humancortical/dropseq_humancortical_formatted.RData")

library(eSVD2)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

print("Loading in data")
dat <- anndata::read_h5ad("../../../../data/dropseq_mouselung/lung_regeneration_after_bleo")

print("Starting Seurat ")
tmp <- Matrix::t(dat$X); tmp <- as.matrix(tmp); tmp <- Matrix::Matrix(tmp, sparse = T)
lung <- Seurat::CreateSeuratObject(counts = tmp)
rm(list = "tmp")
lung[["celltype"]] <- dat$obs$clusters
lung <- Seurat::NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung <-  Seurat::FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)
rm(list = "dat")

mat <- lung[["RNA"]]@counts[Seurat::VariableFeatures(lung),]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

binary_mat <- mat
binary_mat[binary_mat > 0] <- 1
tmp <- matrixStats::colSums2(binary_mat)
idx <- which(tmp < nrow(binary_mat)/100)
if(length(idx) > 0){
  mat <- mat[,-idx]
}
rm(list = c("binary_mat", "tmp", "idx"))

# run eSVD
print("Initializing Poisson")
set.seed(10)
n <- nrow(mat)
library_size_vec <- rowSums(mat)
covariates <- cbind(matrix(1, nrow = n, ncol = 1), log(library_size_vec))
colnames(covariates) <- c("Intercept", "Log-UMI")

K <- 30
time_start1 <- Sys.time()
init <- eSVD2::initialize_esvd(mat, k = K, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = 1,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
time_end1 <- Sys.time()

print("Estimating Poisson")
time_start2 <- Sys.time()
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)
time_end2 <- Sys.time()

save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb.RData")

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
print("Estimating nusiance")
time_start3 <- Sys.time()
nuisance_vec <- eSVD2::initialize_nuisance_param(mat, nat_mat, family = "neg_binom",
                                                 library_size_vec = 1)
time_end3 <- Sys.time()
rm(list = "nat_mat")
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb.RData")

print("Initializing NB")
time_start4 <- Sys.time()
set.seed(10)
init2 <- eSVD2::initialize_esvd(mat, k = K, family = "neg_binom",
                                nuisance_param_vec = nuisance_vec,
                                library_size_vec = 1,
                                covariates = covariates,
                                check_rank = F,
                                config = eSVD2::initialization_options(), verbose = 1)
time_end4 <- Sys.time()
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb.RData")

print("Estimating NB")
time_start5 <- Sys.time()
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(init2$x_mat, init2$y_mat, mat, family = "neg_binom",
                             nuisance_param_vec = nuisance_vec,
                             library_size_vec = 1,
                             b_init = init2$b_mat,
                             covariates = covariates,
                             max_iter = 100,
                             verbose = 1)
time_end5 <- Sys.time()
print("Finished")
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb.RData")


