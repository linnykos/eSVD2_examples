rm(list=ls())

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

##############
# initialization
K <- 30
n <- nrow(mat)
p <- ncol(mat)

covariates <- cbind(1, log(matrixStats::rowMeans2(mat)))
colnames(covariates) <- c("Intercept", "Log-UMI")
b_init <- cbind(log(matrixStats::colMeans2(mat)), 1)
colnames(b_init) <- c("Intercept", "Log-UMI")
nuisance_vec <- rep(100, p)

set.seed(10)
nat_mat <- log1p(mat)
residual_mat <- nat_mat - tcrossprod(covariates, b_init)
svd_res <- irlba::irlba(residual_mat, nv = K)
x_init <- eSVD2:::.mult_mat_vec(svd_res$u, sqrt(svd_res$d))
y_init <- eSVD2:::.mult_mat_vec(svd_res$v, sqrt(svd_res$d))

# tmp <- tcrossprod(x_init, y_init) + tcrossprod(covariates, b_init)
################

print("Estimating NB via eSVD, w/ reparameterization and gene-specific nuisance")
time_start1 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(x_init, y_init, mat,
                            family = "neg_binom2",
                            nuisance_param_vec = nuisance_vec,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = b_init,
                            covariates = covariates,
                            reestimate_nuisance = T,
                            global_estimate = F,
                            reparameterize = T,
                            max_iter = 100,
                            tol = 1e-8,
                            verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup8b/writeup8b_dropseq_mouselung_esvd2.RData")
