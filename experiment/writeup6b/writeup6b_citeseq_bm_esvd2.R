rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
bm <- SeuratData::LoadData(ds = "bmcite")

Seurat::DefaultAssay(bm) <- "RNA"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm <-  Seurat::FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)

mat <- bm[["RNA"]]@counts[Seurat::VariableFeatures(bm),]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

# run eSVD
print("Estimating Poisson")
set.seed(10)
n <- nrow(mat)
covariates <- matrix(1, nrow = n, ncol = 1)
K <- 30
init <- eSVD2::initialize_esvd(mat, k = K, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = NA,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = NA,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)

save.image("../../../../out/writeup6b/writeup6b_citeseq_bm_esvd2.RData")
