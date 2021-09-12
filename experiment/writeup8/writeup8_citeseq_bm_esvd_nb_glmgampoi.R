rm(list=ls())

library(eSVD2)
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

binary_mat <- mat
binary_mat[binary_mat > 0] <- 1
tmp <- matrixStats::colSums2(binary_mat)
idx <- which(tmp < nrow(binary_mat)/100)
if(length(idx) > 0){
  mat <- mat[,-idx]
}
rm(list = c("binary_mat", "tmp", "idx"))

print("Initializing via glmGamPoi")
time_start1 <- Sys.time()
set.seed(10)
glmGamPoi_res <- glmGamPoi::glm_gp(t(mat),
                                   design = ~1,
                                   col_data = NULL,
                                   size_factors = "deconvolution",
                                   overdispersion = TRUE,
                                   overdispersion_shrinkage = FALSE,
                                   do_cox_reid_adjustment = TRUE,
                                   subsample = FALSE,
                                   on_disk = NULL,
                                   verbose = TRUE)
time_end1 <- Sys.time()

quantile(glmGamPoi_res$overdispersions)
quantile(glmGamPoi_res$size_factors)

nuisance_vec <- 1/glmGamPoi_res$overdispersions
nuisance_vec <- pmin(nuisance_vec, 1e5)
library_size_vec <- glmGamPoi_res$size_factors
save.image("../../../../out/writeup8/writeup8_citeseq_bm_esvd_nb_glmgampoi.RData")

K <- 30
n <- nrow(mat)
covariates <- cbind(1, log(library_size_vec))
colnames(covariates) <- c("Intercept", "Log-UMI")

print("Initializing NB")
time_start2 <- Sys.time()
set.seed(10)
init <- eSVD2::initialize_esvd(mat, k = K,
                               family = "neg_binom",
                               nuisance_param_vec = nuisance_vec,
                               library_size_vec = 1,
                               covariates = covariates,
                               check_rank = F,
                               config = eSVD2::initialization_options(),
                               verbose = 1)
time_end2 <- Sys.time()
save.image("../../../../out/writeup8/writeup8_citeseq_bm_esvd_nb_glmgampoi.RData")

print("Estimating NB")
time_start3 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat,
                            family = "neg_binom",
                            nuisance_param_vec = nuisance_vec,
                            library_size_vec = 1,
                            b_init = init$b_mat,
                            covariates = covariates,
                            max_iter = 100,
                            verbose = 2)
time_end3 <- Sys.time()
print("Finished")
save.image("../../../../out/writeup8/writeup8_citeseq_bm_esvd_nb_glmgampoi.RData")

