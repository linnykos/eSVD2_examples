rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

library(Seurat)

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_large.RData")

print("Forming covariates")
mat <- sns[["RNA"]]@counts[Seurat::VariableFeatures(sns), which(sns@meta.data$celltype == "L2/3")]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

set.seed(10)
n <- nrow(mat)
library_size_vec <- rowSums(mat)
covariates <- cbind(matrix(1, nrow = n, ncol = 1), log(library_size_vec))
metadata <- sns@meta.data
metadata <- metadata[which(metadata$celltype == "L2/3"),]
uniq_diagnos <- unique(metadata$diagnosis)
uniq_sex <- unique(metadata$sex)
for(i in uniq_diagnos[-1]){
  tmp <- rep(0, n)
  tmp[which(metadata$diagnosis == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
for(i in uniq_sex[-1]){
  tmp <- rep(0, n)
  tmp[which(metadata$sex == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
colnames(covariates) <- c("Intercept", "Log-library", uniq_diagnos[-1], uniq_sex[-1])

K <- 30
print("Starting initialization")
time_start1 <- Sys.time()
set.seed(10)
init <- eSVD2::initialize_esvd(mat, k = K, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = 1,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates_large2.RData")

print("Starting estimation")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)
time_end2 <- Sys.time()

save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates_large2.RData")


