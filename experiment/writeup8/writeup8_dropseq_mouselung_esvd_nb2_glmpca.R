rm(list=ls())
load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")

library(eSVD2); library(glmGamPoi); library(scran); library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(Matrix::t(mat))
nuisance_vec <- rep(glmpca_res$glmpca_family$nb_theta[1], ncol(mat))
library_size_vec <- matrixStats::rowSums2(mat)
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb2_glmpca.RData")

K <- 30
n <- nrow(mat)
covariates <- cbind(1, log(library_size_vec))
colnames(covariates) <- c("Intercept", "Log-UMI")

print("Initializing NB")
time_start2 <- Sys.time()
set.seed(10)
init <- eSVD2::initialize_esvd(mat, k = K,
                               family = "neg_binom2",
                               nuisance_param_vec = nuisance_vec,
                               library_size_vec = 1,
                               covariates = covariates,
                               check_rank = F,
                               config = eSVD2::initialization_options(),
                               verbose = 1)
time_end2 <- Sys.time()
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb2_glmpca.RData")

print("Estimating NB")
time_start3 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat,
                            family = "neg_binom2",
                            nuisance_param_vec = nuisance_vec,
                            library_size_vec = 1,
                            b_init = init$b_mat,
                            covariates = covariates,
                            max_iter = 100,
                            tol = 1e-8,
                            verbose = 2)
time_end3 <- Sys.time()
print("Finished")
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb2_glmpca.RData")

