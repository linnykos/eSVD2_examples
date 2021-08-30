rm(list=ls())

library(Seurat); library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_36501genes.RData")
mat <- as.matrix(mat)

print("Forming covariates")
set.seed(10)
n <- nrow(mat)
library_size_vec <- rowSums(mat)
covariates <- cbind(matrix(1, nrow = n, ncol = 1), log(library_size_vec))
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

print("Shuffle covariates")
uniq_indiv <- unique(metadata$individual)
autism_indiv <- unique(metadata$individual[metadata$diagnosis == "ASD"])
set.seed(10)
autism_indiv2 <- sample(uniq_indiv, size = length(autism_indiv))
covariates[,"ASD"] <- rep(0, n)
covariates[which(metadata$individual %in% autism_indiv2),"ASD"] <- 1

K <- 5
print("Starting initialization")
time_start1 <- Sys.time()
set.seed(10)
init <- eSVD2::initialize_esvd(mat, k = K, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = 1,
                               covariates = covariates,
                               check_rank = F,
                               config = eSVD2::initialization_options(), verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_shuffled.RData")

print("Starting estimation")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)
time_end2 <- Sys.time()

save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_shuffled.RData")


