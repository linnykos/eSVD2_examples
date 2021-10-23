rm(list=ls())
library(Seurat); library(eSVD2)
load("sns_layer23_example.RData")

mat <- as.matrix(Matrix::t(sns_layer23[["RNA"]]@counts))
print(dim(mat))

categorical_var <- c("diagnosis", "individual", "region", "sex", "Seqbatch") #, "individual")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
n <- ncol(sns_layer23)
covariates <- as.matrix(sns_layer23@meta.data[,numerical_var])
covariates <- cbind(1, log(matrixStats::rowMeans2(mat)), covariates)
colnames(covariates)[1:2] <- c("Intercept", "Log-UMI")

for(variable in categorical_var){
  vec <- sns_layer23@meta.data[,variable]
  uniq_level <- unique(vec)
  for(i in uniq_level[-1]){
    tmp <- rep(0, n)
    tmp[which(vec == i)] <- 1

    var_name <- paste0(variable, "_", i)
    covariates <- cbind(covariates, tmp)
    colnames(covariates)[ncol(covariates)] <- var_name
  }
}
head(covariates)
dim(covariates)

##############
# initialization
K <- 10
n <- nrow(mat)
p <- ncol(mat)

# This takes about 5 minutes to complete -- I know a better
#   way to speed this up since the bottleneck is the repeated
#   QR factorization to regress out each gene individually.
#   I'll fix this by next time
time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat,
                                   k = K,
                                   family = "neg_binom2",
                                   covariates = covariates,
                                   column_set_to_one = "Log-UMI",
                                   offset_vec = rep(0, nrow(mat)),
                                   verbose = 1)
time_end1 <- Sys.time()

###################3

print("Estimating NB via eSVD, round 1")
time_start2 <- Sys.time()
# The following line is where the method crashes due to Heissian singularity
# If we did not include "individual" in Line 7, then the resulting covariate matrix
#   would not contain individual indicators, and then the method would run smoothly
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
                            init_res$y_mat,
                            mat,
                            family = "neg_binom2",
                            nuisance_param_vec = init_res$nuisance_param_vec,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = init_res$b_mat,
                            covariates = init_res$covariates,
                            offset_vec = init_res$offset_vec,
                            reestimate_nuisance = T,
                            global_estimate = T,
                            reparameterize = T,
                            max_iter = 50,
                            verbose = 2)
time_end2 <- Sys.time()

print("Estimating NB via eSVD, round 2")
time_start3 <- Sys.time()
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(esvd_res$x_mat,
                             esvd_res$y_mat, mat,
                             family = "neg_binom2",
                             nuisance_param_vec = esvd_res$nuisance_param_vec,
                             library_size_vec = 1,
                             method = "newton",
                             b_init = esvd_res$b_mat,
                             covariates = esvd_res$covariates,
                             offset_vec = esvd_res$offset_vec,
                             reestimate_nuisance = T,
                             global_estimate = F,
                             reparameterize = T,
                             max_iter = 50,
                             tol = 1e-8,
                             verbose = 1)
time_end3 <- Sys.time()

