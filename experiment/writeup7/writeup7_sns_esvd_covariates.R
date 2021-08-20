rm(list=ls())
load("../../../../data/sns_autism/sns_formatted.RData")

mat <- sns[["RNA"]]@counts[Seurat::VariableFeatures(sns),]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

# run eSVD
print("Estimating Poisson")
set.seed(10)
n <- nrow(mat)
library_size_vec <- rowSums(mat)
covariates <- cbind(matrix(1, nrow = n, ncol = 1), log(library_size_vec))

# uniq_indiv <- unique(sns@meta.data$individual)
uniq_diagnos <- unique(sns@meta.data$diagnosis)
uniq_sex <- unique(sns@meta.data$sex)
# for(i in uniq_indiv){
#   tmp <- rep(0, n)
#   tmp[which(sns@meta.data$individual == i)] <- 1
#   covariates <- cbind(covariates, tmp)
# }
for(i in uniq_diagnos[-1]){
  tmp <- rep(0, n)
  tmp[which(sns@meta.data$diagnosis == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
for(i in uniq_sex[-1]){
  tmp <- rep(0, n)
  tmp[which(sns@meta.data$sex == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
#colnames(covariates) <- c("Intercept", "Log-library", uniq_indiv, uniq_diagnos, uniq_sex)
colnames(covariates) <- c("Intercept", "Log-library", uniq_diagnos[-1], uniq_sex[-1])


K <- 30
print("Starting initialization")
init <- eSVD2::initialize_esvd(mat, k = K, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = 1,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates.RData")

print("Starting estimation")
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)

save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates.RData")
