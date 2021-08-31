rm(list=ls())

library(Seurat); library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_36501genes.RData")
mat <- as.matrix(mat)

library_size_vec <- rowSums(mat)
covariates <- cbind(1, log(library_size_vec))

n <- nrow(mat)
uniq_diagnos <- unique(metadata$diagnosis)
uniq_sex <- unique(metadata$sex)
uniq_indiv <- unique(metadata$individual)
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
for(i in uniq_indiv[-1]){
  tmp <- rep(0, n)
  tmp[which(metadata$individual == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
colnames(covariates) <- c("Intercept", "Log-library", uniq_diagnos[-1], uniq_sex[-1], paste0("i", uniq_indiv[-1]))

# now fix the covariates
protected_columns <- 1:4
unprotected_columns <- c(1:ncol(covariates))[-protected_columns]

tol <- 1e-4
for(i in unprotected_columns){
  print(i)
  tmp <- stats::lm.fit(x = covariates[,protected_columns,drop=F],
                       y = covariates[,i,drop=F])
  if(max(abs(tmp$residuals)) <= tol) break()

  covariates[,i] <- tmp$residuals
  colnames(covariates)[i] <- paste0("adjusted_", colnames(covariates)[i])
  protected_columns <- c(protected_columns, i)
}
if(length(protected_columns) != ncol(covariates)){
  covariates <- covariates[,protected_columns,drop = F]
}

#########

# now kick out the genes that are too-few expressed
binary_mat <- mat
binary_mat[binary_mat > 0] <- 1
tmp <- matrixStats::colSums2(binary_mat)
idx <- which(tmp < n/100)
if(length(idx) > 0){
  mat <- mat[,-idx]
}

rm(list = c("binary_mat", "i", "idx", "library_size_vec",
            "protected_columns", "tmp", "tol", "uniq_diagnos",
            "uniq_indiv", "uniq_sex", "unprotected_columns", "zz"))

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
save.image("../../../../out/writeup7/writeup7_sns_esvd_nb_part1.RData")

print("Starting estimation")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)
time_end2 <- Sys.time()

save.image("../../../../out/writeup7/writeup7_sns_esvd_nb_part1.RData")


