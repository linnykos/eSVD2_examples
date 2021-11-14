rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_layer23_processed.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

gene_frequency <- apply(mat, 2, function(x){length(which(x > 0))})/nrow(mat)
quantile(gene_frequency)
gene_threshold <- 0.8
gene_idx <- which(gene_frequency > gene_threshold)
length(intersect(colnames(mat)[gene_idx], de_genes))
length(de_genes)
mat_subset <- mat[,gene_idx]

K <- 15
n <- nrow(mat_subset)
p1 <- ncol(mat_subset)

time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat_subset,
                                   k = K,
                                   family = "poisson",
                                   covariates = covariates,
                                   column_set_to_one = "Log-UMI",
                                   offset_vec = rep(0, nrow(mat)),
                                   verbose = 1)
time_end1 <- Sys.time()

time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
                            init_res$y_mat,
                            mat_subset,
                            family = "poisson",
                            nuisance_param_vec = NA,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = init_res$b_mat,
                            covariates = init_res$covariates,
                            offset_vec = init_res$offset_vec,
                            reestimate_nuisance = F,
                            global_estimate = F,
                            reparameterize = T,
                            max_iter = 50,
                            l2pen = 0.1,
                            verbose = 1)
time_end2 <- Sys.time()

#############

# assess the fit
nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates, esvd_res$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
angle_vec <- sapply(1:p1, function(j){
  tmp <- cbind(mat_subset[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)
deviation_vec <- sapply(1:p1, function(j){
  residual_vec <- mat_subset[,j] - mean_mat[,j]
  sd_vec <- 3*sqrt(mean_mat[,j])
  length(which(abs(residual_vec) <= sd_vec))/length(residual_vec)
})
quantile(deviation_vec)

esvd_res$b_mat[which(rownames(esvd_res$b_mat) %in% de_genes), "diagnosis_ASD"]
quantile(esvd_res$b_mat[,"diagnosis_ASD"])

# nuisance_param_vec <- glmGamPoi::overdispersion_mle(y = t(mat_subset),
#                                                     mean = t(mean_mat),
#                                                     global_estimate = F)
nuisance_param_vec <- sapply(1:p1, function(j){
  # MASS::theta.ml(y = mat_subset[,j], mu = mean_mat[,j])
  MASS::theta.mm(y = mat_subset[,j], mu = mean_mat[,j], dfr = nrow(mat_subset)-1)
})
nuisance_param_vec <- pmax(pmin(nuisance_param_vec, 500), 1)
quantile(nuisance_param_vec)

##########

# initialize the vectors for the other cells
colnames(esvd_res$b_mat)[which(colnames(esvd_res$b_mat) == "Log-UMI")] <- "Log_UMI"
colnames(esvd_res$covariates)[which(colnames(esvd_res$covariates) == "Log-UMI")] <- "Log_UMI"

tmp <- eSVD2:::.regress_out_matrix(log(mat[,-gene_idx] + 1e-3),
                                   covariates = cbind(esvd_res$x_mat, esvd_res$covariates),
                                   verbose = 1)
colnames(tmp$b_mat) <- colnames(cbind(esvd_res$x_mat, esvd_res$covariates))
rownames(tmp$b_mat) <- colnames(mat[,-gene_idx])
rownames(tmp$residual_mat) <- rownames(mat[,-gene_idx])
colnames(tmp$residual_mat) <- colnames(mat[,-gene_idx])

init_x_mat <- esvd_res$x_mat
init_y_mat <- matrix(0, nrow = ncol(mat), ncol = ncol(esvd_res$y_mat))
init_y_mat[gene_idx,] <- esvd_res$y_mat
init_y_mat[-gene_idx,] <- tmp$b_mat[,1:K]
rownames(init_y_mat) <- colnames(mat)
colnames(init_y_mat) <- colnames(esvd_res$y_mat)
covariates <- esvd_res$covariates
init_b_mat <- matrix(0, nrow = ncol(mat), ncol = ncol(esvd_res$b_mat))
init_b_mat[gene_idx,] <- esvd_res$b_mat
init_b_mat[-gene_idx,] <- tmp$b_mat[,-c(1:K)]
rownames(init_b_mat) <- colnames(mat)
colnames(init_b_mat) <- colnames(esvd_res$b_mat)

# try adding extra dimensions
zz <- tmp$residual_mat
svd_res <- irlba::irlba(zz, nv = 7)

nat_mat1 <- tcrossprod(init_x_mat, init_y_mat)
nat_mat2 <- tcrossprod(covariates, init_b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
angle_vec <- sapply(1:ncol(mean_mat), function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)
quantile(angle_vec[gene_idx])

j <- which.min(abs(angle_vec - 89)) # which.max(angle_vec)
quantile(mat[,j])
stats::cor(mean_mat[,j], mat[,j])
tmp <- cbind(mat[,j], mean_mat[,j])
eSVD2:::.compute_principal_angle(tmp)
idx <- which(mat[,j] != 0)
stats::cor(mean_mat[idx,j], mat[idx,j])
tmp <- cbind(mat[idx,j], mean_mat[idx,j])
eSVD2:::.compute_principal_angle(tmp)
quantile(mean_mat[,j])
quantile(mean_mat[idx,j])

