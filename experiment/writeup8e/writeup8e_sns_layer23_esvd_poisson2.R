rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_layer23_processed.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

K <- 30
n <- nrow(mat)
p <- ncol(mat)

offset_vec <- covariates[,"Log_UMI"]
covariates <- covariates[,-which(colnames(covariates) == "Log_UMI")]

time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat,
                                   k = K,
                                   family = "poisson",
                                   covariates = covariates,
                                   column_set_to_one = NULL,
                                   offset_vec = offset_vec,
                                   verbose = 1)
time_end1 <- Sys.time()

print("Starting Poisson fit")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
                            init_res$y_mat,
                            mat,
                            family = "poisson",
                            nuisance_param_vec = NA,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = init_res$b_mat,
                            covariates = init_res$covariates,
                            offset_vec = init_res$offset_vec,
                            reestimate_nuisance = F,
                            global_estimate = F,
                            reparameterize = F,
                            max_iter = 50,
                            l2pen = 0.1,
                            verbose = 1)
time_end2 <- Sys.time()

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates, esvd_res$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)

gene_frequency <- apply(mat, 2, function(x){length(which(x > 0))})/nrow(mat)
quantile(gene_frequency)
gene_threshold <- 0.8
gene_idx <- which(gene_frequency > gene_threshold)
quantile(angle_vec[gene_idx])

nuisance_param_vec <- sapply(1:p, function(j){
  MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
})

load("../../../../data/sns_autism/sns_formatted.RData")

metadata <- sns@meta.data
cell_names <- rownames(esvd_res$x_mat)
cell_idx <- sapply(cell_names, function(name){
  which(rownames(metadata) == name)
})
metadata <- metadata[cell_idx,]

save(covariates, esvd_res, mat, date_of_run, session_info, metadata,
     file = "../../../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson2.RData")

