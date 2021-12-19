rm(list=ls())
load("../../../../out/writeup8e/writeup8e_sns_pseudoreal_data.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

gene_frequency <- apply(mat, 2, function(x){length(which(x > 0))})/nrow(mat)
quantile(gene_frequency)
gene_threshold <- 0.2
gene_maxval <- apply(mat, 2, function(x){
  quantile(x[x>0], probs = 0.95)
})
gene_threshold2 <- 30
gene_idx <- intersect(which(gene_frequency < gene_threshold), which(gene_maxval < gene_threshold2))
length(gene_idx)
mat2 <- mat
zero_substitute <- 0.25
mat2[mat2 == 0] <- zero_substitute
mat_subset <- mat2[,gene_idx]

K <- 50
n <- nrow(mat_subset)
p1 <- ncol(mat_subset)

time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat_subset,
                                   k = K,
                                   family = "poisson",
                                   covariates = covariates,
                                   column_set_to_one = "Log_UMI",
                                   offset_vec = rep(0, nrow(mat)),
                                   verbose = 1)
time_end1 <- Sys.time()

print("Starting Poisson fit")
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
                            global_estimate = F,
                            l2pen = 0.01,
                            max_iter = 50,
                            reparameterize = T,
                            reestimate_nuisance = F,
                            verbose = 1)
time_end2 <- Sys.time()

#############################

# initialize the vectors for the other cells
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

zz <- tmp$residual_mat
svd_res <- irlba::irlba(zz, nv = 20)

init_x_mat <- cbind(init_x_mat, eSVD2:::.mult_mat_vec(svd_res$u[,1:5], sqrt(svd_res$d[1:5])))
tmp <- matrix(0, nrow = nrow(init_y_mat), ncol = ncol(init_y_mat) + 5)
tmp[,1:ncol(init_y_mat)] <- init_y_mat
tmp[-gene_idx,(ncol(init_y_mat)+1):(ncol(init_y_mat)+5)] <-  eSVD2:::.mult_mat_vec(svd_res$v[,1:5], sqrt(svd_res$d[1:5]))
init_y_mat <- tmp

##############

print("Starting large Poisson fit")
time_start3 <- Sys.time()
set.seed(10)
esvd_res_full <- eSVD2::opt_esvd(init_x_mat,
                                 init_y_mat,
                                 mat2,
                                 family = "poisson",
                                 nuisance_param_vec = NA,
                                 library_size_vec = 1,
                                 method = "newton",
                                 b_init = init_b_mat,
                                 covariates = init_res$covariates,
                                 global_estimate = F,
                                 l2pen = 0.01,
                                 max_iter = 50,
                                 reparameterize = T,
                                 reestimate_nuisance = F,
                                 verbose = 1)
time_end3 <- Sys.time()
save(date_of_run, session_info, gene_idx, de_genes,
     metadata, mat, esvd_res_full, time_start3, time_end3,
     covariates,
     file = "../../../../out/writeup8e/writeup8e_sns_pseudoreal_esvd_poisson.RData")
