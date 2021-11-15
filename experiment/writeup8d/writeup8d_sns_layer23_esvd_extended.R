# rm(list=ls())
# load("../../../../out/writeup8d/writeup8d_sns_layer23_processed.RData")
#
# library(Seurat)
# set.seed(10)
# date_of_run <- Sys.time()
# session_info <- devtools::session_info()
#
# gene_frequency <- apply(mat, 2, function(x){length(which(x > 0))})/nrow(mat)
# quantile(gene_frequency)
# gene_threshold <- 0.8
# gene_idx <- which(gene_frequency > gene_threshold)
# length(intersect(colnames(mat)[gene_idx], de_genes))
# length(de_genes)
# mat_subset <- mat[,gene_idx]
#
# K <- 15
# n <- nrow(mat_subset)
# p1 <- ncol(mat_subset)
#
# time_start1 <- Sys.time()
# init_res <- eSVD2::initialize_esvd(mat_subset,
#                                    k = K,
#                                    family = "poisson",
#                                    covariates = covariates,
#                                    column_set_to_one = "Log-UMI",
#                                    offset_vec = rep(0, nrow(mat)),
#                                    verbose = 1)
# time_end1 <- Sys.time()
#
# print("Starting Poisson fit")
# time_start2 <- Sys.time()
# set.seed(10)
# esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
#                             init_res$y_mat,
#                             mat_subset,
#                             family = "poisson",
#                             nuisance_param_vec = NA,
#                             library_size_vec = 1,
#                             method = "newton",
#                             b_init = init_res$b_mat,
#                             covariates = init_res$covariates,
#                             offset_vec = init_res$offset_vec,
#                             reestimate_nuisance = F,
#                             global_estimate = F,
#                             reparameterize = T,
#                             max_iter = 50,
#                             l2pen = 0.1,
#                             verbose = 1)
# time_end2 <- Sys.time()
# save.image("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended.RData")
#
# #############
#
# print("Computing NB initialization")
# # assess the fit
# nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
# nat_mat2 <- tcrossprod(esvd_res$covariates, esvd_res$b_mat)
# nat_mat <- nat_mat1 + nat_mat2
# mean_mat <- exp(nat_mat)
# nuisance_param_vec <- sapply(1:p1, function(j){
#   MASS::theta.mm(y = mat_subset[,j], mu = mean_mat[,j], dfr = nrow(mat_subset)-1)
# })
# nuisance_param_vec <- pmax(pmin(nuisance_param_vec, 1000), 1)
#
# rm(list = c("nat_mat1", "nat_mat2", "mean_mat", "mat_subset"))
#
# ##########
#
# # initialize the vectors for the other cells
# tmp <- eSVD2:::.regress_out_matrix(log(mat[,-gene_idx] + 1e-3),
#                                    covariates = cbind(esvd_res$x_mat, esvd_res$covariates),
#                                    verbose = 1)
# colnames(tmp$b_mat) <- colnames(cbind(esvd_res$x_mat, esvd_res$covariates))
# rownames(tmp$b_mat) <- colnames(mat[,-gene_idx])
# rownames(tmp$residual_mat) <- rownames(mat[,-gene_idx])
# colnames(tmp$residual_mat) <- colnames(mat[,-gene_idx])
#
# init_x_mat <- esvd_res$x_mat
# init_y_mat <- matrix(0, nrow = ncol(mat), ncol = ncol(esvd_res$y_mat))
# init_y_mat[gene_idx,] <- esvd_res$y_mat
# init_y_mat[-gene_idx,] <- tmp$b_mat[,1:K]
# rownames(init_y_mat) <- colnames(mat)
# colnames(init_y_mat) <- colnames(esvd_res$y_mat)
# covariates <- esvd_res$covariates
# init_b_mat <- matrix(0, nrow = ncol(mat), ncol = ncol(esvd_res$b_mat))
# init_b_mat[gene_idx,] <- esvd_res$b_mat
# init_b_mat[-gene_idx,] <- tmp$b_mat[,-c(1:K)]
# rownames(init_b_mat) <- colnames(mat)
# colnames(init_b_mat) <- colnames(esvd_res$b_mat)
#
# zz <- tmp$residual_mat
# svd_res <- irlba::irlba(zz, nv = 20)
#
# init_x_mat <- cbind(init_x_mat, eSVD2:::.mult_mat_vec(svd_res$u[,1:5], sqrt(svd_res$d[1:5])))
# tmp <- matrix(0, nrow = nrow(init_y_mat), ncol = ncol(init_y_mat) + 5)
# tmp[,1:ncol(init_y_mat)] <- init_y_mat
# tmp[-gene_idx,(ncol(init_y_mat)+1):(ncol(init_y_mat)+5)] <-  eSVD2:::.mult_mat_vec(svd_res$v[,1:5], sqrt(svd_res$d[1:5]))
# init_y_mat <- tmp
#
# rm(list = c("tmp", "zz", "svd_res"))
# save.image("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended.RData")

##########################
load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended.RData")

gene_group_factor <- rep("overdispersed_2", ncol(mat))
nuisance_param_vec_full <- rep(0.5, ncol(mat))
gene_group_factor[gene_idx] <- paste0("normal_", 1:length(gene_idx))
nuisance_param_vec_full[gene_idx] <- nuisance_param_vec
gene_frequency <- apply(mat, 2, function(x){length(which(x > 0))})/nrow(mat)
gene_group_factor[which(gene_frequency < 0.4)] <- "overdispersed_1"
nuisance_param_vec_full[which(gene_frequency < 0.4)] <- 0.1
gene_group_factor <- factor(gene_group_factor)
nuisance_value_lower <- c(rep(1, length(gene_idx)), 0.1, 0.1)
nuisance_value_upper <- c(rep(1000, length(gene_idx)), 1, 1)

print("Starting NB fit")
time_start3 <- Sys.time()
set.seed(10)
esvd_res_nb1 <- eSVD2::opt_esvd(init_x_mat,
                                init_y_mat,
                                mat,
                                family = "neg_binom2",
                                nuisance_param_vec = nuisance_param_vec_full,
                                offset_vec = init_res$offset_vec,
                                library_size_vec = 1,
                                method = "newton",
                                b_init = init_b_mat,
                                covariates = init_res$covariates,
                                bool_run_cpp = F,
                                gene_group_factor = gene_group_factor,
                                gene_ignore_excessive_zero = T,
                                nuisance_value_lower = nuisance_value_lower,
                                nuisance_value_upper = nuisance_value_upper,
                                reestimate_nuisance = T,
                                reestimate_nuisance_per_iteration = 10,
                                reparameterize = T,
                                max_iter = 50,
                                l2pen = 0.1,
                                verbose = 1)
time_end3 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended.RData")

gene_group_factor <- rep(NA, ncol(mat))
gene_group_factor[gene_idx] <- paste0("normal_", 1:length(gene_idx))
gene_group_factor[-gene_idx] <- paste0("overdispersed_", 1:(ncol(mat)-length(gene_idx)))
nuisance_value_lower <- rep(0.1, ncol(mat))
nuisance_value_lower[1:length(gene_idx)] <- 1
nuisance_value_upper <- rep(1000, ncol(mat))
time_start4 <- Sys.time()
set.seed(10)
esvd_res_nb2 <- eSVD2::opt_esvd(esvd_res_nb1$x_mat,
                                esvd_res_nb1$y_mat,
                                mat,
                                family = "neg_binom2",
                                nuisance_param_vec = esvd_res_nb1$nuisance_param_vec,
                                offset_vec = esvd_res_nb1$offset_vec,
                                library_size_vec = 1,
                                method = "newton",
                                b_init = esvd_res_nb1$b_mat,
                                covariates = esvd_res_nb1$covariates,
                                bool_run_cpp = F,
                                gene_group_factor = gene_group_factor,
                                gene_ignore_excessive_zero = T,
                                nuisance_value_lower = nuisance_value_lower,
                                nuisance_value_upper = nuisance_value_upper,
                                reestimate_nuisance = T,
                                reestimate_nuisance_per_iteration = 10,
                                reparameterize = T,
                                max_iter = 50,
                                l2pen = 0.1,
                                verbose = 1)
time_end4 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended.RData")


