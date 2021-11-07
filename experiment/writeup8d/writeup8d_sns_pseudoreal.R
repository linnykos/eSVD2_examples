# rm(list=ls())
# library(Seurat)
#
# load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd.RData")
# ls_vec <- ls()
# ls_vec <- ls_vec[!ls_vec %in% c("esvd_res2")]
# rm(list = ls_vec)
#
# n <- nrow(esvd_res2$x_mat)
# p <- nrow(esvd_res2$y_mat)
# autism_idx <- which(colnames(esvd_res2$covariates) == "diagnosis_ASD")
# esvd_res2$b_mat[,] <- 0
# nat_mat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(esvd_res2$covariates, esvd_res2$b_mat)
#
# set.seed(10)
# autism_gene_idx <- sample(1:p, size = round(p/100))
# multiplier <- 10
# for(j in autism_gene_idx){
#   vec <- nat_mat[,j]
#   vec_autism <- vec[which(esvd_res2$covariates[,"diagnosis_ASD"] > 0.5)]
#   vec_control <- vec[which(esvd_res2$covariates[,"diagnosis_ASD"] < 0.5)]
#
#   target <- multiplier * max(mean(exp(vec_control)), mean(exp(vec_autism)))
#   current <- mean(exp(vec_autism))
#   value <-  log(target/current)
#   esvd_res2$b_mat[j,autism_idx] <- value
# }
#
# nat_mat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(esvd_res2$covariates, esvd_res2$b_mat)
# quantile(nat_mat[which(esvd_res2$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx])
# quantile(nat_mat[which(esvd_res2$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx])
#
# #################
#
# # now simulate data
# mat <- nat_mat
# for(j in 1:ncol(mat)){
#   set.seed(j)
#   mat[,j] <- stats::rnbinom(n, size = esvd_res2$nuisance_param_vec[j], mu = exp(nat_mat[,j]))
# }
# tmp <- mat[which(esvd_res2$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx]
# length(which(tmp == 0))/prod(dim(tmp))
# quantile(tmp[tmp > 0])
#
# tmp <- mat[which(esvd_res2$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx]
# length(which(tmp == 0))/prod(dim(tmp))
# quantile(tmp[tmp > 0])
#
# ######################
#
# # now fit
#
# # initialization
# K <- 10
# n <- nrow(mat)
# p <- ncol(mat)
# covariates <- esvd_res2$covariates
# esvd_res_truth <- esvd_res2
#
# time_start1 <- Sys.time()
# init_res <- eSVD2::initialize_esvd(mat,
#                                    k = K,
#                                    family = "neg_binom2",
#                                    covariates = covariates,
#                                    column_set_to_one = "Log-UMI",
#                                    offset_vec = rep(0, nrow(mat)),
#                                    verbose = 1)
# time_end1 <- Sys.time()
# save.image("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")
#
# quantile(init_res$b_mat[autism_gene_idx,autism_idx])
# quantile(init_res$b_mat[-autism_gene_idx,autism_idx])

###################3
load("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

print("Estimating NB via eSVD, round 1")
time_start2 <- Sys.time()
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
                            bool_run_cpp = F,
                            max_iter = 50,
                            l2pen = 0.1,
                            verbose = 1)
time_end2 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

quantile(esvd_res$b_mat[autism_gene_idx,autism_idx])
quantile(esvd_res$b_mat[-autism_gene_idx,autism_idx])

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
                             bool_run_cpp = F,
                             max_iter = 50,
                             tol = 1e-8,
                             l2pen = 0.1,
                             verbose = 1)
time_end3 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

quantile(esvd_res2$b_mat[autism_gene_idx,autism_idx])
quantile(esvd_res2$b_mat[-autism_gene_idx,autism_idx])




