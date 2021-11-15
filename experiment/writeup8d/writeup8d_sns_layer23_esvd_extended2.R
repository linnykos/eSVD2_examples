# try something similar, but with gene_ignore_excessive_zero turned on
rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended_intermediary.RData")
date_of_run <- Sys.time()
session_info <- devtools::session_info()

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
                                reestimate_nuisance_per_iteration = 5,
                                reparameterize = F,
                                max_iter = 50,
                                l2pen = 0.1,
                                verbose = 1)
time_end3 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended2.RData")

gene_group_factor <- rep(NA, ncol(mat))
gene_group_factor[gene_idx] <- paste0("normal_", 1:length(gene_idx))
gene_group_factor[-gene_idx] <- paste0("overdispersed_", 1:(ncol(mat)-length(gene_idx)))
gene_group_factor <- factor(gene_group_factor)
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
                                reestimate_nuisance_per_iteration = 5,
                                reparameterize = F,
                                max_iter = 50,
                                l2pen = 0.1,
                                verbose = 1)
time_end4 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_extended2.RData")


