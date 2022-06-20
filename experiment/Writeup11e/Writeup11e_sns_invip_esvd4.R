rm(list=ls())
library(Seurat)

load("../../../../out/Writeup11e/Writeup11e_sns_invip_esvd.RData")
# load("../../../../out/Writeup10/Writeup10_sns_invip_esvd2.RData")
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")
source("initialization.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# doing the initialization as in Writeup10, but now with an intercept
dat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
init_res <- initialize_esvd2(dat = dat,
                             k = 50,
                             covariates = covariates,
                             verbose = 1)

# replace eSVD_obj with the appropriate things
eSVD_obj$param$init_k <- 50
eSVD_obj$fit_First <- NULL
eSVD_obj$fit_Second <- NULL
covariates <- cbind(init_res$covariates, init_res$offset_vec)
colnames(covariates)[ncol(covariates)] <- "Log_UMI"
eSVD_obj$covariates <- covariates
eSVD_obj$fit_Init$x_mat <- init_res$x_mat
eSVD_obj$fit_Init$y_mat <- init_res$y_mat
eSVD_obj$fit_Init$z_mat <- cbind(init_res$b_mat, 1)
colnames(eSVD_obj$fit_Init$z_mat)[ncol(eSVD_obj$fit_Init$z_mat)] <- "Log_UMI"

print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.01,
                             max_iter = 100,
                             offset_variables = "Log_UMI",
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_First",
                             fit_previous = "fit_Init")
time_end2 <- Sys.time()

print("Second fit")
time_start3 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.01,
                             max_iter = 100,
                             offset_variables = NULL,
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_Second",
                             fit_previous = "fit_First")
time_end3 <- Sys.time()

print("Nuisance estimation")
time_start4 <- Sys.time()
eSVD_obj <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                      bool_covariates_as_library = T,
                                      bool_library_includes_interept = T,
                                      verbose = 1)
time_end4 <- Sys.time()


eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      bool_covariates_as_library = T)
metadata <- sns@meta.data
metadata[,"individual"] <- as.factor(metadata[,"individual"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "individual",
                                           metadata = metadata,
                                           verbose = 1)
time_end5 <- Sys.time()

save(date_of_run, session_info, sns, covariate_df,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../out/Writeup11e/Writeup11e_sns_invip_esvd4.RData")




