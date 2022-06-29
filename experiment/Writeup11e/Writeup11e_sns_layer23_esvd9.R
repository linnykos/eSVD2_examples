# try no intercept, but rescaling many covariates besides nFeature_RNA
# but not Log_UMI

rm(list=ls())
library(Seurat)

load("../../../../out/Writeup11c/Writeup11c_sns_layer23_esvd.RData")
load("../../../../out/Writeup10/Writeup10_sns_layer23_processed2.RData")
source("initialization.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# doing the initialization as in Writeup10, but now with an intercept
dat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
init_res <- initialize_esvd2(dat = dat,
                             k = 50,
                             covariates = covariates,
                             rescale_vars = c("age", "nFeature_RNA", "percent.mt", "post.mortem.hours"),
                             bool_intercept = F,
                             verbose = 1)

# replace eSVD_obj with the appropriate things
eSVD_obj$param$init_bool_intercept <- F
eSVD_obj$param$init_k <- 50

eSVD_obj$param$init_lambda <- NULL
eSVD_obj$param$init_mixed_effect_variables <- NULL
eSVD_obj$param$fit_First_l2pen <- NULL
eSVD_obj$param$fit_First_max_iter <- NULL
eSVD_obj$param$fit_First_offset_variables <- NULL
eSVD_obj$param$fit_First_tol <- NULL
eSVD_obj$param$fit_Second_family <- NULL
eSVD_obj$param$fit_Second_l2pen <- NULL
eSVD_obj$param$fit_Second_max_iter <- NULL
eSVD_obj$param$fit_Second_offset_variables <- NULL
eSVD_obj$param$fit_Second_tol <- NULL
eSVD_obj$param$nuisance_bool_library_includes_interept <- NULL
eSVD_obj$param$nuisance_bool_covariates_as_library <- NULL
eSVD_obj$param$posterior_alpha_max <- NULL
eSVD_obj$param$posterior_bool_adjust_covariates <- NULL
eSVD_obj$param$posterior_bool_covariates_as_library <- NULL
eSVD_obj$param$posterior_nuisance_lower_quantile <- NULL
eSVD_obj$param$test_case_individuals <- NULL
eSVD_obj$param$test_control_individuals <- NULL
eSVD_obj$fit_First <- NULL
eSVD_obj$fit_Second <- NULL

covariates <- cbind(1, init_res$covariates, init_res$offset_vec)
colnames(covariates)[1] <- "Intercept"
colnames(covariates)[ncol(covariates)] <- "Log_UMI"
eSVD_obj$covariates <- covariates
eSVD_obj$fit_Init$x_mat <- init_res$x_mat
eSVD_obj$fit_Init$y_mat <- init_res$y_mat
eSVD_obj$fit_Init$z_mat <- cbind(0, init_res$z_mat, init_res$offset_coef)
colnames(eSVD_obj$fit_Init$z_mat)[1] <- "Intercept"
colnames(eSVD_obj$fit_Init$z_mat)[ncol(eSVD_obj$fit_Init$z_mat)] <- "Log_UMI"

print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.01,
                             max_iter = 100,
                             offset_variables = c("Intercept", "Log_UMI"),
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

save(date_of_run, session_info, sns,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../out/Writeup11e/Writeup11e_sns_layer23_esvd9.RData")




