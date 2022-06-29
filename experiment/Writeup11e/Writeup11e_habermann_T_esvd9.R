rm(list=ls())

load("../../../../out/Writeup11d/Writeup11d_habermann_T_esvd.RData")
source("initialization.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

var_features <- colnames(eSVD_obj$dat)
mat <- Matrix::t(habermann[["RNA"]]@counts[var_features,])
categorical_var <- c("Diagnosis", "Sample_Name", "Gender", "Tobacco")
numerical_var <- c("percent.mt", "Age", "nFeature_RNA")
n <- nrow(mat)
metadata <- habermann@meta.data
covariates <- as.matrix(metadata[,numerical_var])
covariates <- cbind(1, log(Matrix::rowSums(mat)), covariates)
colnames(covariates)[1:2] <- c("Intercept", "Log_UMI")

for(variable in categorical_var){
  vec <- metadata[,variable]
  uniq_level <- unique(vec)
  for(i in uniq_level[-1]){
    tmp <- rep(0, n)
    tmp[which(vec == i)] <- 1

    var_name <- paste0(variable, "_", i)
    covariates <- cbind(covariates, tmp)
    colnames(covariates)[ncol(covariates)] <- var_name
  }
}

dat <- as.matrix(Matrix::t(habermann[["RNA"]]@counts[var_features,]))
init_res <- initialize_esvd2(dat = dat,
                             k = 30,
                             covariates = covariates,
                             rescale_vars = c("Age", "nFeature_RNA", "percent.mt"),
                             batches_names = NULL,
                             bool_intercept = F,
                             subject_names = "Sample_Name_",
                             verbose = 2)

# replace eSVD_obj with the appropriate things
eSVD_obj$param$init_bool_intercept <- F
eSVD_obj$param$init_k <- 30

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
metadata <- habermann@meta.data
metadata[,"Sample_Name"] <- as.factor(metadata[,"Sample_Name"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Sample_Name",
                                           metadata = metadata,
                                           verbose = 1)
time_end5 <- Sys.time()

save(date_of_run, session_info, habermann,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../out/Writeup11e/Writeup11e_habermann_T_esvd9.RData")





