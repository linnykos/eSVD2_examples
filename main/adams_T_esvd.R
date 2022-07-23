rm(list=ls())

library(Seurat)
load("../../../out/main/adams_T_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

var_features <- Seurat::VariableFeatures(adams)
mat <- Matrix::t(adams[["RNA"]]@counts[var_features,])
covariate_dat <- adams@meta.data[,c("Disease_Identity", "Subject_Identity", "Gender",
                                  "Tobacco","percent.mt", "Age")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"Disease_Identity"] <- as.factor(covariate_df[,"Disease_Identity"])
covariate_df[,"Gender"] <- as.factor(covariate_df[,"Gender"])
covariate_df[,"Tobacco"] <- as.factor(covariate_df[,"Tobacco"])
covariate_df[,"Subject_Identity"] <- factor(covariate_df[,"Subject_Identity"], levels = names(sort(table(covariate_df[,"Subject_Identity"]), decreasing = F)))
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("percent.mt", "Age"))

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariates,
                                    case_control_variable = "Disease_Identity_IPF",
                                    k = 15,
                                    lambda = 0.1,
                                    verbose = 1)
time_end1 <- Sys.time()


omitted_variables <- colnames(eSVD_obj$covariates)[grep("Subject_Identity", colnames(eSVD_obj$covariates))]
eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  eSVD_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = c("Log_UMI", omitted_variables)
)

print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = setdiff(colnames(eSVD_obj$covariates), "Disease_Identity_IPF"),
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_First",
                             fit_previous = "fit_Init")
time_end2 <- Sys.time()

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  eSVD_obj = eSVD_obj,
  fit_name = "fit_First",
  omitted_variables = c("Log_UMI", omitted_variables)
)

print("Second fit")
time_start3 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = NULL,
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_Second",
                             fit_previous = "fit_First")
time_end3 <- Sys.time()

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  eSVD_obj = eSVD_obj,
  fit_name = "fit_Second",
  omitted_variables = omitted_variables
)

print("Nuisance estimation")
time_start4 <- Sys.time()
eSVD_obj2 <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                       bool_covariates_as_library = T,
                                       bool_library_includes_interept = T,
                                       bool_use_log = T,
                                       verbose = 1)
log_nuisance <- eSVD_obj$fit_Second$nuisance_vec
time_end4 <- Sys.time()

time_start4b <- Sys.time()
eSVD_obj <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                      bool_covariates_as_library = T,
                                      bool_library_includes_interept = T,
                                      bool_use_log = F,
                                      verbose = 1)
time_end4b <- Sys.time()

eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      bool_covariates_as_library = T)
metadata <- adams@meta.data
metadata[,"Subject_Identity"] <- as.factor(metadata[,"Subject_Identity"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Subject_Identity",
                                           metadata = metadata,
                                           verbose = 1)
time_end5 <- Sys.time()

save(date_of_run, session_info, adams,
     eSVD_obj, log_nuisance,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start4b, time_end4b,
     time_start5, time_end5,
     file = "../../../out/main/adams_T_esvd.RData")



