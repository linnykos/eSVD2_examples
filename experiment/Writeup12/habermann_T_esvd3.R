rm(list=ls())

library(Seurat)
load("../../../../out/main/habermann_T_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

var_features <- Seurat::VariableFeatures(habermann)
mat <- Matrix::t(habermann[["RNA"]]@counts[var_features,])
covariate_dat <- habermann@meta.data[,c("Diagnosis", "Sample_Name", "Gender",
                                        "Tobacco","percent.mt", "Age", "Ethnicity", "Sample_Source")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"Diagnosis"] <- as.factor(covariate_df[,"Diagnosis"])
covariate_df[,"Gender"] <- droplevels(factor(covariate_df[,"Gender"], levels = names(sort(table(covariate_df[,"Gender"]), decreasing = T))))
covariate_df[,"Tobacco"] <- droplevels(factor(covariate_df[,"Tobacco"], levels = names(sort(table(covariate_df[,"Tobacco"]), decreasing = T))))
covariate_df[,"Ethnicity"] <- droplevels(factor(covariate_df[,"Ethnicity"], levels = names(sort(table(covariate_df[,"Ethnicity"]), decreasing = T))))
# covariate_df[,"Sample_Source"] <- factor(covariate_df[,"Sample_Source"], levels = names(sort(table(covariate_df[,"Sample_Source"]), decreasing = T)))
covariate_df[,"Sample_Name"] <- factor(covariate_df[,"Sample_Name"], levels = names(sort(table(covariate_df[,"Sample_Name"]), decreasing = T)))
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("percent.mt", "Age"))

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariates[,-grep("Sample_Name", colnames(covariates))],
                                    case_control_variable = "Diagnosis_IPF",
                                    bool_intercept = T,
                                    k = 15,
                                    lambda = 0.1,
                                    metadata_case_control = covariates[,"Diagnosis_IPF"],
                                    metadata_individual = covariate_df[,"Sample_Name"],
                                    verbose = 1)
time_end1 <- Sys.time()

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  eSVD_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = "Log_UMI"
)

print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = setdiff(colnames(eSVD_obj$covariates), "Diagnosis_IPF"),
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_First",
                             fit_previous = "fit_Init")
time_end2 <- Sys.time()

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  eSVD_obj = eSVD_obj,
  fit_name = "fit_First",
  omitted_variables = "Log_UMI"
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
  omitted_variables = numeric(0)
)

print("Nuisance estimation")
time_start4 <- Sys.time()
eSVD_obj <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                      bool_covariates_as_library = T,
                                      bool_library_includes_interept = T,
                                      bool_use_log = F,
                                      verbose = 1)
time_end4 <- Sys.time()

eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      alpha_max = NULL,
                                      bool_covariates_as_library = T,
                                      bool_stabilize_underdispersion = F,
                                      library_min = 1,
                                      pseudocount = 1)

time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           verbose = 1)
time_end5 <- Sys.time()

save(date_of_run, session_info, habermann,
     eSVD_obj, covariate_df,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../out/Writeup12/habermann_T_esvd3.RData")



