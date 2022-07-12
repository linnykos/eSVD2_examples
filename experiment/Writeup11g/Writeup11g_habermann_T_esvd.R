rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/habermann_T_preprocessed.RData")
source("reparameterization.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

var_features <- Seurat::VariableFeatures(habermann)
mat <- Matrix::t(habermann[["RNA"]]@counts[var_features,])
covariate_dat <- habermann@meta.data[,c("Diagnosis", "Sample_Name", "Gender",
                                    "Tobacco","percent.mt", "Age")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"Diagnosis"] <- as.factor(covariate_df[,"Diagnosis"])
covariate_df[,"Gender"] <- as.factor(covariate_df[,"Gender"])
covariate_df[,"Tobacco"] <- droplevels(as.factor(covariate_df[,"Tobacco"]))
covariate_df[,"Sample_Name"] <- factor(covariate_df[,"Sample_Name"], levels = names(sort(table(covariate_df[,"Sample_Name"]), decreasing = F)))
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("percent.mt", "Age"))

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariates,
                                    case_control_variable = "Diagnosis_IPF",
                                    k = 15,
                                    lambda = 0.1,
                                    verbose = 1)
time_end1 <- Sys.time()

omitted_variables <- colnames(eSVD_obj$covariates)[grep("Sample_Name", colnames(eSVD_obj$covariates))]
eSVD_obj <- reparameterization2(eSVD_obj,
                                fit_name = "fit_Init",
                                omitted_variables = omitted_variables)

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

eSVD_obj <- reparameterization2(eSVD_obj,
                                fit_name = "fit_First",
                                omitted_variables = omitted_variables)

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

eSVD_obj <- reparameterization2(eSVD_obj,
                                fit_name = "fit_Second",
                                omitted_variables = omitted_variables)

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
     file = "../../../out/Writeup11g/Writeup11g_habermann_T_esvd.RData")



