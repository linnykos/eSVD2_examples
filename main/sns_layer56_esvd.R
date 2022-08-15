rm(list=ls())
load("../../../out/main/sns_layer56_processed.RData")

library(Seurat)
library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##########

mat <- Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])
covariate_dat <- sns@meta.data[,c("percent.mt", "individual", "region", "age", "sex",
                                  "RNA.Integrity.Number", "post.mortem.hours",
                                  "diagnosis", "Seqbatch", "Capbatch")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"individual"] <- factor(covariate_df[,"individual"], levels = names(sort(table(covariate_df[,"individual"]), decreasing = T)))
covariate_df[,"region"] <- factor(covariate_df[,"region"], levels = names(sort(table(covariate_df[,"region"]), decreasing = T)))
covariate_df[,"diagnosis"] <- factor(covariate_df[,"diagnosis"], levels = c("Control", "ASD"))
covariate_df[,"sex"] <- factor(covariate_df[,"sex"], levels = names(sort(table(covariate_df[,"sex"]), decreasing = T)))
covariate_df[,"Seqbatch"] <- factor(covariate_df[,"Seqbatch"], levels = names(sort(table(covariate_df[,"Seqbatch"]), decreasing = T)))
covariate_df[,"Capbatch"] <- factor(covariate_df[,"Capbatch"], levels = names(sort(table(covariate_df[,"Capbatch"]), decreasing = T)))
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours"))

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariates,
                                    case_control_variable = "diagnosis_ASD",
                                    bool_intercept = T,
                                    k = 30,
                                    lambda = 0.1,
                                    verbose = 1)
time_end1 <- Sys.time()

omitted_variables <- colnames(eSVD_obj$covariates)[c(grep("individual", colnames(eSVD_obj$covariates)),
                                                     grep("Seqbatch", colnames(eSVD_obj$covariates)),
                                                     grep("Capbatch", colnames(eSVD_obj$covariates)))]
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
                             offset_variables = setdiff(colnames(eSVD_obj$covariates), "diagnosis_ASD"),
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
log_nuisance <- eSVD_obj2$fit_Second$nuisance_vec
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
                                      alpha_max = NULL,
                                      bool_covariates_as_library = T,
                                      library_min = 1e-4)
metadata <- sns@meta.data
metadata[,"individual"] <- as.factor(metadata[,"individual"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "individual",
                                           metadata = metadata,
                                           verbose = 1)
time_end5 <- Sys.time()

save(date_of_run, session_info, sns,
     eSVD_obj, log_nuisance,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start4b, time_end4b,
     time_start5, time_end5,
     file = "../../../out/main/sns_layer56_esvd.RData")



