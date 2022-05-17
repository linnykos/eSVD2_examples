rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_layer4_processed2.RData")

library(Seurat)
library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])
covariate_dat <- sns@meta.data[,c("percent.mt", "individual", "region", "age", "sex",
                                  "RNA.Integrity.Number", "post.mortem.hours",
                                  "diagnosis", "Seqbatch")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"individual"] <- as.factor(covariate_df[,"individual"])
covariate_df[,"region"] <- as.factor(covariate_df[,"region"])
covariate_df[,"diagnosis"] <- factor(covariate_df[,"diagnosis"], levels = c("Control", "ASD"))
covariate_df[,"sex"] <- as.factor(covariate_df[,"sex"])
covariate_df[,"Seqbatch"] <- as.factor(covariate_df[,"Seqbatch"])
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        mixed_effect_variables = c("individual", "Seqbatch"))
case_control_variable <- "diagnosis_ASD"
mixed_effect_variables <- c(colnames(covariates)[grep("^individual", colnames(covariates))],
                            colnames(covariates)[grep("^Seqbatch", colnames(covariates))])

print("Initializing")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                   bool_intercept = T,
                                   covariates = covariates,
                                   case_control_variable = case_control_variable,
                                   k = 30,
                                   lambda = 0.1,
                                   mixed_effect_variables = mixed_effect_variables,
                                   offset_variables = "Log_UMI",
                                   verbose = 1)
time_end1 <- Sys.time()

# pval_thres <- min(10^(quantile(eSVD_obj$initial_Reg$log_pval, probs = 0.05)), 1e-5)
pval_thres <- 0.001
eSVD_obj <- eSVD2:::apply_initial_threshold(eSVD_obj = eSVD_obj,
                                            pval_thres = pval_thres,
                                            verbose = 1)

offset_variables <- setdiff(colnames(eSVD_obj$covariates), case_control_variable)
print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = offset_variables,
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_First",
                             fit_previous = "fit_Init")
time_end2 <- Sys.time()

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

print("Nuisance estimation")
time_start4 <- Sys.time()
eSVD_obj <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                      verbose = 1)
time_end4 <- Sys.time()

save(date_of_run, session_info, sns, covariate_df,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     file = "../../../../out/Writeup11c/Writeup11c_sns_layer4_esvd.RData")

eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj)
metadata <- sns@meta.data
metadata[,"individual"] <- as.factor(metadata[,"individual"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "individual",
                                           metadata = metadata,
                                           verbose = 2)
time_end5 <- Sys.time()

save(date_of_run, session_info, sns, covariate_df,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../out/Writeup11c/Writeup11c_sns_layer4_esvd.RData")



