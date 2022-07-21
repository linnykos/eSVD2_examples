rm(list=ls())
load("../../../out/main/regevEpi_entprog_preprocessed.RData")

library(Seurat)
library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(regevEpi))
keep_vec[which(regevEpi$Sample_Health == "Inflamed")] <- 0
regevEpi$keep <- keep_vec
regevEpi <- subset(regevEpi, keep == 1)
regevEpi[["percent.mt"]] <- Seurat::PercentageFeatureSet(regevEpi, pattern = "^MT-")

##########

mat <- Matrix::t(regevEpi[["RNA"]]@counts[regevEpi[["RNA"]]@var.features,])
covariate_dat <- regevEpi@meta.data[,c("percent.mt", "Sample", "Subject_Disease", "Subject_Gender",
                                       "Subject_Location", "Subject_Smoking")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"Subject_Disease"] <- factor(covariate_df[,"Subject_Disease"], levels = c("HC", "Colitis"))
covariate_df[,"Sample"] <- factor(covariate_df[,"Sample"], levels = names(sort(table(covariate_df[,"Sample"]), decreasing = T)))
covariate_df[,"Subject_Gender"] <- factor(covariate_df[,"Subject_Gender"], levels = names(sort(table(covariate_df[,"Subject_Gender"]), decreasing = T)))
covariate_df[,"Subject_Location"] <- factor(covariate_df[,"Subject_Location"], levels = names(sort(table(covariate_df[,"Subject_Location"]), decreasing = T)))
covariate_df[,"Subject_Smoking"] <- factor(covariate_df[,"Subject_Smoking"], levels = names(sort(table(covariate_df[,"Subject_Smoking"]), decreasing = T)))

covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("percent.mt"))

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariates,
                                    case_control_variable = "Subject_Disease_Colitis",
                                    bool_intercept = T,
                                    k = 15,
                                    lambda = 0.1,
                                    verbose = 1)
time_end1 <- Sys.time()

omitted_variables <- colnames(eSVD_obj$covariates)[c(grep("^Sample_N", colnames(eSVD_obj$covariates)),
                                                     grep("^Subject_Location", colnames(eSVD_obj$covariates)))]
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
                             offset_variables = setdiff(colnames(eSVD_obj$covariates), "Subject_Disease_Colitis"),
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
metadata <- regevEpi@meta.data
metadata[,"Sample"] <- as.factor(metadata[,"Sample"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Sample",
                                           metadata = metadata,
                                           verbose = 1)
time_end5 <- Sys.time()

save(date_of_run, session_info, regevEpi,
     eSVD_obj, log_nuisance,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start4b, time_end4b,
     time_start5, time_end5,
     file = "../../../out/main/regevEpi_entprog-noninflamed_esvd.RData")



