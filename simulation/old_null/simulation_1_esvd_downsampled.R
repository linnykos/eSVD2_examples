rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../out/simulation/simulation_1.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##########

mat_tmp <- Matrix::t(seurat_obj[["RNA"]]@counts)
covariate_dat <- seurat_obj@meta.data[,c("cc", "age", "gender", "tobacco", "individual")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"cc"] <- as.factor(covariate_df[,"cc"])
covariate_df[,"tobacco"] <- factor(covariate_df[,"tobacco"], levels = names(sort(table(covariate_df[,"tobacco"]), decreasing = T)))
covariate_df[,"gender"] <- factor(covariate_df[,"gender"], levels = names(sort(table(covariate_df[,"gender"]), decreasing = T)))
covariate_df[,"individual"] <- factor(covariate_df[,"individual"], levels = names(sort(table(covariate_df[,"individual"]), decreasing = T)))
covariates <- eSVD2:::format_covariates(dat = mat_tmp,
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("age"))

downsample_values <- seq(0.9, 0.2, by = -0.1)
for(downsample_value in downsample_values){
  print("==========================")
  print("==========================")
  print(paste0("Working on downsample: ", downsample_value))
  if("mat" %in% ls()) rm(list = "mat")
  gc()

  load(paste0("../../out/simulation/simulation_1_downsampled-", downsample_value, ".RData"))
  set.seed(10)
  date_of_run <- Sys.time()
  session_info <- devtools::session_info()

  print("Initialization")
  time_start1 <- Sys.time()
  eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                      covariates = covariates[,-grep("individual", colnames(covariates))],
                                      case_control_variable = "cc_1",
                                      bool_intercept = T,
                                      k = 10,
                                      lambda = 0.1,
                                      metadata_case_control = covariates[,"cc_1"],
                                      metadata_individual = covariate_df[,"individual"],
                                      verbose = 1)
  time_end1 <- Sys.time()

  eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_Init",
    omitted_variables = "Log_UMI"
  )

  print("First fit")
  time_start2 <- Sys.time()
  eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                               l2pen = 0.1,
                               max_iter = 100,
                               offset_variables = setdiff(colnames(eSVD_obj$covariates), "cc_1"),
                               tol = 1e-6,
                               verbose = 1,
                               fit_name = "fit_First",
                               fit_previous = "fit_Init")
  time_end2 <- Sys.time()

  eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
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
    input_obj = eSVD_obj,
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
                                        bool_stabilize_underdispersion = T,
                                        library_min = 1,
                                        pseudocount = 1)

  time_start5 <- Sys.time()
  eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                             verbose = 1)
  time_end5 <- Sys.time()

  save(date_of_run, session_info, seurat_obj,
       eSVD_obj, downsample_value,
       file = paste0("../../out/simulation/simulation_1_esvd_downsampled-", downsample_value, ".RData"))
}

