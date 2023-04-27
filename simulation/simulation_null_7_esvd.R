rm(list=ls())
library(Seurat)
library(eSVD2)

load("../eSVD2_examples/simulation/simulation_null_7.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##########

mat <- Matrix::t(seurat_obj[["RNA"]]@counts)
covariate_dat <- seurat_obj@meta.data[,c("CC", "Sex", "Age", "Individual")]
covariate_df <- data.frame(covariate_dat)
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df)
covariate_df[,"Individual"] <- as.factor(covariate_df[,"Individual"])

print("Initialization")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    covariates = covariate[,c("Intercept", "CC", "Log_UMI", "Sex", "Age")],
                                    case_control_variable = "CC",
                                    bool_intercept = T,
                                    k = 2,
                                    lambda = 0.1,
                                    metadata_case_control = covariate[,"CC"],
                                    metadata_individual = covariate_df[,"Individual"],
                                    verbose = 1)
time_end1 <- Sys.time()

eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
  fit_name = "fit_Init",
  omitted_variables = "Log_UMI"
)

print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(
  input_obj = eSVD_obj,
  l2pen = 0.1,
  max_iter = 100,
  offset_variables = setdiff(colnames(eSVD_obj$covariates), "CC"),
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
eSVD_obj <- eSVD2:::opt_esvd(
  input_obj = eSVD_obj,
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
  omitted_variables = NULL
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

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
teststat_vec <- eSVD_obj$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

multtest_res <- eSVD2:::multtest(gaussian_teststat)

multtest_res$method
c(multtest_res$null_mean, multtest_res$null_sd)
# hist(gaussian_teststat[-c(1:10)])
# quantile(gaussian_teststat[1:10])
plot(sort(multtest_res$pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(multtest_res$pvalue_vec[-c(1:10)])), asp = T)
lines(c(0,1), c(0,1), col = 2, lty = 2)

save(date_of_run, session_info,
     seurat_obj, eSVD_obj,
     gaussian_teststat, multtest_res,
     covariate,
     df,
     gene_plot_idx,
     x_mat,
     y_mat,
     y_block_assignment,
     z_mat,
     file = "../eSVD2_examples/simulation/simulation_null_7_esvd.RData")
