rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/regevImm_macro_preprocessed.RData")
# table(regevImm$Subject, regevImm$Sample_Health)
# table(regevImm$Subject, regevImm$Subject_Disease)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(regevImm))
keep_vec[which(regevImm$Sample_Health == "Inflamed")] <- 0
regevImm$keep <- keep_vec
regevImm <- subset(regevImm, keep == 1)

# take only half of the healthy subjects
tab <- table(regevImm$Subject, regevImm$Subject_Disease)
healthy_subj <- rownames(tab[tab[,"HC"] != 0,])
set.seed(10)
split1 <- sample(healthy_subj, size = round(length(healthy_subj)/2), replace = F)
split2 <- setdiff(healthy_subj, split1)
keep_vec <- rep(1, ncol(regevImm))
# Non-inflamed analysis uses split1, Inflamed uses split2
if(any(regevImm$Sample_Health == "Non-inflamed")){
  keep_vec[which(regevImm$Subject %in% split2)] <- 0
} else {
  keep_vec[which(regevImm$Subject %in% split1)] <- 0
}
regevImm$keep <- keep_vec
regevImm <- subset(regevImm, keep == 1)

regevImm[["percent.mt"]] <- Seurat::PercentageFeatureSet(regevImm, pattern = "^MT-")

##########

mat <- Matrix::t(regevImm[["RNA"]]@counts[regevImm[["RNA"]]@var.features,])
covariate_dat <- regevImm@meta.data[,c("percent.mt", "Sample", "Subject_Disease", "Subject_Gender",
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
                                    metadata_case_control = covariates[,"Subject_Disease_Colitis"],
                                    metadata_individual = covariate_df[,"Sample"],
                                    verbose = 1)
time_end1 <- Sys.time()

omitted_variables <- colnames(eSVD_obj$covariates)[c(grep("^Sample_N", colnames(eSVD_obj$covariates)),
                                                     grep("^Subject_Location", colnames(eSVD_obj$covariates)))]
eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
  input_obj = eSVD_obj,
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
  input_obj = eSVD_obj,
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
  input_obj = eSVD_obj,
  fit_name = "fit_Second",
  omitted_variables = omitted_variables
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

save(date_of_run, session_info, regevImm,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     time_start5, time_end5,
     file = "../../../../out/Writeup13/regevImm_macro-noninflamed_esvd.RData")



