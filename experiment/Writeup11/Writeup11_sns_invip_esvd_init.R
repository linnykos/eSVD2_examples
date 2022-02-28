rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
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

#####################

esvd_init <- eSVD2:::initialize_esvd(dat = dat,
                                     covariates = covariates,
                                     case_control_variable = "diagnosis_ASD",
                                     offset_variables = "Log_UMI",
                                     verbose = 2,
                                     tmp_path = "../../../../out/Writeup10/Writeup10_sns_invip_esvd_coef_tmp.RData")

save(dat, esvd_init,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_esvd_coef.RData")
