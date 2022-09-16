rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##########

mat_tmp <- sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]
covariate_dat <- sns@meta.data[,c("percent.mt", "region", "age", "sex",
                                  "RNA.Integrity.Number", "post.mortem.hours",
                                  "diagnosis", "Seqbatch", "Capbatch")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"region"] <- factor(covariate_df[,"region"], levels = names(sort(table(covariate_df[,"region"]), decreasing = T)))
covariate_df[,"diagnosis"] <- factor(covariate_df[,"diagnosis"], levels = c("Control", "ASD"))
covariate_df[,"sex"] <- factor(covariate_df[,"sex"], levels = names(sort(table(covariate_df[,"sex"]), decreasing = T)))
covariate_df[,"Seqbatch"] <- factor(covariate_df[,"Seqbatch"], levels = names(sort(table(covariate_df[,"Seqbatch"]), decreasing = T)))
covariate_df[,"Capbatch"] <- factor(covariate_df[,"Capbatch"], levels = names(sort(table(covariate_df[,"Capbatch"]), decreasing = T)))
covariates <- eSVD2:::format_covariates(dat = Matrix::t(mat_tmp),
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours"))
covariates <- covariates[,which(!colnames(covariates) %in% "Intercept")]

downsample_values <- seq(0.95, 0.6, by = -0.05)
for(downsample_value in downsample_values){
  print(paste0("Working on downsample: ", downsample_value))
  if("mat" %in% ls()) rm(list = "mat")
  gc()

  load(paste0("../../../out/main/sns_layer23_processed_downsampled-", downsample_value, ".RData"))
  set.seed(10)
  date_of_run <- Sys.time()
  session_info <- devtools::session_info()

  mat <- Matrix::t(mat)
  mat <- mat[Matrix::rowSums(mat)!=0,]
  set.seed(10)
  K <- 30
  time_start <- Sys.time()
  glmpca_res <- glmpca::glmpca(mat, L = K, fam = "poi",
                               X = covariates,
                               ctl = list(verbose = T),
                               minibatch = "stochastic")
  time_end <- Sys.time()
  pred_mat <- glmpca:::predict.glmpca(glmpca_res)

  save(date_of_run, session_info, sns,
       glmpca_res, pred_mat, downsample_value,
       time_start, time_end,
       file = paste0("../../../out/main/sns_layer23_glmpca_downsampled-", downsample_value, ".RData"))
}




