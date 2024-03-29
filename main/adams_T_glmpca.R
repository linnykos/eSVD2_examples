rm(list=ls())
library(Seurat)
library(glmpca)

load("../../../out/main/adams_T_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

var_features <- Seurat::VariableFeatures(adams)
mat <- adams[["RNA"]]@counts[var_features,]
covariate_dat <- adams@meta.data[,c("Disease_Identity", "Subject_Identity", "Gender",
                                    "Tobacco","percent.mt", "Age")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"Disease_Identity"] <- as.factor(covariate_df[,"Disease_Identity"])
covariate_df[,"Gender"] <- as.factor(covariate_df[,"Gender"])
covariate_df[,"Tobacco"] <- as.factor(covariate_df[,"Tobacco"])
covariate_df[,"Subject_Identity"] <- factor(covariate_df[,"Subject_Identity"], levels = names(sort(table(covariate_df[,"Subject_Identity"]), decreasing = F)))
covariates <- eSVD2:::format_covariates(dat = Matrix::t(mat),
                                        covariate_df = covariate_df,
                                        rescale_numeric_variables = c("percent.mt", "Age"))

covariates <- covariates[,which(!colnames(covariates) %in% "Intercept")]
set.seed(10)
while(TRUE){
  rank_res <- Matrix::rankMatrix(covariates)
  if(ncol(covariates) > rank_res){
    rm_num <- ncol(covariates)-rank_res
    col_idx <- grep("Subject_Identity", colnames(covariates))
    covariates <- covariates[,-sample(col_idx, rm_num),drop = T]
  } else {
    break()
  }
}

mat <- mat[Matrix::rowSums(mat)!=0,]
set.seed(10)
K <- 15
time_start <- Sys.time()
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "poi",
                             X = covariates,
                             ctl = list(verbose = T),
                             minibatch = "stochastic")
time_end <- Sys.time()
pred_mat <- glmpca:::predict.glmpca(glmpca_res)

save(date_of_run, session_info, adams,
     glmpca_res, pred_mat, covariates,
     time_start, time_end,
     file = "../../../out/main/adams_T_glmpca.RData")




