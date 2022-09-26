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
covariates <- covariates[,-grep("individual", colnames(covariates))]

downsample_values <- seq(0.95, 0.6, by = -0.05)
downsample_coef_list <- vector("list", length = length(downsample_values)+1)

p <- 1:ncol(obs_mat)
downsample_coef_list[[1]] <- t(sapply(1:p, function(j){
  y <- obs_mat[,j]
  df <- data.frame(cbind(y, covariates))
  colnames(df)[1] <- "y"

  glm_res <- MASS::glm.nb(y ~ . - 1, data = df)
  stats::coef(glm_res)
}))

downsample_values <- seq(0.95, 0.6, by = -0.05)
for(kk in 1:length(downsample_values)){
  downsample_value <- downsample_values[kk]
  print("==========================")
  print(paste0("Working on downsample: ", downsample_value))
  if("mat" %in% ls()) rm(list = "mat")
  gc()

  load(paste0("../../out/simulation/simulation_1_downsampled-", downsample_value, ".RData"))

  p <- ncol(mat)
  downsample_coef_list[[kk+1]] <- t(sapply(1:p, function(j){
    if(j %% floor(p/10) == 0) cat('*')
    y <- mat[,j]
    df <- data.frame(cbind(y, covariates))
    colnames(df)[1] <- "y"

    glm_res <- MASS::glm.nb(y ~ . - 1, data = df)
    stats::coef(glm_res)
  }))
}

save(date_of_run, session_info, seurat_obj,
     downsample_coef_list, downsample_values,
     file = paste0("../../out/simulation/simulation_1_nbreg_downsampled.RData"))

