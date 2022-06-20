rm(list=ls())
library(Seurat)

load("../../../../out/Writeup11e/Writeup11e_sns_invip_esvd.RData")
# load("../../../../out/Writeup10/Writeup10_sns_invip_esvd2.RData")
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")
source("initialization.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# doing the initialization as in Writeup10, but now with an intercept
dat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
init_res <- initialize_esvd2(dat = dat,
                             k = 50,
                             covariates = covariates,
                             verbose = 1)

# replace eSVD_obj with the appropriate things
eSVD_obj$param$init_k <- 50
eSVD_obj$fit_First <- NULL
eSVD_obj$fit_Second <- NULL
covariates <- cbind(init_res$covariates, init_res$offset_vec)
colnames(covariates)[ncol(covariates)] <- "Log_UMI"
eSVD_obj$covariates <- covariates
eSVD_obj$fit_Init$x_mat <- init_res$x_mat
eSVD_obj$fit_Init$y_mat <- init_res$y_mat
eSVD_obj$fit_Init$z_mat <- cbind(init_res$b_mat, 1)
colnames(eSVD_obj$fit_Init$z_mat)[ncol(eSVD_obj$fit_Init$z_mat)] <- "Log_UMI"


