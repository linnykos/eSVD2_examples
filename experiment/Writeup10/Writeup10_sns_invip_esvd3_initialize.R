rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")
source("initialization_lme4.R")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))

K <- min(30, round(min(dim(mat))*.5))
n <- nrow(mat)
p <- ncol(mat)

time_start1 <- Sys.time()
init_res <- initialize_lme4(mat,
                            k = K,
                            metadata = sns@meta.data,
                            pval_thres = 0.05,
                            tmp_path = "../../../../out/Writeup10/Writeup10_sns_invip_esvd3_tmp.RData",
                            verbose = 2)
time_end1 <- Sys.time()

save(date_of_run, session_info,
     sns, init_res, time_start1, time_end1,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_esvd3_initialize.RData")
