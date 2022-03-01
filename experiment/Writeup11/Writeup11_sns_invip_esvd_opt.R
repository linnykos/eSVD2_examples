rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_esvd_coef.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
