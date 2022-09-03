rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)

load("../../../../out/Writeup12/adams_T_esvd.RData")
date_of_run

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

################################
