rm(list=ls())
library(Seurat)
library(eSVD2)
library(glmpca)

load("../../../out/main/sns_layer23_glmpca.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
