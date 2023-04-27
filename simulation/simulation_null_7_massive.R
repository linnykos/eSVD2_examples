rm(list=ls())
library(Seurat)
library(eSVD2)

source("simulation_null_7_massive_functions.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

trials <- 1000
result_list <- rep(NA, trials)

for(i in 1:trials){
  print(i)
  dat <- simulation_generate_data(seed = i)
  result_list[[i]] <- simulation_run_esvd(
    covariate = dat$covariate,
    seurat_obj = dat$seurat_obj
  )

  save(result_list, date_of_run, session_info,
       file = "../../../out/simulation/simulation_null_7_tmp.RData")
}
