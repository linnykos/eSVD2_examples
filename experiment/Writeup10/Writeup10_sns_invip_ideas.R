rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")
library(MAST)
library(ideas)

# see https://github.com/Sun-lab/ideas_pipeline/blob/main/simulation/step2_evaluate_methods.R
# following the analysis in https://github.com/himelmallick/BenchmarkSingleCell/blob/master/Library/run_MAST.R
# and https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])

# from https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1c_ideas.R
# https://github.com/Sun-lab/ideas_pipeline/blob/main/simulation/step2_evaluate_methods.R
dist1 <- ideas::ideas_dist(count_matrix = mat,
                           meta_cell = , meta_ind,
                           var_per_cell, var2test, var2adjust,
                           var2test_type, d_metric = d_metric,
                           fit_method = "nb")
