rm(list=ls())
set.seed(10)
library(Seurat)
library(eSVD2)
source("../eSVD2_examples/simulation/simulation_power_generator.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

res <- form_simulation_data(
  gene_num_mixed_membership = 100,
  gene_num_null = 200,
  gene_num_per_topic = 100,
  individual_cell_size_factor = 1,
  num_topics = 4
)

seurat_obj <- res$seurat_obj
case_individuals <- res$case_individuals
control_individuals <- res$control_individuals
covariates <- res$covariates
gene_labeling <- res$gene_labeling
gene_labeling2 <- res$gene_labeling2
gene_library_vec <- res$gene_library_vec
individual_vec <- res$individual_vec
nuisance_vec <- res$nuisance_vec
nat_mat <- res$nat_mat
obs_mat <- res$obs_mat
true_fdr_vec <- res$true_fdr_vec
true_logpvalue_vec <- res$true_logpvalue_vec
true_null_mean <- res$true_null_mean
true_null_sd <- res$true_null_sd
true_teststat_vec <- res$true_teststat_vec
x_mat <- res$x_mat
y_mat <- res$y_mat
z_mat <- res$z_mat

save(seurat_obj,
     case_individuals,
     control_individuals,
     covariates,
     gene_labeling,
     gene_labeling2,
     gene_library_vec,
     individual_vec,
     nuisance_vec,
     nat_mat,
     obs_mat,
     true_fdr_vec,
     true_logpvalue_vec,
     true_null_mean,
     true_null_sd,
     true_teststat_vec,
     x_mat,
     y_mat,
     z_mat,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_power.RData")

Seurat::DimPlot(seurat_obj, reduction = "isomap", group.by = "individual")
Seurat::DimPlot(seurat_obj, reduction = "isomap", group.by = "cc")
Seurat::DimPlot(seurat_obj, reduction = "isomap", group.by = "gender")
Seurat::FeaturePlot(seurat_obj, reduction = "isomap", features = "age")
