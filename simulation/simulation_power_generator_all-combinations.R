rm(list=ls())
library(Seurat)
library(eSVD2)

source("simulation_power_generator_functions.R")
set.seed(10)

gene_number_list <- list(
  gene100 = list(gene_num_mixed_membership = 100-14*4-28,
                 gene_num_null = 28,
                 gene_num_per_topic = 14),
  gene700 = list(gene_num_mixed_membership = 100,
                 gene_num_null = 200,
                 gene_num_per_topic = 100),
  gene5000 = list(gene_num_mixed_membership = 5000-710*4-1420,
                  gene_num_null = 1420,
                  gene_num_per_topic = 710)
)

individual_list <- list(
  size_factor1 = 1,
  size_factor2 = 2,
  size_factor4 = 4
)

gene_number = gene_number_list[[1]]
size_factor = individual_list[[1]]

for(ii in 1:length(gene_number_list)){
  print(paste0("Working on gene setting: ", ii))

  for(jj in 1:length(individual_list)){

    print(paste0("Working on size_factor setting: ", jj))

    gene_number <- gene_number_list[[ii]]
    size_factor <- individual_list[[jj]]

    set.seed(10)
    res <- form_simulation_data(
      gene_num_mixed_membership = gene_number$gene_num_mixed_membership,
      gene_num_null = gene_number$gene_num_null,
      gene_num_per_topic = gene_number$gene_num_per_topic,
      individual_cell_size_factor = size_factor,
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

    date_of_run <- Sys.time()
    session_info <- devtools::session_info()

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
         file = paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power_geneSetting",
                       ii, "_individualSetting", jj,
                       ".RData"))
  }
}

print("Done! :)")

