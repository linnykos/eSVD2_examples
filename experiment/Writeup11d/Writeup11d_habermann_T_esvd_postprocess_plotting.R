rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/Writeup11d/Writeup11d_habermann_T_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# eSVD_obj$teststat_vec <- NULL
# eSVD_obj$fit_Second$posterior_mean_mat <- NULL
# eSVD_obj$fit_Second$posterior_var_mat <- NULL
#
# eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
#                                       bool_adjust_covariates = F,
#                                       bool_covariates_as_library = T)

metadata <- habermann@meta.data
metadata[,"Sample_Name"] <- as.factor(metadata[,"Sample_Name"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Sample_Name",
                                           metadata = metadata,
                                           verbose = 2)
time_end5 <- Sys.time()

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("B", "Macrophage", "Macrophage Alveolar", "NK"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("Macrophages_disease_vs_control_.csv", "Monocytes_disease_vs_control_.csv",
              "B_Cells_disease_vs_control_.csv", "NK_Cells_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
de_idx <- which(gene_names %in% de_genes)

tab <- table(habermann$Sample_Name, habermann$Diagnosis)
indiv_cases <- rownames(tab)[which(tab[,"IPF"] != 0)]
indiv_controls <- rownames(tab)[which(tab[,"Control"] != 0)]
indiv_vec <- factor(as.character(habermann$Sample_Name))

round(apply(eSVD_obj$fit_Second$z_mat, 2, quantile), 2)

#########################################

png(paste0("../../../../out/fig/Writeup11d/Writeup11d_habermann_T_diagnostic_gene.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "sparsity",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "Diagnosis_IPF",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "sparsity",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))

graphics.off()


#########################################

# select 25 cell-cycling genes
set.seed(10)
cycling_idx_subset <- sample(cycling_idx, 25)

for(k in 1:length(cycling_idx_subset)){
  print(k)

  idx <- cycling_idx_subset[k]
  tmp <- eSVD_obj$dat[,idx]

  png(paste0("../../../../out/fig/Writeup11d/Writeup11d_habermann_T_diagnostic_cell-cycling-",
             gene_names[idx], ".png"),
      height = 2500, width = 2500,
      units = "px", res = 300)
  par(mfrow = c(2,2), mar = c(4,4,4,0.5))

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "fit",
                    what_2 = "dat",
                    bool_jitter_y = T,
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0(gene_names[idx], ": Zero% = ", round(100*length(which(tmp == 0))/length(tmp)))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean_nonadjusted",
                    what_2 = "adjusted_relative_expression",
                    bool_jitter_y = T,
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Nuisance = ", round(eSVD_obj$fit_Second$nuisance_vec[idx], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "fit",
                    what_2 = "fit",
                    fit_included_covariates_1 = c("Diagnosis_IPF"),
                    fit_included_covariates_2 = c("Diagnosis_IPF", colnames(eSVD_obj$covariates)[grep("Sample_Name", colnames(eSVD_obj$covariates))]),
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("CC coefficient: ", round(eSVD_obj$fit_Second$z_mat[idx,"Diagnosis_IPF"], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean",
                    what_2 = "posterior_sd",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Test statistic = ", round(eSVD_obj$teststat_vec[idx], 2))
  )

  graphics.off()
}


##############################

# select 25 signal genes
set.seed(10)
de_idx_subset <- sample(de_idx, 25)

for(k in 1:length(de_idx_subset)){
  print(k)

  idx <- de_idx_subset[k]
  tmp <- eSVD_obj$dat[,idx]

  png(paste0("../../../../out/fig/Writeup11d/Writeup11d_habermann_T_diagnostic_cell-de-",
             gene_names[idx], ".png"),
      height = 2500, width = 2500,
      units = "px", res = 300)
  par(mfrow = c(2,2), mar = c(4,4,4,0.5))


  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "fit",
                    what_2 = "dat",
                    bool_jitter_y = T,
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0(gene_names[idx], ": Zero% = ", round(100*length(which(tmp == 0))/length(tmp)))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean_nonadjusted",
                    what_2 = "adjusted_relative_expression",
                    bool_jitter_y = T,
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Nuisance = ", round(eSVD_obj$fit_Second$nuisance_vec[idx], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "fit",
                    what_2 = "fit",
                    fit_included_covariates_1 = c("Diagnosis_IPF"),
                    fit_included_covariates_2 = c("Diagnosis_IPF", colnames(eSVD_obj$covariates)[grep("Sample_Name", colnames(eSVD_obj$covariates))]),
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("CC coefficient: ", round(eSVD_obj$fit_Second$z_mat[idx,"Diagnosis_IPF"], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean",
                    what_2 = "posterior_sd",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Test statistic = ", round(eSVD_obj$teststat_vec[idx], 2))
  )

  graphics.off()
}

############################

png(paste0("../../../../out/fig/Writeup11d/Writeup11d_habermann_T_diagnostic_gene_histogram.png"),
    height = 1500, width = 2500,
    units = "px", res = 300)
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))
graphics.off()

eSVD_obj_habermann <- eSVD_obj

load("../../../../out/Writeup11d/Writeup11d_adams_T_esvd.RData")

eSVD_obj$teststat_vec <- NULL
metadata <- adams@meta.data
metadata[,"Subject_Identity"] <- as.factor(metadata[,"Subject_Identity"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Subject_Identity",
                                           metadata = metadata,
                                           verbose = 2)
time_end5 <- Sys.time()
eSVD_obj_adams <- eSVD_obj

png(paste0("../../../../out/fig/Writeup11d/Writeup11d_adams_T_diagnostic_gene_histogram.png"),
    height = 1500, width = 2500,
    units = "px", res = 300)
eSVD2:::gene_plot(eSVD_obj_adams,
                  what_1 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(3,2))
graphics.off()

png(paste0("../../../../out/fig/Writeup11d/Writeup11d_adams_habermann_T_teststatistic.png"),
    height = 1500, width = 3500,
    units = "px", res = 300)
par(mfrow = c(1,3))
plot(x = eSVD_obj_adams$teststat_vec,
     y = eSVD_obj_habermann$teststat_vec,
     xlab = "Adams", ylab = "Habermann", pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
     main = "All genes", asp = T)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[cycling_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[cycling_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 3,
     main = "Cycling genes", asp = T)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)

plot(x = eSVD_obj_adams$teststat_vec[gene_names[de_idx]],
     y = eSVD_obj_habermann$teststat_vec[gene_names[de_idx]],
     xlab = "Adams", ylab = "Habermann", pch = 16, col = 2,
     main = "DE genes", asp = T)
lines(c(-100,100), rep(0,2), col = 2, lwd = 2, lty = 2)
lines(rep(0,2), c(-100,100), col = 2, lwd = 2, lty = 2)
graphics.off()
