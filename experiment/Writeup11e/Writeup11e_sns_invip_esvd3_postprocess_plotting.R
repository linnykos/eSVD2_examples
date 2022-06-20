rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/Writeup11e/Writeup11e_sns_invip_esvd3.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "IN-VIP"),]
de_gene_specific <- tmp[,"Gene name"]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
de_idx <- which(gene_names %in% de_gene_specific)
cycling_idx <- setdiff(cycling_idx, de_idx)

tab <- table(sns$individual, sns$diagnosis)
indiv_cases <- rownames(tab)[which(tab[,"ASD"] != 0)]
indiv_controls <- rownames(tab)[which(tab[,"Control"] != 0)]
indiv_vec <- factor(as.character(sns$individual))

# round(apply(eSVD_obj$fit_Second$z_mat, 2, quantile), 2)

#########################################

png(paste0("../../../../out/fig/Writeup11e/Writeup11e_sns_invip_esvd3_diagnostic_gene.png"),
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
                  what_1 = "diagnosis_ASD",
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


png(paste0("../../../../out/fig/Writeup11e/Writeup11e_sns_invip_esvd3_diagnostic_gene2.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "diagnosis_ASD",
                  what_2 = "diagnosis_ASD",
                  which_fit_1 = "fit_Init",
                  which_fit_2 = "fit_First",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  xlab = "Initial diagnosis_ASD",
                  ylab = "First diagnosis_ASD",
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "diagnosis_ASD",
                  what_2 = "diagnosis_ASD",
                  which_fit_1 = "fit_First",
                  which_fit_2 = "fit_Second",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  xlab = "First diagnosis_ASD",
                  ylab = "Second diagnosis_ASD",
                  color_palette = c(3,2))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "diagnosis_ASD",
                  what_2 = "diagnosis_ASD",
                  which_fit_1 = "fit_Init",
                  which_fit_2 = "fit_Second",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  xlab = "Initial diagnosis_ASD",
                  ylab = "Second diagnosis_ASD",
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

  png(paste0("../../../../out/fig/Writeup11e/Writeup11e_sns_invip_esvd3_diagnostic_cell-cycling-",
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
                    fit_included_covariates_1 = c("diagnosis_ASD"),
                    fit_included_covariates_2 = c("diagnosis_ASD", colnames(eSVD_obj$covariates)[grep("individual", colnames(eSVD_obj$covariates))]),
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec
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

  png(paste0("../../../../out/fig/Writeup11e/Writeup11e_sns_invip_esvd3_diagnostic_cell-de-",
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
                    fit_included_covariates_1 = c("diagnosis_ASD"),
                    fit_included_covariates_2 = c("diagnosis_ASD", colnames(eSVD_obj$covariates)[grep("individual", colnames(eSVD_obj$covariates))]),
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec
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
