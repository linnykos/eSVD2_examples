rm(list=ls())
load("../../../../out/Writeup11c/Writeup11c_adams_T_esvd.RData")

library(Seurat)
library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
eSVD_obj$param$init_library_size_variable <- "Log_UMI"

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

tab <- table(adams$Subject_Identity, adams$Disease_Identity)
indiv_cases <- rownames(tab)[which(tab[,"IPF"] != 0)]
indiv_controls <- rownames(tab)[which(tab[,"Control"] != 0)]
indiv_vec <- factor(as.character(adams$Subject_Identity))

round(apply(eSVD_obj$fit_Second$z_mat, 2, quantile), 2)

#########################################

png(paste0("../../../../out/fig/Writeup11c/Writeup11c_adams_T_gene.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(2,3))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "Log_UMI",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(2,3))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "Disease_Identity_IPF",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(2,3))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "Intercept",
                  what_2 = "teststat",
                  gene_list = list(gene_names[cycling_idx],
                                   gene_names[de_idx]),
                  color_palette = c(2,3))

graphics.off()

#########################################

# select 25 cell-cycling genes
set.seed(10)
cycling_idx_subset <- sample(cycling_idx, 25)

for(k in 1:length(cycling_idx_subset)){
  print(k)

  idx <- cycling_idx_subset[k]
  tmp <- eSVD_obj$dat[,idx]

  png(paste0("../../../../out/fig/Writeup11c/Writeup11c_adams_T_cell-cycling-",
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
                    what_2 = "relative_expression",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Nuisance = ", round(eSVD_obj$fit_Second$nuisance_vec[idx], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean",
                    what_2 = "posterior_mean_nonadjusted",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec)

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean",
                    what_2 = "posterior_variance",
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

  png(paste0("../../../../out/fig/Writeup11c/Writeup11c_adams_T_cell-de-",
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
                    what_2 = "relative_expression",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Nuisance = ", round(eSVD_obj$fit_Second$nuisance_vec[idx], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean",
                    what_2 = "posterior_mean_nonadjusted",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec)

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean",
                    what_2 = "posterior_variance",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Test statistic = ", round(eSVD_obj$teststat_vec[idx], 2))
  )

  graphics.off()
}
