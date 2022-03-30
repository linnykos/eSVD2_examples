rm(list=ls())
library(Seurat)
load("../../../../out/Writeup11/Writeup11_adams_ciliated_esvd.RData")
source("../Writeup10/multiple_testing.R")

mat <- as.matrix(Matrix::t(adams[["RNA"]]@counts[adams[["RNA"]]@var.features,]))
case_control_variable <- "Disease_Identity_IPF"
offset_var <- setdiff(colnames(esvd_init$covariates), case_control_variable)

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,case_control_variable,drop = F],
                       esvd_res_full$b_mat[,case_control_variable,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,offset_var],
  esvd_res_full$b_mat[,offset_var]
))

nuisance_vec <- pmax(nuisance_vec, quantile(nuisance_vec, probs = 0.01))
Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
Alpha <- pmin(Alpha, 50)
AplusAlpha <- mat + Alpha
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2

#################

metadata <- adams@meta.data
case_individuals <- unique(metadata[which(metadata$Disease_Identity == "IPF"),"Subject_Identity"])
control_individuals <- unique(metadata[which(metadata$Disease_Identity == "Control"),"Subject_Identity"])
case_idx <- which(metadata[,"Disease_Identity"] == "IPF")
control_idx <- which(metadata[,"Disease_Identity"] == "Control")

individual_stats <- lapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  # next find the cells, then compute one gaussian per individual
  case_gaussians <- sapply(case_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$Subject_Identity == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$Subject_Identity == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  list(case_gaussians = case_gaussians,
       control_gaussians = control_gaussians)
})

# see https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
group_stats <- lapply(1:length(individual_stats), function(j){
  case_gaussians <- individual_stats[[j]]$case_gaussians
  control_gaussians <- individual_stats[[j]]$control_gaussians

  case_gaussian <- list(mean_val = mean(case_gaussians[1,]),
                        var_val = mean(case_gaussians[2,]) + mean(case_gaussians[1,]^2) - (mean(case_gaussians[1,]))^2,
                        n = ncol(case_gaussians))
  control_gaussian <- list(mean_val = mean(control_gaussians[1,]),
                           var_val = mean(control_gaussians[2,]) + mean(control_gaussians[1,]^2) - (mean(control_gaussians[1,]))^2,
                           n = ncol(control_gaussians))

  list(case_gaussian = case_gaussian,
       control_gaussian = control_gaussian)
})

teststat_vec <- sapply(1:length(group_stats), function(j){
  case_gaussian <- group_stats[[j]]$case_gaussian
  control_gaussian <- group_stats[[j]]$control_gaussian

  n1 <- control_gaussian$n; n2 <- case_gaussian$n
  mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
  cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$var_val

  combined_cov <- cov1/n1 + cov2/n2
  (mean2 - mean1)/sqrt(combined_cov)
})

##########

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "Ciliated")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("AT1", "AT2", "Basal", "Club", "Goblet", "Mesothelial"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/Ciliated_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("AT1_disease_vs_control_.csv", "AT2_disease_vs_control_.csv",
              "Basal_disease_vs_control_.csv", "Differentiating_Ciliated_disease_vs_control_.csv",
              "KRT5-KRT17+_disease_vs_control_.csv", "MUC5AC+_High_disease_vs_control_.csv",
              "MUC5B+_disease_vs_control_.csv", "Proliferating_Epithelial_Cells_disease_vs_control_.csv",
              "SCGB3A2+_disease_vs_control_.csv", "SCGB3A2+_SCGB1A1+_disease_vs_control_.csv",
              "Transitional_AT2_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
other_genes <- unique(c(adams_df_genes_others, habermann_df_genes_others))

# hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(colnames(mat) %in% cycling_genes)
de_idx <- which(colnames(mat) %in% de_genes)
other_idx <- which(colnames(mat) %in% other_genes)

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

teststat_vec <- pmax(pmin(teststat_vec, 30), -30)
max_val <- max(abs(teststat_vec))
png("../../../../out/fig/Writeup11/adams_ciliated_esvd_teststat_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05
hist(teststat_vec, breaks = break_vec,
     xlim = c(-max_val, max_val),
     main = "Adams: Histogram of test statistic",
     xlab = "Z-score", ylab = "Frequency", freq = T)
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(teststat_vec[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Other interest gene", "Housekeeping gene"),
       fill = c(2,4,3), cex = 0.6)
graphics.off()

png("../../../../out/fig/Writeup11/adams_ciliated_esvd_teststat_histogram_separate.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05
main_vec <- c("Published DE", "Other interest\n(DE other celltypes)", "Cell-cycle")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency", freq = T)
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()

#################

names(teststat_vec) <- colnames(mat)
save(teststat_vec, session_info, date_of_run,
     file = "../../../../out/Writeup11/Writeup11_adams_ciliated_esvd_teststat.RData")
