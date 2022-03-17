rm(list=ls())
library(Seurat)
load("../../../../out/Writeup10/Writeup10_sns_invip_processed.RData")
load("../../../../out/Writeup11/Writeup11_sns_invip_esvd_fit.RData")
source("twosampletest.R")

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
case_control_variable = "diagnosis_ASD"
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

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
AplusAlpha <- mat + Alpha
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2
# tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
# quantile(tmp)
# quantile(mean_mat_nolib)
# quantile(Alpha)

#################

metadata <- sns@meta.data
case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
case_idx <- which(metadata[,"diagnosis"] == "ASD")
control_idx <- which(metadata[,"diagnosis"] == "Control")

individual_stats <- lapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  # next find the cells, then compute one gaussian per individual
  case_gaussians <- sapply(case_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
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

teststat_vec1 <- sapply(1:length(group_stats), function(j){
  case_gaussian <- group_stats[[j]]$case_gaussian
  control_gaussian <- group_stats[[j]]$control_gaussian

  n1 <- control_gaussian$n; n2 <- case_gaussian$n
  mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
  cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$var_val

  combined_cov <- cov1/n1 + cov2/n2
  (mean2 - mean1)/sqrt(combined_cov)
})


factor_vec <- factor(c(rep("case", length(case_individuals)),
                       rep("control", length(control_individuals))))
teststat_vec2 <- sapply(1:length(individual_stats), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  mean_vec <- c(individual_stats[[j]]$case_gaussians["mean_val",],
                individual_stats[[j]]$control_gaussians["mean_val",])
  var_vec <- c(individual_stats[[j]]$case_gaussians["var_val",],
               individual_stats[[j]]$control_gaussians["var_val",])

  compute_twosample_pvalue(factor_vec = factor_vec,
                           mean_vec = mean_vec,
                           var_vec = var_vec,
                           k = 3, verbose = 0)$test_stat
})

##########

load("../../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "IN-VIP"),]
de_gene_specific <- tmp[,"Gene name"]
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list[-1], function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
de_genes <- de_genes[!de_genes %in% de_gene_specific]
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(colnames(mat) %in% c(hk_genes, cycling_genes))
de_idx <- which(colnames(mat) %in% de_gene_specific)
other_idx <- which(colnames(mat) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(teststat_vec1))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

##################################

png("../../../../out/fig/Writeup11/sns_invip_esvd_pvalue_comparison.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(NA, xlim = range(teststat_vec1), ylim = range(teststat_vec2),
     xlab = "Mixture Gaussian method",
     ylab = "Graph-based method")
lines(x = c(-1e4, 1e4), y = rep(0,2), col = "red", lty = 2)
lines(x = rep(0,2), y = c(-1e4, 1e4), col = "red", lty = 2)
points(x = teststat_vec1[-shuf_idx], y = teststat_vec2[-shuf_idx],
       col = col_vec[-shuf_idx], pch = 16)
points(x = teststat_vec1[shuf_idx], y = teststat_vec2[shuf_idx],
       col = col_vec[shuf_idx], pch = 16)
graphics.off()
