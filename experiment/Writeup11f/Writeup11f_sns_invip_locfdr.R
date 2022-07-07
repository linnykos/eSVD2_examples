rm(list=ls())
load("../../../out/main/sns_invip_processed.RData")

library(Seurat)
library(eSVD2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# compute degrees of freedom
input_obj <- eSVD_obj
metadata <- sns@meta.data
covariate_individual <- "individual"

case_control_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
covariates <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
cc_vec <- covariates[,case_control_variable]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])

latest_Fit <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
posterior_mean_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_mean_mat", which_fit = latest_Fit)
posterior_var_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_var_mat", which_fit = latest_Fit)

individual_vec <- metadata[,covariate_individual]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])
tmp <- eSVD2:::.determine_individual_indices(case_individuals = case_individuals,
                                     control_individuals = control_individuals,
                                     covariate_individual = covariate_individual,
                                     metadata = metadata)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- eSVD2:::.construct_averaging_matrix(idx_list = all_indiv_idx,
                                       n = nrow(posterior_mean_mat))
avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

case_row_idx <- 1:length(case_individuals)
control_row_idx <- (length(case_individuals)+1):nrow(avg_posterior_mean_mat)
case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_row_idx,,drop = F])
control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_row_idx,,drop = F])
case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[case_row_idx,,drop = F],
  avg_posterior_var_mat = avg_posterior_var_mat[case_row_idx,,drop = F]
)
control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[control_row_idx,,drop = F],
  avg_posterior_var_mat = avg_posterior_var_mat[control_row_idx,,drop = F]
)

n1 <- length(case_individuals)
n2 <- length(control_individuals)

# see https://www.theopeneducator.com/doe/hypothesis-Testing-Inferential-Statistics-Analysis-of-Variance-ANOVA/Two-Sample-T-Test-Unequal-Variance
df_vec <- (case_gaussian_var/n1 + control_gaussian_var/n2)^2/((case_gaussian_var/n1)^2/(n1-1) + (control_gaussian_var/n2)^2/(n2-1))

teststat_vec <- eSVD_obj$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

####################

locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
quantile(locfdr_res$fdr)

load("../../../data/sns_autism/velmeshev_genes.RData")
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
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(names(eSVD_obj$teststat_vec) %in% c(hk_genes, cycling_genes))
de_idx <- which(names(eSVD_obj$teststat_vec) %in% de_gene_specific)
other_idx <- which(names(eSVD_obj$teststat_vec) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj$teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

png("../../../out/fig/Writeup11f/Writeup11f_sns_invip_esvd5_locfdr.png",
    height = 1200, width = 1200,
    units = "px", res = 300)
plot(jitter(gaussian_teststat), jitter(-log10(locfdr_res$fdr)), pch = 16, col = rgb(0.5, 0.5, 0.5, 0.3),
     xlab = "Gaussian test-statistic", ylab = "-Log10 localFDR")
points(jitter(gaussian_teststat[shuf_idx]),
       jitter(-log10(locfdr_res$fdr[shuf_idx])),
       pch = 16,
       col = col_vec[shuf_idx], cex = 2)
graphics.off()

fdr_vec <- -log10(locfdr_res$fdr)
selected_genes <- names(teststat_vec)[which(fdr_vec >= 2)]
length(selected_genes)
length(which(selected_genes %in% de_gene_specific))
length(which(selected_genes %in% de_genes))
length(which(selected_genes %in% sfari_genes))
length(which(selected_genes %in% de_gene_specific))
length(which(selected_genes %in% hk_genes))
length(which(selected_genes %in% cycling_genes))

sink("../../../out/main/sns_invip_deGenes.txt")
for(i in 1:length(selected_genes)){
  cat(selected_genes[i])
  cat("\n")
}
sink()

