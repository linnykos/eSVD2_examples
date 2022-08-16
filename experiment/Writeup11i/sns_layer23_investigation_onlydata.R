rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

######3

col_idx <- which(colnames(eSVD_obj$covariates) == "diagnosis_ASD")
lib_mat <- exp(tcrossprod(eSVD_obj$covariates[,-col_idx], eSVD_obj$fit_Second$z_mat[,-col_idx]))
lib_mat <- pmax(lib_mat, 1)

dat <- as.matrix(eSVD_obj$dat)+.5

rel_expr_mat <- dat/lib_mat

nat_mat1 <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,"diagnosis_ASD"], eSVD_obj$fit_Second$z_mat[,"diagnosis_ASD"])
fit_mat <- exp(nat_mat1 + nat_mat2)
numerator_mat <- sweep(fit_mat, MARGIN = 2, STATS = eSVD_obj$fit_Second$nuisance_vec, FUN = '*')
quantile(fit_mat)
quantile(numerator_mat)
quantile(eSVD_obj$fit_Second$nuisance_vec)

###############

metadata <- sns@meta.data
metadata[,"individual"] <- as.factor(metadata[,"individual"])
covariate_individual = "individual"

case_control_variable <- eSVD2:::.get_object(eSVD_obj = eSVD_obj, what_obj = "init_case_control_variable", which_fit = "param")
covariates <- eSVD2:::.get_object(eSVD_obj = eSVD_obj, what_obj = "covariates", which_fit = NULL)
cc_vec <- covariates[,case_control_variable]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])

individual_vec <- metadata[,covariate_individual]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])

tmp <- eSVD2:::.determine_individual_indices(case_individuals = case_individuals,
                                     control_individuals = control_individuals,
                                     covariate_individual = covariate_individual,
                                     metadata = metadata)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- eSVD2:::.construct_averaging_matrix(idx_list = all_indiv_idx,
                                       n = nrow(rel_expr_mat))
avg_posterior_mean_mat <- as.matrix(avg_mat %*% rel_expr_mat)
avg_posterior_var_mat <- t(sapply(1:nrow(avg_mat), function(i){
  idx <- eSVD2:::.nonzero_col(Matrix::t(avg_mat), col_idx = i, bool_value = F)
  matrixStats::colVars(rel_expr_mat[idx,])
}))

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
teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
  (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))
names(teststat_vec) <- colnames(rel_expr_mat)
eSVD_obj$teststat_vec <- teststat_vec

#############################

load("../../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L2/3"),]
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

hk_idx <- which(names(eSVD_obj$teststat_vec) %in% c(hk_genes, cycling_genes))
de_idx <- which(names(eSVD_obj$teststat_vec) %in% de_gene_specific)
other_idx <- which(names(eSVD_obj$teststat_vec) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj$teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

teststat_vec <- pmax(pmin(eSVD_obj$teststat_vec, 30), -30)
max_val <- max(abs(eSVD_obj$teststat_vec))
png("../../../../out/fig/Writeup11i/Writeup11i_sns_layer23_esvd_onlydata_teststat_histogram.png",
    height = 1200, width = 1200,
    units = "px", res = 300)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
hist(teststat_vec, breaks = break_vec,
     xlim = c(-max_val, max_val),
     main = "Histogram of test statistic",
     xlab = "Z-score", ylab = "Frequency", freq = T)
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(teststat_vec[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Other interest gene", "Housekeeping gene"),
       fill = c(2,4,3), cex = 0.6)
graphics.off()

png("../../../../out/fig/Writeup11i/Writeup11i_sns_layer23_esvd_onlydata_teststat_histogram_separate.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
main_vec <- c("Published DE", "Other interest\n(SFARI, DE other region)", "Housekeeping+Cell-cycle")
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






