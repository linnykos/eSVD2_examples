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

factor_vec <- factor(c(rep("case", length(case_individuals)),
                       rep("control", length(control_individuals))))
teststat_vec <- sapply(1:length(individual_stats), function(j){
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

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

quantile(teststat_vec)
quantile(teststat_vec, probs = seq(0.9, 1, length.out=11))
colnames(mat)[order(abs(teststat_vec), decreasing = T)[1:20]]

teststat_vec <- pmax(pmin(teststat_vec, 30), -30)
max_val <- max(abs(teststat_vec))
png("../../../../out/fig/Writeup11/sns_invip_esvd_graphbased_teststat_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05
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

png("../../../../out/fig/Writeup11/sns_invip_esvd_graphbased_teststat_histogram_separate.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05
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

###########################

set.seed(10)
null_res <- logcondens::logConDens(teststat_vec[hk_idx],
                                   smoothed = T,
                                   print = FALSE,
                                   xs = seq(1.5*min(teststat_vec),
                                            1.5*max(teststat_vec),
                                            length.out = 1000))
dens_val <- null_res$f.smoothed
dens_val <- dens_val * 150/max(dens_val)
max_val <- max(abs(teststat_vec))
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05
png("../../../../out/fig/Writeup11/sns_invip_esvd_graphbased_teststat_histogram_logconcave.png", height = 1200, width = 1200,
    units = "px", res = 300)
hist(teststat_vec, breaks = break_vec,
     xlim = c(-max_val, max_val),
     main = "Histogram of two-sided Z-scores",
     xlab = "Z-score", ylab = "Frequency", freq = T)
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(teststat_vec[i], col = col_vec[i], lwd = 2)
}
lines(null_res$xs, dens_val, lwd = 2, lty = 2, col = 3)
legend("topright", c("Published DE gene", "Other interest gene", "Housekeeping gene"),
       fill = c(2,4,3), cex = 0.6)
graphics.off()

names(teststat_vec) <- colnames(mat)
multtest_res <- multttest_calibrate(teststat_vec = teststat_vec,
                                    null_dens = null_res$f.smoothed,
                                    null_x = null_res$xs,
                                    fdr_cutoff = 0.05,
                                    two_sided = T)

#########

indiv_list <- lapply(unique(metadata$individual), function(indiv){
  which(metadata$individual == indiv)
})
names(indiv_list) <- unique(metadata$individual)
case_individuals <- as.character(unique(metadata[which(metadata$diagnosis == "ASD"),"individual"]))
control_individuals <- as.character(unique(metadata[which(metadata$diagnosis == "Control"),"individual"]))
mat_avg <- t(sapply(indiv_list, function(idx_vec){
  matrixStats::colMeans2(mat[idx_vec,])
}))
rownames(mat_avg) <- names(indiv_list)
x_vec <- sapply(1:ncol(mat_avg), function(j){
  # log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))
  log2(mean(mat_avg[case_individuals,j])+1) - log2(mean(mat_avg[control_individuals,j])+1)
})

### let's draw it nicer
quantile(multtest_res$neglog_p_val)
y_max <- 20
x_max <- ceiling(max(abs(x_vec)))
png("../../../../out/fig/Writeup10/sns_invip_esvd_graphbased_volcano_calibrate.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-x_max, x_max), ylim = range(0, y_max), bty = "n",
     main = "Volcano plot for IN-VIP",
     xlab = "Log2 fold change (i.e., Log2 mean difference)", ylab = "-Log10(P value)")
for(x in seq(-x_max, x_max,by=0.5)){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
lines(rep(0,2), c(-1e5,1e5), col = "gray")
for(y in seq(0,max(multtest_res$neglog_p_val),by=2)){
  lines(c(-1e5,1e5), rep(y,2), lty = 2, col = "gray", lwd = 0.5)
}
points(x = x_vec[-unique(c(hk_idx,de_idx))],
       y = multtest_res$neglog_p_val[-unique(c(hk_idx,de_idx))],
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = x_vec[shuf_idx],
       y = multtest_res$neglog_p_val[shuf_idx],
       pch = 16, col = "white", cex = 1.5)
points(x = x_vec[shuf_idx],
       y = multtest_res$neglog_p_val[shuf_idx],
       pch = 16, col = col_vec[shuf_idx])
lines(x = c(-2*x_max, 2*x_max), y = rep(min(multtest_res$neglog_p_val[multtest_res$idx]),2),
      col = 2, lwd = 2, lty = 2)
legend("topleft", c("Published DE gene", "Other interest gene", "Housekeeping gene", "Other"),
       fill = c(2,4,3,rgb(0.5,0.5,0.5)), cex = 0.5)
graphics.off()
