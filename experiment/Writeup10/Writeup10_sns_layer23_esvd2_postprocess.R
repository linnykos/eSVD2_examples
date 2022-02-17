rm(list=ls())
library(Seurat)
load("../../../../out/Writeup10/Writeup10_sns_layer23_esvd2.RData")

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
library_idx <- which(colnames(esvd_res_full$covariates) != "diagnosis_ASD")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,"diagnosis_ASD",drop = F],
                       esvd_res_full$b_mat[,"diagnosis_ASD",drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,library_idx],
  esvd_res_full$b_mat[,library_idx]
))

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
AplusAlpha <- mat + Alpha
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2
tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
quantile(tmp)
quantile(mean_mat_nolib)
quantile(Alpha)

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
png("../../../../out/fig/Writeup10/sns_layer23_esvd2_teststat_histogram.png", height = 1200, width = 1200,
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

png("../../../../out/fig/Writeup10/sns_layer23_esvd2_teststat_histogram_separate.png",
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
dens_val <- dens_val * 120/max(dens_val)
max_val <- max(abs(teststat_vec))
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05
png("../../../../out/fig/Writeup10/sns_layer23_esvd2_teststat_histogram_logconcave.png", height = 1200, width = 1200,
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
# gene_metadata <- read.csv("../../../../data/sns_autism/genes.tsv",
#                           sep = "\t", header = F)
# ensg_selected <- gene_metadata[gene_metadata[,2] %in% multtest_res$idx,1]
# ensg_backround <- gene_metadata[gene_metadata[,2] %in% colnames(mat),1]
# ensg_de <- gene_metadata[gene_metadata[,2] %in% de_gene_specific,1]
# write_genes(ensg_selected,
#             file = "../../../../out/Writeup10/Writeup10_sns_layer23_esvd_calibrate_ensg_selected.txt")
# write_genes(ensg_backround,
#             file = "../../../../out/Writeup10/Writeup10_sns_layer23_esvd_ensg_background.txt")
# write_genes(ensg_de,
#             file = "../../../../out/Writeup10/Writeup10_sns_layer23_esvd_ensg_velmeshev.txt")
write_genes(multtest_res$idx,
            file = "../../../../out/Writeup10/Writeup10_sns_layer23_esvd_calibrate_name_selected.txt")

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
y_max <- max(multtest_res$neglog_p_val)
x_max <- ceiling(max(abs(x_vec)))
png("../../../../out/fig/Writeup10/sns_layer23_esvd2_volcano_calibrate.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-x_max, x_max), ylim = range(0, y_max), bty = "n",
     main = "Volcano plot for Layer 2/3",
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

##########################

teststat_vec[multtest_res$idx]
multtest_res$idx
length(multtest_res$idx)
length(intersect(multtest_res$idx, de_gene_specific))
length(intersect(multtest_res$idx, de_genes))
length(intersect(multtest_res$idx, sfari_genes))
length(intersect(multtest_res$idx, hk_genes))
