rm(list=ls())
load("../../../../out/writeup8g/writeup8g_sns_layer23_esvd.RData")
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
source("../writeup8f/fano_nuisance.R")

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

library_idx <- which(colnames(esvd_res_full$covariates) == "Log_UMI")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,-library_idx], esvd_res_full$b_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
mean_mat_nolib <- pmin(mean_mat_nolib, 1e4)

# mean_avg_nolib <- t(sapply(indiv_list, function(idx_vec){
#   matrixStats::colMeans2(mean_mat_nolib[idx_vec,])
# }))

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})
# library_avg <- t(sapply(indiv_list, function(idx_vec){
#   matrixStats::colMeans2(library_mat[idx_vec,])
# }))

res_list <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')
  # print(j)
  calculate_fano_parameter(y = mat[,j],
                           mu = mean_mat_nolib[,j],
                           sf = library_mat[,j],
                           max_val = 1e4,
                           min_val = 1)
})
quantile(res_list[1,])

save("mat", "mean_mat_nolib", "library_mat", "res_list", "esvd_res_full", "de_genes",
     file = "../../../../out/writeup8g/tmp.RData")

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = res_list[1,], FUN = "*")
AplusAlpha <- mat + Alpha
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = res_list[1,], FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2
tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
quantile(tmp)
quantile(mean_mat_nolib)
quantile(Alpha)

###########

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
tmp_idx <- c(which(colnames(esvd_res_full$covariates) %in% c("Intercept", "diagnosis_ASD")))
# tmp_idx <- c(tmp_idx, grep("individual", colnames(esvd_res_full$covariates)))
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,tmp_idx], esvd_res_full$b_mat[,tmp_idx])
nat_mat_clean <- nat_mat1 + nat_mat2
mean_mat_clean <- exp(nat_mat_clean)

ratio_mat <- mean_mat_clean/mean_mat_nolib
quantile(ratio_mat)
posterior_mean_mat2 <- posterior_mean_mat * ratio_mat
# posterior_var_mat <- posterior_var_mat * ratio_mat

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

    mean_val <- mean(posterior_mean_mat2[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat2[cell_idx,j])
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

p_val_vec <- sapply(1:length(group_stats), function(j){
  case_gaussian <- group_stats[[j]]$case_gaussian
  control_gaussian <- group_stats[[j]]$control_gaussian

  n1 <- control_gaussian$n; n2 <- case_gaussian$n
  mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
  cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$var_val

  combined_cov <- cov1/n1 + cov2/n2
  test_stat <- abs(mean1 - mean2)/sqrt(combined_cov)
  df <- (cov1/n1 + cov2/n2)^2/((cov1/n1)^2/(n1-1) + (cov2/n2)^2/(n2-1))

  stats::pt(test_stat, df = df, lower.tail = F, log.p = T)/log(10) + log10(2)
})

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
  log2(mean(mat_avg[case_individuals,j])) - log2(mean(mat_avg[control_individuals,j]))
})

hk_idx <- which(colnames(mat) %in% hk_genes)
de_idx <- which(colnames(mat) %in% de_genes)
sfari_idx <- which(colnames(mat) %in% sfari_genes)

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(p_val_vec))
col_vec[sfari_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, sfari_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

### let's draw it nicer
y_max <- ceiling(max(-p_val_vec))
x_max <- ceiling(max(abs(x_vec)))
png("../../../../out/fig/writeup8g/sns_layer23_esvd_volcano.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-x_max, x_max), ylim = range(0, y_max), bty = "n",
     main = "Volcano plot for Layer 2/3",
     xlab = "Log2 fold change (i.e., Log2 mean difference)", ylab = "-Log10(P value)")
for(x in seq(-x_max, x_max,by=0.5)){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
lines(rep(0,2), c(-1e5,1e5), col = "gray")
for(y in seq(0,max(-p_val_vec),by=2)){
  lines(c(-1e5,1e5), rep(y,2), lty = 2, col = "gray", lwd = 0.5)
}
points(x = x_vec[-unique(c(hk_idx,de_idx))],
       y = -p_val_vec[-unique(c(hk_idx,de_idx))],
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = x_vec[shuf_idx],
       y = -p_val_vec[shuf_idx],
       pch = 16, col = "white", cex = 1.5)
points(x = x_vec[shuf_idx],
       y = -p_val_vec[shuf_idx],
       pch = 16, col = col_vec[shuf_idx])
legend("topright", c("Published DE gene", "SFARI gene", "Housekeeping gene", "Other"),
       fill = c(2,4,3,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()

max_val <- 5
zz <- 10^p_val_vec/2
zz[x_vec > 0] <- .5 + (.5-zz[x_vec > 0])
zz <- pmax(pmin(stats::qnorm(zz), max_val), -max_val)
png("../../../../out/fig/writeup8g/sns_layer23_esvd_zscore_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
hist(zz, breaks = seq(-max_val-0.05, max_val+0.05, by = 0.1), xlim = c(-5,5),
     main = "Histogram of two-sided Z-scores",
     xlab = "Z-score", ylab = "Frequency")
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(zz[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "SFARI gene", "Housekeeping gene"),
       fill = c(2,4,3), cex = 0.6)
graphics.off()

png("../../../../out/fig/writeup8g/sns_layer23_esvd_zscore_histogram_separate.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
main_vec <- c("Published DE", "SFARI", "Housekeeping")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(zz[idx], breaks = seq(-max_val-0.05, max_val+0.05, by = 0.1), xlim = c(-5,5),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency")
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(zz[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()


#############################
zero_prop <- apply(mat, 2, function(x){length(which(x == 0))/length(x)})

case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
all_individuals <- c(case_individuals, control_individuals)
col_vec_individuals <- c(rep(2, length(case_individuals)), rep(3, length(case_individuals)))
shuf_idx <- sample(1:length(col_vec_individuals))
all_individuals <- all_individuals[shuf_idx]; col_vec_individuals <- col_vec_individuals[shuf_idx]
gene_indices <- c(which(colnames(mat) == "SAT2"),
                  which(colnames(mat) == "DEXI"),
                  which.min(p_val_vec),
                  which.max(abs(x_vec)),
                  which.max(zero_prop),
                  which.min(zero_prop),
                  which.min(abs(zero_prop - 0.9)),
                  which.min(abs(zero_prop - 0.2)))
filename_vec <- c("truede", "nonde", "mostsigp", "largestx",
                  "maxzero", "leastzero", "0.9zero", "0.2zero")

for(kk in 1:length(gene_indices)){
  idx <- gene_indices[kk]
  png(paste0("../../../../out/fig/writeup8g/sns_layer23_esvd_gene_", filename_vec[kk], ".png"),
      height = 2500, width = 2500,
      units = "px", res = 300)
  par(mfrow = c(2,2), mar = c(4,4,4,0.5))

  tmp <- cbind(mean_mat[,idx], mat[,idx])
  xlim <- range(tmp)
  plot(tmp[,1], tmp[,2], asp = T,
       xlim = xlim, ylim = xlim,
       xlab = "Predicted mean (full)", ylab = "Observed value",
       main = paste0("eSVD fit for ", colnames(mat)[idx], "\n0-percentage: ",
                     round(zero_prop[idx], 2), ", Beta: ", round(res_list[1,idx], 2)),
       pch = 16, col = rgb(0.5,0.5,0.5,0.1))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = "white", pch = 16, cex = 3)
  }
  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = col_vec_individuals[i], pch = 16, cex = 2)
  }

  ###

  tmp <- cbind(posterior_mean_mat[,idx], mat[,idx]/library_mat[,idx])
  xlim <- range(tmp)
  plot(tmp[,1], tmp[,2], asp = T,
       xlim = xlim, ylim = xlim,
       xlab = "Posterior mean (full)", ylab = "Observed relative expression",
       main = paste0("Posterior mean fit\nSize coefficient: ", round(esvd_res_full$b_mat[idx,"Log_UMI"], 2)),
       pch = 16, col = rgb(0.5,0.5,0.5,0.1))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = "white", pch = 16, cex = 3)
  }
  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = col_vec_individuals[i], pch = 16, cex = 2)
  }

  ###

  tmp <- cbind(mean_mat_clean[,idx], mean_mat_nolib[,idx])
  xlim <- range(tmp)
  plot(tmp[,1], tmp[,2], asp = T,
       xlim = xlim, ylim = xlim,
       xlab = "Predicted mean (Clean)", ylab = "Predicted mean (No lib)",
       main = "Posterior mean with and without covariates",
       pch = 16, col = rgb(0.5,0.5,0.5,0.1))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = "white", pch = 16, cex = 3)
  }
  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = col_vec_individuals[i], pch = 16, cex = 2)
  }

  ###

  xlim <- range(posterior_mean_mat2[,idx], sqrt(posterior_var_mat[,idx]))
  plot(posterior_mean_mat2[,idx], sqrt(posterior_var_mat[,idx]),
       xlab = "Posterior mean (Clean)", ylab = "Posterior standard deviation",
       pch = 16, col = rgb(0.5,0.5,0.5,0.1),
       xlim = c(0, max(posterior_mean_mat2[,idx])),
       ylim = c(0, max(sqrt(posterior_var_mat[,idx]))),
       main = paste0("Posterior standard deviation vs. mean\n-Log10(P value): ",
                     round(-p_val_vec[idx], 2)))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  tmp <- cbind(individual_stats[[idx]]$case_gaussians, individual_stats[[idx]]$control_gaussians)
  col_vec <- c(rep(2, ncol(individual_stats[[idx]]$case_gaussians)),
               rep(3, ncol(individual_stats[[idx]]$control_gaussians)))
  shuf_idx <- sample(1:length(col_vec))
  tmp <- tmp[,shuf_idx]; col_vec <- col_vec[shuf_idx]
  points(tmp[1,], sqrt(tmp[2,]), pch = 16, col = "white", cex = 3)
  points(tmp[1,], sqrt(tmp[2,]), pch = 16, col = col_vec, cex = 2)

  graphics.off()
}
