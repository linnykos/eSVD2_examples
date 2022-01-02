rm(list=ls())
load("../../../../out/writeup8f/writeup8f_sns_pseudoreal_esvd_poisson.RData")
de_idx <- true_objects$autism_gene_idx

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

library_idx <- which(colnames(esvd_res_full$covariates) == "Log_UMI")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,-library_idx], esvd_res_full$b_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
mean_mat_nolib <- pmin(mean_mat_nolib, 1e4)

indiv_list <- lapply(unique(metadata$individual), function(indiv){
  which(metadata$individual == indiv)
})
mat_avg <- t(sapply(indiv_list, function(idx_vec){
  matrixStats::colMeans2(mat[idx_vec,])
}))
mean_avg_nolib <- t(sapply(indiv_list, function(idx_vec){
  matrixStats::colMeans2(mean_mat_nolib[idx_vec,])
}))

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})
library_avg <- t(sapply(indiv_list, function(idx_vec){
  matrixStats::colMeans2(library_mat[idx_vec,])
}))

res_list <- sapply(1:ncol(mat), function(j){
  print(j)
  calculate_fano_parameter(y = mat_avg[,j],
                           mu = mean_avg_nolib[,j],
                           sf = library_avg[,j],
                           max_val = 1e4)
})
quantile(res_list[1,])

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = res_list[1,], FUN = "*")
AplusAlpha <- mat + Alpha
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = res_list[1,], FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2
tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
quantile(tmp)

###########

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
tmp_idx <- c(which(colnames(covariates) %in% c("Intercept", "diagnosis_ASD")))
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,tmp_idx], esvd_res_full$b_mat[,tmp_idx])
nat_mat_clean <- nat_mat1 + nat_mat2
mean_mat_clean <- exp(nat_mat_clean)

ratio_mat <- mean_mat_clean/mean_mat_nolib
quantile(ratio_mat)
posterior_mean_mat2 <- posterior_mean_mat * ratio_mat

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

x_vec <- sapply(1:ncol(mat), function(j){
  # log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))
  log2(mean(mat[case_idx,j])) - log2(mean(mat[control_idx,j]))
})

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(p_val_vec))
col_vec[de_idx] <- 2
shuf_idx <- c(de_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

### let's draw it nicer
x_max <- max(abs(x_vec))
y_max <- 40
png("../../../../out/fig/writeup8f/sns_pseudoreal_volcano_fano.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-x_max, x_max), ylim = range(0, y_max), bty = "n",
     main = "Volcano plot for Layer 2/3  (Alt)",
     xlab = "Log2 fold change (i.e., Log2 mean difference)", ylab = "-Log10(P value)")
for(x in c(seq(0, x_max, by = 0.5), seq(0, -x_max, by = -0.5))){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
lines(rep(0,2), c(-1e5,1e5), col = "gray")
for(y in seq(0,max(-p_val_vec),by=1)){
  lines(c(-1e5,1e5), rep(y,2), lty = 2, col = "gray", lwd = 0.5)
}
points(x = x_vec[-unique(c(de_idx))],
       y = -p_val_vec[-unique(c(de_idx))],
       pch = 16, col = col_vec[-unique(c(de_idx))])
points(x = x_vec[shuf_idx],
       y = -p_val_vec[shuf_idx],
       pch = 16, col = "white", cex = 1.5)
points(x = x_vec[shuf_idx],
       y = -p_val_vec[shuf_idx],
       pch = 16, col = col_vec[shuf_idx])
legend("topleft", c("Published DE gene", "Other"),
       fill = c(2,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()

max_val <- 10
zz <- 10^p_val_vec/2
zz[x_vec > 0] <- .5 + (.5-zz[x_vec > 0])
zz <- pmax(pmin(stats::qnorm(zz), max_val), -max_val)
png("../../../../out/fig/writeup8f/sns_pseudoreal_zscore_histogram_fano.png", height = 1200, width = 1200,
    units = "px", res = 300)
hist(zz, breaks = seq(-max_val-0.05, max_val+0.05, by = 0.1), xlim = c(-max_val,max_val),
     main = "Histogram of two-sided Z-scores (Alt)",
     xlab = "Z-score", ylab = "Frequency")
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(zz[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene"),
       fill = c(2), cex = 0.6)
graphics.off()

#############################


gene_list <- list()
gene_list[[1]] <- true_objects$up_idx[order(p_val_vec[true_objects$up_idx], decreasing = F)[1:2]]
gene_list[[2]] <- true_objects$up_idx[order(p_val_vec[true_objects$up_idx], decreasing = T)[1:2]]
gene_list[[3]] <- true_objects$down_idx[order(p_val_vec[true_objects$down_idx], decreasing = F)[1:2]]
gene_list[[4]] <- true_objects$down_idx[order(p_val_vec[true_objects$down_idx], decreasing = T)[1:2]]

non_de_idx <- c(1:ncol(mat))[-true_objects$autism_gene_idx]
low_lib_idx <- intersect(which(abs(x_vec) <= 0.05), non_de_idx)
high_lib_idx <- intersect(which(abs(x_vec) >= 0.05), non_de_idx)
gene_list[[5]] <- low_lib_idx[order(p_val_vec[low_lib_idx], decreasing = F)[1:2]]
gene_list[[6]] <- high_lib_idx[order(p_val_vec[high_lib_idx], decreasing = F)[1:2]]
gene_list[[7]] <- non_de_idx[order(p_val_vec[non_de_idx], decreasing = T)[1:2]]
names(gene_list) <- c("Upexpressed_Sig",
                      "Upexpressed_NotSig",
                      "Downexpressed_Sig",
                      "Downexpressed_NotSig",
                      "NotDE_SmallLib_Sig",
                      "NotDE_LargeLib_Sig",
                      "NotDE_NotSig")

zero_prop <- apply(mat, 2, function(x){length(which(x == 0))/length(x)})
case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
all_individuals <- c(case_individuals, control_individuals)
col_vec_individuals <- c(rep(2, length(case_individuals)), rep(3, length(control_individuals)))
shuf_idx <- sample(1:length(col_vec_individuals))
all_individuals <- all_individuals[shuf_idx]; col_vec_individuals <- col_vec_individuals[shuf_idx]

for(k in 1:length(gene_list)){
  for(j in 1:length(gene_list[[k]])){
    png(paste0("../../../../out/fig/writeup8f/sns_pseudoreal_gene_fano_", names(gene_list)[k], "_", j, ".png"),
        height = 2500, width = 2500,
        units = "px", res = 300)
    par(mfrow = c(2,2), mar = c(4,4,4,0.5))
    idx <- gene_list[[k]][j]

    tmp <- cbind(mean_mat[,idx], mat[,idx])
    xlim <- range(tmp)
    plot(tmp[,1], tmp[,2], asp = T,
         xlim = xlim, ylim = xlim,
         xlab = "Predicted mean (full)", ylab = "Observed value",
         main = paste0(names(gene_list)[k], "(", j, ", Idx: ", idx, ")", "\n0-percentage: ",
                       round(zero_prop[idx], 2)),
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
         main = paste0("Posterior mean fit\nSize coefficient: ", round(esvd_res_full$b_mat[idx,"Log_UMI"], 2),
                       ", Original LFC: ", round(x_vec[idx], 2)),
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
}
