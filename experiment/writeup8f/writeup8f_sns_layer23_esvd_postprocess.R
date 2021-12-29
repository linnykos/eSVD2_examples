rm(list=ls())
load("../../../../out/writeup8f/writeup8f_sns_layer23_esvd.RData")
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]

#########

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

library_idx <- which(colnames(esvd_res_full$covariates) == "Log_UMI")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,-library_idx], esvd_res_full$b_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
mean_mat_nolib <- pmin(mean_mat_nolib, 1e4)

zero_prob <- apply(mat, 2, function(x){
  length(which(x == 0))/nrow(mat)
})
nuisance_param_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  val1 <- MASS::theta.ml(y = mat[,j], mu = mean_mat[,j])
  val2 <- MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
  val3 <- glmGamPoi::overdispersion_mle(y = mat[,j], mean = mean_mat[,j])$estimate

  vec <- c(val1, val2, val3)
  vec <- vec[!is.na(vec)]
  if(length(vec) == 1) return(vec[1])
  vec <- pmax(pmin(vec, 1e5), 1)
  vec <- c(vec, 1)

  obs_prob <- length(which(mat[,j] == 0))/nrow(mat)
  target_prob_vec <- sapply(vec, function(val){
    mean((1+mean_mat[,j]/val)^(-val))
  })
  return(vec[which.min(abs(target_prob_vec - obs_prob))])
})
quantile(nuisance_param_vec)
quantile(nuisance_param_vec[nuisance_param_vec>1])

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})
library_mat <- pmin(library_mat, 10000)
AplusR <- sweep(mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(mean_mat_nolib, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
RoverMuplusS <- RoverMu + library_mat
posterior_mean_mat <- AplusR/RoverMuplusS
posterior_var_mat <- AplusR/RoverMuplusS^2
tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
quantile(tmp)
tmp <- posterior_mean_mat/mean_mat_nolib
round(quantile(tmp),4)

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,c("Intercept", "diagnosis_ASD")], esvd_res_full$b_mat[,c("Intercept", "diagnosis_ASD")])
nat_mat_clean <- nat_mat1 + nat_mat2
mean_mat_clean <- exp(nat_mat_clean)

ratio_mat <- mean_mat_clean/mean_mat_nolib
quantile(ratio_mat)
posterior_mean_mat2 <- posterior_mean_mat * ratio_mat

case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
case_idx <- which(esvd_res_full$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res_full$covariates[,"diagnosis_ASD"] == 0)

individual_stats <- lapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')
  r_val <- nuisance_param_vec[j]

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
hk_idx <- which(colnames(mat) %in% hk_genes)
de_idx <- which(colnames(mat) %in% de_genes)

x_vec <- sapply(1:ncol(mat), function(j){
  # log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))
  log2(mean(mat[case_idx,j])) - log2(mean(mat[control_idx,j]))
})

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(p_val_vec))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

### let's draw it nicer
png("../../../../out/fig/writeup8f/sns_layer23_volcano.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-2,2), ylim = range(0, 10), bty = "n",
     main = "Volcano plot for Layer 2/3",
     xlab = "Log2 fold change (i.e., Log2 mean difference)", ylab = "-Log10(P value)")
for(x in seq(-2,2,by=0.5)){
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
legend("topright", c("Published DE gene", "Housekeeping gene", "Other"),
       fill = c(2,3,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()

max_val <- 10
zz <- 10^p_val_vec/2
zz[x_vec > 0] <- .5 + (.5-zz[x_vec > 0])
zz <- pmax(pmin(stats::qnorm(zz), max_val), -max_val)
png("../../../../out/fig/writeup8f/sns_layer23_zscore_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
hist(zz, breaks = seq(-max_val-0.05, max_val+0.05, by = 0.1), xlim = c(-5,5),
     main = "Histogram of two-sided Z-scores",
     xlab = "Z-score", ylab = "Frequency")
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(zz[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Housekeeping gene"),
       fill = c(2,3), cex = 0.6)
graphics.off()
