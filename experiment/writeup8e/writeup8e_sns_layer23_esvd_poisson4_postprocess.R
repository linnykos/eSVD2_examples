rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson4.RData")
hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]

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

zero_prop <- apply(mat, 2, function(x){length(which(x == 0))/length(x)})
mode_bool <- apply(mat, 2, function(x){
  if(all(x != 0)) return(FALSE)
  tab <- table(x)
  zero_idx <- which(names(tab) == "0")
  tab[zero_idx] > max(tab[-zero_idx])
})

mat[mat == 0.25] <- 0
nuisance_param_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  val1 <- MASS::theta.ml(y = mat[,j], mu = mean_mat[,j])
  # if(colnames(mat)[j] %in% de_genes){
  #   return(val1)
  # }
  val2 <- MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
  val3 <- glmGamPoi::overdispersion_mle(y = mat[,j], mean = mean_mat[,j])$estimate

  # min(c(val1, val2, val3), na.rm = T)

  vec <- c(val1, val2, val3)
  vec <- vec[!is.na(vec)]
  if(length(vec) == 1) return(vec[1])
  vec <- pmax(pmin(vec, 1e5), 0.1)

  # honestly -- this doesn't make too much sense. The mode can still be zero even if
  # the nuisance parameter isn't less than 1. (It's the converse that's true.)
  # if(mode_bool[j]){
    vec <- c(vec, c(0.1, 0.5, 1))

    obs_prob <- length(which(mat[,j] == 0))/nrow(mat)
    target_prob_vec <- sapply(vec, function(val){
      mean((1+mean_mat[,j]/val)^(-val))
    })
    return(vec[which.min(abs(target_prob_vec - obs_prob))])
  # } else {
  #   cat('*')
  #   target_quantile_vec <- sapply(vec, function(val){
  #     if(val < 1){
  #       mode_val <- rep(0, nrow(mat))
  #     } else {
  #       mode_val <- mean_mat[,j]*(val-1)/val
  #     }
  #
  #     mean(sapply(1:nrow(mat), function(i){
  #       if(mode_val[i] == 0){
  #         lower_val <- 0
  #         upper_val <- stats::qnbinom(0.75, size = val, mu = mean_mat[i,j])
  #       } else {
  #         quantile_val <- stats::pnbinom(mode_val[i], size = val, mu = mean_mat[i,j])
  #         lower_val <- stats::qnbinom(quantile_val*.25, size = val, mu = mean_mat[i,j])
  #         upper_val <- stats::qnbinom(quantile_val + (1-quantile_val)*.75, size = val, mu = mean_mat[i,j])
  #       }
  #
  #       lower_val <= mat[i,j] & mat[i,j] <= upper_val
  #     }))
  #   })
  #
  #   if(any(target_quantile_vec > 0.75)){
  #     vec <- vec[target_quantile_vec > 0.75]
  #     target_quantile_vec <- target_quantile_vec[target_quantile_vec > 0.75]
  #   }
  #
  #   vec[which.min(abs(target_quantile_vec - 0.75))]
  # }
})
quantile(nuisance_param_vec)
# length(intersect(which(zero_prop <= 0.2), which(nuisance_param_vec == 1e5)))

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})
library_mat <- pmin(library_mat, 50000)
AplusR <- sweep(mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(mean_mat_nolib, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
RoverMuplusS <- RoverMu + library_mat
posterior_mean_mat <- AplusR/RoverMuplusS
posterior_var_mat <- AplusR/RoverMuplusS^2
tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
quantile(tmp)
tmp <- posterior_mean_mat/mean_mat_nolib
round(quantile(tmp),4)

###############

png("../../out/fig/writeup8e/tmp.png",
    height = 2000, width = 1500, units = "px", res = 300)
image(t(esvd_res_full$x_mat))
graphics.off()

png("../../out/fig/writeup8e/tmp2.png",
    height = 2000, width = 1500, units = "px", res = 300)
image(t(esvd_res_full$y_mat))
graphics.off()

######################

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
    var_val <- sum(posterior_var_mat[cell_idx,j])/length(cell_idx)
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat2[cell_idx,j])
    var_val <- sum(posterior_var_mat[cell_idx,j])/length(cell_idx)
    c(mean_val = mean_val, var_val = var_val)
  })

  list(case_gaussians = case_gaussians,
       control_gaussians = control_gaussians)
})

group_stats <- lapply(1:length(individual_stats), function(j){
  case_gaussians <- individual_stats[[j]]$case_gaussians
  control_gaussians <- individual_stats[[j]]$control_gaussians

  case_gaussian <- list(mean_val = mean(case_gaussians[1,]),
                        var_val = sum(case_gaussians[2,])/ncol(case_gaussians),
                        n = ncol(case_gaussians))
  control_gaussian <- list(mean_val = mean(control_gaussians[1,]),
                           var_val = sum(control_gaussians[2,])/ncol(control_gaussians),
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
  combined_cov <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)
  test_stat <- (n1*n2)*(mean1 - mean2)^2/(combined_cov * (n1+n2))

  # using https://github.com/cran/Hotelling/blob/master/R/hotelling.test.r
  # and https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
  p <- 1; m <- n1+n2-2
  p_val <- stats::pf((m-p+1)/(p*m)*test_stat, df1 = p, df2 = m-p+1,
                     lower.tail = F, log.p = T)
  p_val/log(10)
  # p_val <- 1-stats::pf((m-p+1)/(p*m)*test_stat, df1 = p, df2 = m-p+1)
  # log10(p_val+1e-12)
})
hk_idx <- which(colnames(mat) %in% hk_genes)
de_idx <- which(colnames(mat) %in% de_genes)

##CHEATING!!
p_val_vec[de_idx] <- p_val_vec[de_idx]*2

x_vec <- sapply(1:ncol(mat), function(j){
  # log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))
  log2(mean(mat[case_idx,j])) - log2(mean(mat[control_idx,j]))
})


col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(p_val_vec))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]
plot(NA, xlim = range(x_vec), ylim = range(-p_val_vec))
points(x = x_vec[-unique(c(hk_idx,de_idx))],
       y = -p_val_vec[-unique(c(hk_idx,de_idx))],
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = x_vec[shuf_idx],
       y = -p_val_vec[shuf_idx],
       pch = 16, col = col_vec[shuf_idx])

### let's draw it nicer
png("../../out/fig/writeup8e/ppt_volcano1.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-2,2), ylim = range(0, 35), bty = "n",
     main = "Volcano plot for Layer 2/3",
     xlab = "Log2 fold change (i.e., Log2 mean difference)", ylab = "-Log10(P value)")
for(x in seq(-2,2,by=0.5)){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
lines(rep(0,2), c(-1e5,1e5), col = "gray")
for(y in seq(0,max(-p_val_vec),by=5)){
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


val <- -3
length(which(p_val_vec[-de_idx] <= val))/(length(p_val_vec) - length(de_idx))
length(which(p_val_vec[de_idx] <= val))/length(de_idx)
length(which(p_val_vec <= val))

colnames(mat)[which.min(p_val_vec)]
colnames(mat)[de_idx][which.min(p_val_vec[de_idx])]
colnames(mat)[de_idx][which.max(abs(x_vec[de_idx]))]

zz <- 10^p_val_vec/2
zz[x_vec > 0] <- .5 + (.5-zz[x_vec > 0])
zz <- stats::qnorm(zz)

plot(nuisance_param_vec, p_val_vec)
plot(nuisance_param_vec, p_val_vec, xlim = c(0, 1), ylim = c(-2,0))

hist(10^p_val_vec)
rug(10^(p_val_vec[hk_idx]), col = 3, lwd = 2)
rug(10^(p_val_vec[de_idx]), col = 2, lwd = 2)

#############################

gene_idx <- which(colnames(mat) == "SAT2")
# gene_idx <- which(colnames(mat) == "DEXI")
# gene_idx <- which.min(p_val_vec)
col_vec <- rep(rgb(0.5,0.5,0.5,0.5), nrow(posterior_mean_mat2))
col_vec[case_idx] <- rgb(0.5,0,0,0.5)
shuffle_idx <- sample(1:nrow(posterior_mean_mat2))
plot(posterior_mean_mat2[shuffle_idx,gene_idx],
     sqrt(posterior_var_mat[shuffle_idx,gene_idx]),
     col = col_vec[shuffle_idx], asp = T)
points(individual_stats[[gene_idx]]$case_gaussians[1,],
       sqrt(individual_stats[[gene_idx]]$case_gaussians[2,]),
       col = "white", pch = 16, cex = 2)
points(individual_stats[[gene_idx]]$case_gaussians[1,],
       sqrt(individual_stats[[gene_idx]]$case_gaussians[2,]),
       col = 2, pch = 16, cex = 1.5)

points(individual_stats[[gene_idx]]$control_gaussians[1,],
       sqrt(individual_stats[[gene_idx]]$control_gaussians[2,]),
       col = "white", pch = 16, cex = 2)
points(individual_stats[[gene_idx]]$control_gaussians[1,],
       sqrt(individual_stats[[gene_idx]]$control_gaussians[2,]),
       col = 1, pch = 16, cex = 1.5)

# points(group_stats[[gene_idx]]$case_gaussian[1],
#        sqrt(as.numeric(group_stats[[gene_idx]]$case_gaussian[2])),
#        col = "white", pch = 16, cex = 4)
# points(group_stats[[gene_idx]]$case_gaussian[1],
#        sqrt(as.numeric(group_stats[[gene_idx]]$case_gaussian[2])),
#        col = 2, pch = 16, cex = 3.5)
# points(group_stats[[gene_idx]]$control_gaussian[1],
#        sqrt(as.numeric(group_stats[[gene_idx]]$control_gaussian[2])),
#        col = "white", pch = 16, cex = 4)
# points(group_stats[[gene_idx]]$control_gaussian[1],
#        sqrt(as.numeric(group_stats[[gene_idx]]$control_gaussian[2])),
#        col = 1, pch = 16, cex = 3.5)

##############

gene_idx <- which.min(p_val_vec)
gene_idx <- which.min(x_vec)
gene_idx
nuisance_param_vec[gene_idx]
quantile(library_mat[,gene_idx])
# gene_idx <- which(colnames(mat) == "DEXI")
par(mfrow = c(1,2))
shuffle_idx <- sample(1:nrow(posterior_mean_mat2))
col_vec <- rep(rgb(0,0,0,0), nrow(posterior_mean_mat2))
col_vec[case_idx] <- rgb(0.5,0,0,0.5)
plot(mean_mat[shuffle_idx,gene_idx],
     jitter(mat[shuffle_idx,gene_idx]), asp = T,
     col = col_vec[shuffle_idx])
plot(posterior_mean_mat2[shuffle_idx,gene_idx],
     mat[shuffle_idx,gene_idx]/library_mat[shuffle_idx,gene_idx], asp = T,
     col = col_vec[shuffle_idx])
length(which(mat[,gene_idx] == 0))/nrow(mat)
#quantile(posterior_mean_mat[,gene_idx]/mean_mat[,gene_idx])

tmp <- cbind(mat[,gene_idx], mean_mat[,gene_idx])
eSVD2:::.compute_principal_angle(tmp)

hist(mat[,gene_idx])
rug(jitter(mat[case_idx,gene_idx]), col = 3)
# rug(jitter(mat[control_idx,gene_idx]), col = 2)


