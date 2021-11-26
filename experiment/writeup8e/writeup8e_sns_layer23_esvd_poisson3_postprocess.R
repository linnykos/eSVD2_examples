rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson3.RData")

hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]

############

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates, esvd_res$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

nat_mat_wlibrary <- sweep(nat_mat, MARGIN = 1, STATS = esvd_res$offset_vec, FUN = "+")
mean_mat_wlibrary <- exp(nat_mat_wlibrary)
p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)

##########

hist(esvd_res$b_mat[,"diagnosis_ASD"])
rug(esvd_res$b_mat[which(colnames(mat) %in% hk_genes),"diagnosis_ASD"], col = 3)
rug(esvd_res$b_mat[which(colnames(mat) %in% de_genes),"diagnosis_ASD"], col = 2)

##############

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates, esvd_res$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

nuisance_param_vec <- sapply(1:p, function(j){
  MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
})
quantile(nuisance_param_vec)
nuisance_param_vec <- pmin(nuisance_param_vec, 1e5)

# library_mat <- tcrossprod(exp(esvd_res$covariates[,"Log_UMI",drop = F]), exp(esvd_res$b_mat[,"Log_UMI",drop = F]))*ncol(mat)
# AplusR <- sweep(mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
# RoverMu <- 1/sweep(mean_mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
# RoverMuplusS <- RoverMu + library_mat
# posterior_mean_mat <- AplusR/RoverMuplusS
# posterior_var_mat <- AplusR/RoverMuplusS^2
# tmp <- posterior_mean_mat/posterior_var_mat
# quantile(tmp)

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res$covariates[,"Log_UMI",drop = F]*esvd_res$b_mat[j,"Log_UMI"])*ncol(mat)
})
library_mat <- pmin(library_mat, 50000)
A <- mat/library_mat
AplusR <- sweep(A, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(mean_mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
RoverMuplusS <- RoverMu + 1
posterior_mean_mat <- AplusR/RoverMuplusS
posterior_var_mat <- AplusR/RoverMuplusS^2
tmp <- posterior_mean_mat/posterior_var_mat
quantile(tmp)

#######################

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates[,c("Intercept", "diagnosis_ASD")], esvd_res$b_mat[,c("Intercept", "diagnosis_ASD")])
nat_mat_clean <- nat_mat1 + nat_mat2
mean_mat_clean <- exp(nat_mat_clean)

ratio_mat <- mean_mat_clean/mean_mat
posterior_mean_mat2 <- posterior_mean_mat * ratio_mat

case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
case_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 0)

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

x_vec <- sapply(1:ncol(mat), function(j){
  log(mean(posterior_mean_mat2[case_idx,j])) - log(mean(posterior_mean_mat2[control_idx,j]))
})

hk_idx <- which(colnames(mat) %in% hk_genes)
de_idx <- which(colnames(mat) %in% de_genes)
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(p_val_vec))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
plot(NA, xlim = range(x_vec), ylim = range(-p_val_vec))
points(x = x_vec[-unique(c(hk_idx,de_idx))],
       y = -p_val_vec[-unique(c(hk_idx,de_idx))],
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = x_vec[hk_idx],
       y = -p_val_vec[hk_idx],
       pch = 16, col = col_vec[hk_idx])
points(x = x_vec[de_idx],
       y = -p_val_vec[de_idx],
       pch = 16, col = col_vec[de_idx])

#############################3

# gene_idx <- which(colnames(mat) == "SAT2")
# gene_idx <- which(colnames(mat) == "MRPL20")
gene_idx <- which(colnames(mat) == "DEXI")
col_vec <- rep(rgb(0.5,0.5,0.5,0.5), nrow(posterior_mean_mat2))
col_vec[case_idx] <- rgb(0.5,0,0,0.5)
shuffle_idx <- sample(1:nrow(posterior_mean_mat2))
plot(posterior_mean_mat2[shuffle_idx,gene_idx],
     posterior_var_mat[shuffle_idx,gene_idx],
     col = col_vec[shuffle_idx], asp = T)
points(individual_stats[[gene_idx]]$case_gaussians[1,],
       individual_stats[[gene_idx]]$case_gaussians[2,],
       col = "white", pch = 16, cex = 2)
points(individual_stats[[gene_idx]]$case_gaussians[1,],
       individual_stats[[gene_idx]]$case_gaussians[2,],
       col = 2, pch = 16, cex = 1.5)

points(individual_stats[[gene_idx]]$control_gaussians[1,],
       individual_stats[[gene_idx]]$control_gaussians[2,],
       col = "white", pch = 16, cex = 2)
points(individual_stats[[gene_idx]]$control_gaussians[1,],
       individual_stats[[gene_idx]]$control_gaussians[2,],
       col = 1, pch = 16, cex = 1.5)

points(group_stats[[gene_idx]]$case_gaussian[1],
       group_stats[[gene_idx]]$case_gaussian[2],
       col = "white", pch = 16, cex = 4)
points(group_stats[[gene_idx]]$case_gaussian[1],
       group_stats[[gene_idx]]$case_gaussian[2],
       col = 2, pch = 16, cex = 3.5)
points(group_stats[[gene_idx]]$control_gaussian[1],
       group_stats[[gene_idx]]$control_gaussian[2],
       col = "white", pch = 16, cex = 4)
points(group_stats[[gene_idx]]$control_gaussian[1],
       group_stats[[gene_idx]]$control_gaussian[2],
       col = 1, pch = 16, cex = 3.5)

###

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
excluded_idx <- which(colnames(esvd_res$covariates) %in% c("Log_UMI"))
nat_mat2 <- tcrossprod(esvd_res$covariates[,-excluded_idx], esvd_res$b_mat[,-excluded_idx])
nat_mat <- nat_mat1+nat_mat2
mean_mat2 <- exp(nat_mat)

# gene_idx <- which(colnames(mat) == "SAT2")
# gene_idx <- which(colnames(mat) == "MRPL20")
gene_idx <- which(colnames(mat) == "DEXI")
plot(library_mat[,gene_idx]*mean_mat2[,gene_idx]/ncol(mat),
     mat[,gene_idx], asp = T,
     col = rgb(0,0,0,0.1))

zz <- mat/library_mat
quantile(zz)

######################

case_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 0)

p_val_vec <- sapply(1:ncol(mat), function(j){
  stats::wilcox.test(x = posterior_mean_mat2[case_idx,j],
                     y = posterior_mean_mat2[control_idx,j])$p.value
})
