rm(list=ls())
load("../../out/writeup8d/writeup8d_sns_layer23_esvd_extended2.RData")
hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]

nat_mat1 <- tcrossprod(esvd_res_nb2$x_mat, esvd_res_nb2$y_mat)
nat_mat2 <- tcrossprod(esvd_res_nb2$covariates, esvd_res_nb2$b_mat)
nat_mat2b <- tcrossprod(esvd_res_nb2$covariates[,c("Intercept", "diagnosis_ASD")], esvd_res_nb2$b_mat[,c("Intercept", "diagnosis_ASD")])
nat_mat <- nat_mat1 + nat_mat2
nat_matb <- nat_mat1 + nat_mat2b
mean_mat <- exp(nat_mat)
mean_matb <- exp(nat_matb)

p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)
quantile(angle_vec[gene_idx])

# j_vec <- which(angle_vec >= 65); j <- j_vec[which.max(colSums(mat[,j_vec]))]
# j <- which.min(abs(angle_vec - 70))
j <- 28
# j <- which.min(angle_vec)
col_pal <- colorRampPalette(c(2, 3))(10)
col_pal <- sapply(col_pal, function(x){
  paste0(x, "33")
})
library_grid <- seq(min(esvd_res_nb2$covariates[,"Log_UMI"]),
                    max(esvd_res_nb2$covariates[,"Log_UMI"]),
                    length.out = length(col_pal))
col_vec <- sapply(1:nrow(mat), function(i){
  col_pal[which.min(abs(esvd_res_nb2$covariates[i,"Log_UMI"] - library_grid))]
})
col_vec <- esvd_res_nb2$covariates[,"diagnosis_ASD"]+2
par(mfrow = c(1,2))
x_vec <- mean_matb[,j]
y_vec <- jitter(mat[,j])
xlim <- range(c(x_vec, y_vec))
plot(x_vec, y_vec, pch = 16, asp = T,
     col = col_vec, xlim = xlim, ylim = xlim)
hist(mat[,j], col = "gray", probability = T, breaks = seq(-.5, max(mat[,j]+.5), by = 1))
esvd_res_nb2$nuisance_param_vec[j]

par(mfrow = c(1,1))
cell_idx <- which(esvd_res_nb2$covariates[,"diagnosis_ASD"] == 1)
vioplot::vioplot(x_vec[cell_idx], x_vec[-cell_idx])

################

x_l2_vec_init <- apply(init_x_mat, 1, eSVD2:::.l2norm)
x_l2_vec <- apply(esvd_res_nb2$x_mat, 1, eSVD2:::.l2norm)
plot(x_l2_vec_init, x_l2_vec, asp = T, pch = 16)

y_l2_vec_init <- apply(init_y_mat, 1, eSVD2:::.l2norm)
y_l2_vec <- apply(esvd_res_nb2$y_mat, 1, eSVD2:::.l2norm)
plot(y_l2_vec_init, y_l2_vec, asp = T, pch = 16) # hm... it seems like a lot of vectors shrunk, but it's not clear if it's cuz due to the janky initialization
plot(y_l2_vec_init[gene_idx], y_l2_vec[gene_idx], asp = T, pch = 16)

####################3

# j <- which(colnames(mat) == de_genes[which.max(esvd_res_nb2$b_mat[de_genes,"diagnosis_ASD"])])
# j <- 28
# j <- grep("BTF3", colnames(mat))
# esvd_res_nb2$b_mat[j,"diagnosis_ASD"]

# first find the individuals
case_individuals <- unique(sns@meta.data[which(sns@meta.data$diagnosis == "ASD"),"individual"])
control_individuals <- unique(sns@meta.data[which(sns@meta.data$diagnosis == "Control"),"individual"])
case_idx <- which(esvd_res_nb2$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res_nb2$covariates[,"diagnosis_ASD"] == 0)

summary_stats <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')
  r_val <- esvd_res_nb2$nuisance_param_vec[j]

  # next find the cells, then compute one gaussian per individual
  case_gaussians <- sapply(case_individuals, function(indiv){
    cell_names <- rownames(sns@meta.data)[which(sns$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(mean_matb[cell_idx,j])
    var_val <- sum(mean_matb[cell_idx,j] + mean_matb[cell_idx,j]^2/r_val)/length(cell_idx)^2
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(sns@meta.data)[which(sns$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(mean_matb[cell_idx,j])
    var_val <- sum(mean_matb[cell_idx,j] + mean_matb[cell_idx,j]^2/r_val)/length(cell_idx)^2
    c(mean_val = mean_val, var_val = var_val)
  })

  case_gaussian <- list(mean_val = mean(case_gaussians[1,]),
                     var_val = sum(case_gaussians[2,])/ncol(case_gaussians)^2,
                     n = ncol(case_gaussians))
  control_gaussian <- list(mean_val = mean(control_gaussians[1,]),
                        var_val = sum(control_gaussians[2,])/ncol(control_gaussians)^2,
                        n = ncol(control_gaussians))

  n1 <- control_gaussian$n; n2 <- case_gaussian$n
  mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
  cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$mean_val
  combined_cov <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)
  test_stat <- (n1*n2)*(mean1 - mean2)/(combined_cov * (n1+n2))

  # using https://github.com/cran/Hotelling/blob/master/R/hotelling.test.r
  # and https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
  p <- 1; m <- n1+n2-2
  p_val <- 1-stats::pf((m-p+1)/(p*m)*test_stat, df1 = p, df2 = m)

  log_diff <- log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))

  c(p_val = p_val, log_diff = log_diff)
})

hk_idx <- which(colnames(mat) %in% hk_genes)
de_idx <- which(colnames(mat) %in% de_genes)
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), ncol(summary_stats))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
plot(NA, xlim = range(summary_stats[2,]), ylim = range(-log10(summary_stats[1,])))
points(x = summary_stats[2, -unique(c(hk_idx,de_idx))],
       y = -log10(summary_stats[1,-unique(c(hk_idx,de_idx))]),
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = summary_stats[2,hk_idx],
       y = -log10(summary_stats[1,hk_idx]),
       pch = 16, col = col_vec[hk_idx])
points(x = summary_stats[2,de_idx],
       y = -log10(summary_stats[1,de_idx]),
       pch = 16, col = col_vec[de_idx])

