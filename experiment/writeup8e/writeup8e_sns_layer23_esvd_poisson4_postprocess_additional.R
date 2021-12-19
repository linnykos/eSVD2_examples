rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson4.RData")
mat[mat == min(mat)] <- 0

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2_full <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2_full
mean_mat <- exp(nat_mat)

library_idx <- which(colnames(esvd_res_full$covariates) == "Log_UMI")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,-library_idx], esvd_res_full$b_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
mean_mat_nolib <- pmin(mean_mat_nolib, 1e4)

nuisance_param_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  val1 <- MASS::theta.ml(y = mat[,j], mu = mean_mat[,j])
  val2 <- MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
  val3 <- glmGamPoi::overdispersion_mle(y = mat[,j], mean = mean_mat[,j])$estimate

  vec <- c(val1, val2, val3)
  vec <- vec[!is.na(vec)]
  if(length(vec) == 1) return(vec[1])
  vec <- pmax(pmin(vec, 1e5), 0.1)
  vec <- c(vec, c(0.1, 0.5, 1))

  obs_prob <- length(which(mat[,j] == 0))/nrow(mat)
  target_prob_vec <- sapply(vec, function(val){
    mean((1+mean_mat[,j]/val)^(-val))
  })
  return(vec[which.min(abs(target_prob_vec - obs_prob))])
})

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})
library_mat <- pmin(library_mat, 50000)
AplusR <- sweep(mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(mean_mat_nolib, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
RoverMuplusS <- RoverMu + library_mat
posterior_mean_mat <- AplusR/RoverMuplusS

#########################

avg_exp <- matrixStats::colMeans2(mat)
set.seed(10)
gene_grouping <- stats::kmeans(avg_exp, centers = 6)
gene_grouping$size
library_vec_org <- matrixStats::rowSums2(mat)
n <- nrow(mat)

np_list <- lapply(1:6, function(i){
  print(i)
  idx <- which(gene_grouping$cluster == i)

  tmp_list <- lapply(1:length(idx), function(j){
    cell_idx <- sample(1:n, size = 1000)
    cbind(library_vec_org[cell_idx],  mat[cell_idx,idx[j]])
  })
  df <- as.data.frame(do.call(rbind, tmp_list))
  colnames(df) <- c("x", "y")
  np_fit <- npregfast::frfast(y ~ x,
                              data = df,
                              nboot = 100)

  tmp_res <- data.frame(x = df$x, y = df$y,
                        x_vec = np_fit$x, pred = np_fit$p[,1,1])
})

png("../../out/fig/writeup8e/sns_layer23_genes_summary.png",
    height = 2000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(2,3), mar = c(4, 4, 4, 0.5))
for(i in 1:6){
  print(i)
  idx <- which(gene_grouping$cluster == i)

  xlim <- range(np_list[[i]]$x)
  ylim <- c(-0.25, quantile(np_list[[i]]$y, prob = 0.99))
  x <- np_list[[i]]$x
  y <- np_list[[i]]$y
  rm_idx <- which(y >= ylim[2])
  x <- x[-rm_idx]; y <- y[-rm_idx]
  keep_idx <- sample(1:length(x), 1e5)
  x <- x[keep_idx]; y <- y[keep_idx]
  y <- y + runif(length(y), min = -0.25, max = 0.25)

  plot(x, y,
       pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
       xlim = xlim, ylim = ylim, xlab = "Total cell UMI count",
       ylab = "Gene UMI count", main = paste0("Gene group ", i, "\n(", length(idx), " genes)"))
  points(np_list[[i]]$x_vec, np_list[[i]]$pred, col = 2, cex = 2, pch = 16)
}
graphics.off()

