rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta1-inflamed_esvd.RData")
load("../../../out/main/regevEpi_ta1-inflamed_esvd_illustration-npreg.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

total_expression_vec <- log1p(apply(eSVD_obj$dat, 2, sum))
set.seed(10)
kmean_res <- stats::kmeans(total_expression_vec, centers = 6)

cluster_vec <- kmean_res$cluster
cluster_mean <- kmean_res$centers
reorder_idx <- order(cluster_mean, decreasing = F)
p <- length(cluster_vec)
cluster_vec2 <- rep(NA, p)
for(i in 1:length(reorder_idx)){
  cluster_vec2[which(cluster_vec == reorder_idx[i])] <- i
}
cluster_vec <- cluster_vec2

n <- nrow(before_npreg_list)
lib_vec <- eSVD_obj$covariates[,"Log_UMI"]
library_breaks <- seq(min(lib_vec), max(lib_vec), length.out = 1000)
curve_list <- lapply(1:6, function(k){
  print(k)
  idx <- which(cluster_vec == k)
  quantile_mat <- t(sapply(1:n, function(i){
    stats::quantile(before_npreg_list[i,idx], probs = c(0.25, .5, 0.75))
  }))

  curve_mat <- sapply(1:3, function(ell){
    tmp_df <- data.frame(y = quantile_mat[,ell], x = lib_vec)
    reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
    x_vec <- library_breaks
    y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
    y_vec$Estimation[,"Pred"]
  })

  curve_mat
})

png("../../../out/fig/main/regevEpi_ta1-inflamed_npreg-before.png",
    height = 2250, width = 1100,
    units = "px", res = 500)
par(mfrow = c(3,2), mar = c(3,3,3,1))
xlim <- range(library_breaks)
for(i in 1:6){
  ylim <- log1p(c(0, max(curve_list[[i]])))
  ylim[2] <- max(0.01, ylim[2])
  plot(NA, xlim = xlim, ylim = ylim,
       main = paste0("Group ", i, "\n(", length(which(cluster_vec == i)), " genes)"),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  lines(x = c(-100,100), y = rep(log1p(mean(curve_list[[i]][,2])), 2),
        lty = 2)
  polygon(x = c(library_breaks, rev(library_breaks)),
          y = log1p(c(curve_list[[i]][,1], rev(curve_list[[i]][,3]))),
          col = "gray")
  lines(x = library_breaks, y = log1p(curve_list[[i]][,2]),
        lwd = 4, col = "white")
  lines(x = library_breaks, y = log1p(curve_list[[i]][,2]),
        lwd = 3, col = 2)
  axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
}
graphics.off()

############################

n <- nrow(after_npreg_list)
lib_vec <- eSVD_obj$covariates[,"Log_UMI"]
library_breaks <- seq(min(lib_vec), max(lib_vec), length.out = 1000)
curve_list2 <- lapply(1:6, function(k){
  print(k)
  idx <- which(cluster_vec == k)
  quantile_mat <- t(sapply(1:n, function(i){
    stats::quantile(after_npreg_list[i,idx], probs = c(0.25, .5, 0.75))
  }))

  curve_mat <- sapply(1:3, function(ell){
    tmp_df <- data.frame(y = quantile_mat[,ell], x = lib_vec)
    reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
    x_vec <- library_breaks
    y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
    y_vec$Estimation[,"Pred"]
  })

  curve_mat
})

png("../../../out/fig/main/regevEpi_ta1-inflamed_npreg-after.png",
    height = 2250, width = 1100,
    units = "px", res = 500)
par(mfrow = c(3,2), mar = c(3,3,3,1))
xlim <- range(library_breaks)
ylim <- log1p(c(0, max(curve_list[[6]])))
for(i in 1:6){
  plot(NA, xlim = xlim, ylim = ylim,
       main = paste0("Group ", i, "\n(", length(which(cluster_vec == i)), " genes)"),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  lines(x = c(-100,100), y = rep(log1p(mean(curve_list2[[i]][,2])), 2),
        lty = 2)
  polygon(x = c(library_breaks, rev(library_breaks)),
          y = log1p(c(curve_list2[[i]][,1], rev(curve_list2[[i]][,3]))),
          col = "gray")
  lines(x = library_breaks, y = log1p(curve_list2[[i]][,2]),
        lwd = 4, col = "white")
  lines(x = library_breaks, y = log1p(curve_list2[[i]][,2]),
        lwd = 3, col = 2)
  axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
  axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
}
graphics.off()
