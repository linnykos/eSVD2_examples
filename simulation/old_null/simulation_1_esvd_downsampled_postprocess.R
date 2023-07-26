rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../out/simulation/simulation_1.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

true_residual_mat <- apply(nat_mat, 2, function(y){
  df <- data.frame(cbind(y, covariates))
  colnames(df)[1] <- "y"
  lm_res <- stats::lm(y ~ . - 1, data = df)
  stats::residuals(lm_res) + stats::coef(lm_res)["cc"]*covariates[,"cc"]
})
nat_mat_safe <- nat_mat

load("../../out/simulation/simulation_1_esvd.RData")
r <- ncol(eSVD_obj$fit_Second$z_mat)
downsample_values <- seq(0.9, 0.2, by = -0.1)
quantile_list <- vector("list", length(downsample_values))

other_idx <- which(!colnames(eSVD_obj[["covariates"]]) %in% "cc_1")
nat_mat_other <- tcrossprod(eSVD_obj[["covariates"]][,other_idx], eSVD_obj[["fit_Second"]]$z_mat[,other_idx])
nat_mat_alt <- nat_mat_safe - nat_mat_other
n <- nrow(nat_mat_safe); p <- ncol(nat_mat_safe)
quantile_list[[1]] <- sapply(1:n, function(i){
  stats::cor(true_residual_mat[i,], nat_mat_alt[i,])
})
quantile(quantile_list[[1]])

for(k in 1:length(downsample_values)){
  downsample_value <- downsample_values[k]
  print(paste0("Working on ", downsample_value))
  load(paste0("../../out/simulation/simulation_1_esvd_downsampled-", downsample_value, ".RData"))

  nat_mat_other <- tcrossprod(eSVD_obj[["covariates"]][,other_idx], eSVD_obj[["fit_Second"]]$z_mat[,other_idx])
  nat_mat_alt <- nat_mat_safe - nat_mat_other

  quantile_list[[k+1]] <- sapply(1:n, function(i){
    stats::cor(true_residual_mat[i,], nat_mat_alt[i,])
  })
}

sapply(quantile_list, quantile)
eSVD_quantile_mat <- sapply(quantile_list, quantile)

################
covariates <- eSVD_obj[["covariates"]]

load("../../out/simulation/simulation_1_nbreg_downsampled.RData")
nb_quantile_list <- vector("list", length(seq(1, 0.6, by = -.05)))
for(k in 1:length(downsample_coef_list)){

  nat_mat_other <- tcrossprod(covariates[,other_idx], downsample_coef_list[[k]][,other_idx])
  nat_mat_alt <- nat_mat_safe - nat_mat_other

  nb_quantile_list[[k]] <- sapply(1:n, function(i){
    stats::cor(true_residual_mat[i,], nat_mat_alt[i,])
  })
}

sapply(nb_quantile_list, quantile)
nb_quantile_mat <- sapply(nb_quantile_list, quantile)

###############

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
orange_col_trans <- rgb(235, 134, 47, 0.4*255, maxColorValue = 255)
gray_col <- rgb(0.5, 0.5, 0.5)
gray_col_trans <- rgb(0.5, 0.5, 0.5, 0.4)

n <- length(nb_quantile_list)
range_vec <- 4:n
png(paste0("../../out/fig/simulation/simulation_1_signaldropoff.png"),
    height = 1750, width = 1500,
    units = "px", res = 500)
par(mar = c(4,6,1,1))
plot(NA, xlim = range(range_vec), ylim = c(0.2,1),
     xaxt = "n", yaxt = "n", bty = "n", xlab = "Downsample percentage",
     ylab = "Correlation of cell's natural parameters\n after removing confounding effects")

for(j in seq(0,10,by = .5)){
  lines(x = rep(j, 2), y = c(-1e4,1e4), col = "gray", lty = 2, lwd = 1)
}
for(j in seq(0, 1, by = .1)){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}

polygon(x = c(range_vec, rev(range_vec)),
        y = c(eSVD_quantile_mat[2,range_vec], rev(eSVD_quantile_mat[4,range_vec])),
        col = orange_col_trans,
        border = NA)
polygon(x = c(range_vec, rev(range_vec)),
        y = c(nb_quantile_mat[2,range_vec], rev(nb_quantile_mat[4,range_vec])),
        col = gray_col_trans,
        border = NA)

lines(x = range_vec, y = nb_quantile_mat[3,range_vec], col = 1, lwd = 6)
lines(x = range_vec, y = nb_quantile_mat[3,range_vec], col = gray_col, lwd = 4, lty = 2)

lines(x = range_vec, y = eSVD_quantile_mat[3,range_vec], col = 1, lwd = 6)
lines(x = range_vec, y = eSVD_quantile_mat[3,range_vec], col = orange_col, lwd = 4, lty = 2)

axis(1, at = range_vec, labels = paste0(seq(30, 80, by = 10), "%"),
     cex.axis = 1, cex.lab = 1, lwd = 2)
axis(2, cex.axis = 1, cex.lab = 1, lwd = 2)
graphics.off()
