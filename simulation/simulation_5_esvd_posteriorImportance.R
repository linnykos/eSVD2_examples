rm(list=ls())
library(Seurat)

load("../eSVD2_examples/simulation/simulation_5_esvd.RData")

nat_mat1 <- tcrossprod(eSVD_obj[["fit_Second"]]$x_mat, eSVD_obj[["fit_Second"]]$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,"cc_1"], eSVD_obj[["fit_Second"]]$z_mat[,"cc_1"])
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
case_idx <- which(individual_vec %in% case_individuals)
control_idx <- which(individual_vec %in% control_individuals)
teststat_vec <- apply(mean_mat, 2, function(y){
  mean(y[case_idx]) - mean(y[control_idx])
})
max_val <- ceiling(max(abs(eSVD_obj$teststat_vec)))+1
teststat_vec <- pmax(pmin(teststat_vec, max_val), -max_val)

purple_col <- rgb(213, 65, 221, maxColorValue = 255)
green_col <- rgb(70, 179, 70, maxColorValue = 255)
gray_col <- rgb(.8, .8, .8)

col_palette <- grDevices::colorRampPalette(c(purple_col,gray_col,green_col))(100)
true_teststat_vec2 <- true_teststat_vec
idx <- which(true_teststat_vec2 < 0)
true_teststat_vec2[idx] <- true_teststat_vec2[idx]/max(abs(true_teststat_vec2[idx]))
idx <- which(true_teststat_vec2 > 0)
true_teststat_vec2[idx] <- true_teststat_vec2[idx]/max(abs(true_teststat_vec2[idx]))
transformed_vec <- abs(true_teststat_vec2) * sign(true_teststat_vec2)
range_vec <- range(transformed_vec)
col_breaks <- seq(range_vec[1], range_vec[2], length = 100)
col_vec_true <- sapply(transformed_vec, function(x){
  col_palette[which.min(abs(x-col_breaks))]
})
ordering_idx <- order(abs(true_teststat_vec), decreasing = F)
range_vec <- range(c(eSVD_obj$teststat_vec, teststat_vec))

de_idx <- which(true_fdr_vec < 0.0005)
length(de_idx)

den <- stats::density(-teststat_vec)
png("../../out/fig/simulation/simulation_5_raw-density.png",
    height = 250, width = 1250,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = -1*range_vec[c(2,1)],
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = 2, lwd = 1, lty = 2)
for(i in ordering_idx){
  rug(-teststat_vec[i], col = col_vec_true[i], lwd = 3)
}
graphics.off()


den <- stats::density(eSVD_obj$teststat_vec)
png("../../out/fig/simulation/simulation_5_posterior-density.png",
    height = 250, width = 1250,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = range_vec,
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = 2, lwd = 1.5, lty = 2)
for(i in ordering_idx){
  rug(eSVD_obj$teststat_vec[i], col = col_vec_true[i], lwd = 3)
}
graphics.off()

# de_idx <- which(abs(true_teststat_vec) >= 2.25)
# length(de_idx)

y_vec <- teststat_vec
x_vec <- eSVD_obj$teststat_vec
png("../../out/fig/simulation/simulation_5_teststatistic_scatterplot.png",
    height = 1750, width = 1250,
    units = "px", res = 500)
par(mar = c(3,3,0.5,0.5))
plot(NA, asp = T,
     xlim = range(x_vec), ylim = range(y_vec),
     xaxt = "n", yaxt = "n", bty = "n",
     cex.lab = 1.25, type = "n",
     xlab = "", ylab = "")
lines(c(-1e5,1e5), c(-1e5,1e5), col = 1, lty = 2, lwd = 2)
lines(rep(0,2), c(-1e5,1e5), col = 2, lty = 2, lwd = 2)
lines(c(-1e5,1e5), rep(0,2), col = 2, lty = 2, lwd = 2)


points(x_vec[ordering_idx], y_vec[ordering_idx], col = col_vec_true[ordering_idx], pch = 16)
points(x_vec[de_idx], y_vec[de_idx], col = "black", pch = 16, cex = 1.5)
points(x_vec[de_idx], y_vec[de_idx], col = "white", pch = 16, cex = 1.25)
points(x_vec[de_idx], y_vec[de_idx], col = col_vec_true[de_idx], pch = 16)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
