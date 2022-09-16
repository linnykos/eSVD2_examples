rm(list=ls())
library(Seurat)

load("../eSVD2_examples/simulation/simulation_2_esvd.RData")

#########################

# col_palette <- c("none" = rgb(0.5, 0.5, 0.5),
#                  "strong-negative" = rgb(0.75, 0, 0),
#                  "strong-positive" = rgb(0, 0.75, 0),
#                  "weak-negative" = rgb(1, 0.5, 0.9),
#                  "weak-positive" = rgb(0.5, 1, 0.9))
# col_vec <- plyr::mapvalues(gene_labeling2, from = names(col_palette), to = col_palette)
# plot(eSVD_obj$teststat_vec, col = col_vec, pch = 16)

#########################

nat_mat1 <- tcrossprod(eSVD_obj[["fit_Second"]]$x_mat, eSVD_obj[["fit_Second"]]$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,"cc_1"], eSVD_obj[["fit_Second"]]$z_mat[,"cc_1"])
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

var_mat <- sweep(x = mean_mat, MARGIN = 2, STATS = nuisance_vec, FUN = "/")

res <- eSVD2:::compute_test_statistic.default(
  input_obj = mean_mat,
  posterior_var_mat = var_mat,
  case_individuals = case_individuals,
  control_individuals = control_individuals,
  individual_vec = individual_vec
)
teststat_vec <- res$teststat_vec

purple_col <- rgb(213, 65, 221, maxColorValue = 255)
green_col <- rgb(70, 179, 70, maxColorValue = 255)
gray_col <- rgb(.8, .8, .8)

col_palette <- grDevices::colorRampPalette(c(purple_col,gray_col,green_col))(100)
# rank_vec <- rank(true_teststat_vec)-length(true_teststat_vec)/2
# rank_vec <- rank_vec/max(abs(rank_vec))
# rank_vec <- sign(rank_vec)*abs(rank_vec)^25
# range_vec <- range(rank_vec)
# col_breaks <- seq(range_vec[1], range_vec[2], length = 100)
# col_vec_true <- sapply(rank_vec, function(x){
#   col_palette[which.min(abs(x-col_breaks))]
# })
true_teststat_vec2 <- true_teststat_vec
idx <- which(true_teststat_vec2 < 0)
true_teststat_vec2[idx] <- true_teststat_vec2[idx]/max(abs(true_teststat_vec2[idx]))
idx <- which(true_teststat_vec2 > 0)
true_teststat_vec2[idx] <- true_teststat_vec2[idx]/max(abs(true_teststat_vec2[idx]))
transformed_vec <- abs(true_teststat_vec2)^3 * sign(true_teststat_vec2)
range_vec <- range(transformed_vec)
col_breaks <- seq(range_vec[1], range_vec[2], length = 100)
col_vec_true <- sapply(transformed_vec, function(x){
  col_palette[which.min(abs(x-col_breaks))]
})
ordering_idx <- order(abs(true_teststat_vec), decreasing = F)
# plot(true_teststat_vec, col = col_vec_true, pch = 16)
range_vec <- range(c(eSVD_obj$teststat_vec, teststat_vec))

png("../../out/fig/simulation/simulation_2_raw-density.png",
    height = 250, width = 1250,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
den <- stats::density(-teststat_vec)
plot(den, main = "", xlab = "", ylab = "", xlim = -1*range_vec[c(2,1)],
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = 2, lwd = 1, lty = 2)
for(i in ordering_idx){
  rug(-teststat_vec[i], col = col_vec_true[i], lwd = 3)
}
graphics.off()


png("../../out/fig/simulation/simulation_2_posterior-density.png",
    height = 250, width = 1250,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
den <- stats::density(eSVD_obj$teststat_vec)
plot(den, main = "", xlab = "", ylab = "", xlim = range_vec,
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = 2, lwd = 1.5, lty = 2)
for(i in ordering_idx){
  rug(eSVD_obj$teststat_vec[i], col = col_vec_true[i], lwd = 3)
}
graphics.off()

de_idx <- which(abs(true_teststat_vec) >= 2.25)
length(de_idx)

png("../../out/fig/simulation/simulation_2_teststatistic_scatterplot.png",
    height = 1750, width = 1250,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
y_vec <- teststat_vec
x_vec <- eSVD_obj$teststat_vec
plot(NA, asp = T,
     xlim = range(x_vec), ylim = range(y_vec),
     xaxt = "n", yaxt = "n", bty = "n",
     cex.lab = 1.25, type = "n",
     xlab = "", ylab = "")
lines(c(-1e5,1e5), c(-1e5,1e5), col = 1, lty = 2, lwd = 2)
lines(rep(0,2), c(-1e5,1e5), col = 2, lty = 2, lwd = 2)
lines(c(-1e5,1e5), rep(0,2), col = 2, lty = 2, lwd = 2)
points(x_vec[ordering_idx], y_vec[ordering_idx], col = col_vec_true[ordering_idx], pch = 16)
# points(x_vec[de_idx], y_vec[de_idx], col = 1, pch = 1, cex = 2)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
