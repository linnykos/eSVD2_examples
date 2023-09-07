rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/habermann_T_esvd.RData")

y_vec <- eSVD_obj$fit_Second$z_mat[,"Log_UMI"]
x_vec <- -log10(eSVD_obj$fit_Second$nuisance_vec)

cor_val <- stats::cor(x_vec, y_vec)
xlim <- stats::quantile(x_vec, probs = c(0.01, 1))
ylim <- stats::quantile(y_vec, probs = c(0.005, 0.995))

p <- length(y_vec)
cex_vec <- rep(1, p)
png("../../../out/fig/main/habermann_T_scatterplot_logumi-nuisance.png",
    height = 2250, width = 2000,
    units = "px", res = 500)
par(mar = c(4,4,4,0.1))
plot(x_vec, y_vec,
     col = rgb(0.5, 0.5, 0.5, 0.2), pch = 16,
     cex = cex_vec,
     xaxt = "n", yaxt = "n", bty = "n", 
     xlab = "", ylab = "",
     xlim = xlim, ylim = ylim,
     main = paste0("Correlation: ", round(cor_val, 2)))
lines(x = c(-1e4,1e4), y = rep(1,2), col = "red", lwd = 2, lty = 2)
lines(x = rep(0,2), y = c(-1e4,1e4), col = "red", lwd = 2, lty = 2)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
