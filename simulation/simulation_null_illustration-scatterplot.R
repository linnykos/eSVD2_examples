rm(list=ls())
library(Seurat)
library(eSVD2)

load("../eSVD2_examples/simulation/simulation_null.RData")

mat <- t(as.matrix(seurat_obj[["RNA"]]@counts))

case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(10)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(10)
transparent_gray <- rgb(0.5,0.5,0.5,0.4)
two_letters <- substr(transparent_gray, start = 8, stop = 9)
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

col_vec <- c(rep(case_color_trans_palette, each = 100),
             rep(control_color_trans_palette, each = 100))

n <- nrow(mat)
png(paste0("../../out/fig/simulation/simulation_null_scatterplot_trueDE.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.75,0.75))
plot(mat[,4] + stats::runif(n, min = -0.2, max = 0.2),
     ylim = c(0,10),
     pch = 16, col = col_vec, cex = 0.5,
     xlab = "Cells (100 per individual)",
     ylab = "Observed count")
graphics.off()

png(paste0("../../out/fig/simulation/simulation_null_scatterplot_null-large-var.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.75,0.75))
plot(mat[,13] + stats::runif(n, min = -0.2, max = 0.2),
     pch = 16, col = col_vec,  cex = 0.5,
     ylim = c(0,10),
     xlab = "Cells (100 per individual)",
     ylab = "Observed count")
graphics.off()

png(paste0("../../out/fig/simulation/simulation_null_scatterplot_null-interleaved.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.75,0.75))
plot(mat[,12] + stats::runif(n, min = -0.2, max = 0.2),
     pch = 16, col = col_vec, cex = 0.5,
     ylim = c(0,10),
     xlab = "Cells (100 per individual)",
     ylab = "Observed count")
graphics.off()
