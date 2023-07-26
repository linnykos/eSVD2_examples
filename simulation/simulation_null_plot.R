rm(list=ls())
load("../eSVD2_examples/simulation/simulation_null_esvd.RData")

png(paste0("../../out/fig/simulation/simulation_null_esvd.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.5,0.5))
plot(NA,
     xlim = c(0,1), ylim = c(0,1),
     asp = T,
     xlab = "Observed quantile",
     ylab = "Theoretical quantile",
     xaxt = "n", yaxt = "n", bty = "n")
for(x in seq(0,1,length.out=11)){
  lines(rep(x,2), c(0,1), lwd = 0.5, lty = 3, col = "gray")
}
for(y in seq(0,1,length.out=11)){
  lines(c(0,1), rep(y,2), lwd = 0.5, lty = 3, col = "gray")
}
points(sort(multtest_res$pvalue_vec[-c(1:10)]),
       seq(0,1,length.out = length(multtest_res$pvalue_vec[-c(1:10)])),
       pch = 16, col = "black")
lines(c(0,1), c(0,1), col = 2, lty = 2, lwd = 2)
axis(1); axis(2)
graphics.off()

#######################

load("../eSVD2_examples/simulation/simulation_null_deseq2.RData")

png(paste0("../../out/fig/simulation/simulation_null_deseq2.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.5,0.5))
plot(NA,
     xlim = c(0,1), ylim = c(0,1),
     asp = T,
     xlab = "Observed quantile",
     ylab = "Theoretical quantile",
     xaxt = "n", yaxt = "n", bty = "n")
for(x in seq(0,1,length.out=11)){
  lines(rep(x,2), c(0,1), lwd = 0.5, lty = 3, col = "gray")
}
for(y in seq(0,1,length.out=11)){
  lines(c(0,1), rep(y,2), lwd = 0.5, lty = 3, col = "gray")
}
points(sort(pvalue_vec[-c(1:10)]),
       seq(0,1,length.out = length(pvalue_vec[-c(1:10)])),
       pch = 16, col = "black")
lines(c(0,1), c(0,1), col = 2, lty = 2, lwd = 2)
axis(1); axis(2)
graphics.off()

#########################

load("../eSVD2_examples/simulation/simulation_null_mast.RData")

png(paste0("../../out/fig/simulation/simulation_null_mast.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.5,0.5))
plot(NA,
     xlim = c(0,1), ylim = c(0,1),
     asp = T,
     xlab = "Observed quantile",
     ylab = "Theoretical quantile",
     xaxt = "n", yaxt = "n", bty = "n")
for(x in seq(0,1,length.out=11)){
  lines(rep(x,2), c(0,1), lwd = 0.5, lty = 3, col = "gray")
}
for(y in seq(0,1,length.out=11)){
  lines(c(0,1), rep(y,2), lwd = 0.5, lty = 3, col = "gray")
}
points(sort(pvalue_vec[-c(1:10)]),
       seq(0,1,length.out = length(pvalue_vec[-c(1:10)])),
       pch = 16, col = "black")
lines(c(0,1), c(0,1), col = 2, lty = 2, lwd = 2)
axis(1); axis(2)
graphics.off()

#################

load("../eSVD2_examples/simulation/simulation_null_sctransform.RData")

png(paste0("../../out/fig/simulation/simulation_null_sctransform.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.5,0.5))
plot(NA,
     xlim = c(0,1), ylim = c(0,1),
     asp = T,
     xlab = "Observed quantile",
     ylab = "Theoretical quantile",
     xaxt = "n", yaxt = "n", bty = "n")
for(x in seq(0,1,length.out=11)){
  lines(rep(x,2), c(0,1), lwd = 0.5, lty = 3, col = "gray")
}
for(y in seq(0,1,length.out=11)){
  lines(c(0,1), rep(y,2), lwd = 0.5, lty = 3, col = "gray")
}
points(sort(pvalue_vec[-c(1:10)]),
       seq(0,1,length.out = length(pvalue_vec[-c(1:10)])),
       pch = 16, col = "black")
lines(c(0,1), c(0,1), col = 2, lty = 2, lwd = 2)
axis(1); axis(2)
graphics.off()

####################

load("../eSVD2_examples/simulation/simulation_null_glmpca.RData")

png(paste0("../../out/fig/simulation/simulation_null_glmpca.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(4,4,0.5,0.5))
plot(NA,
     xlim = c(0,1), ylim = c(0,1),
     asp = T,
     xlab = "Observed quantile",
     ylab = "Theoretical quantile",
     xaxt = "n", yaxt = "n", bty = "n")
for(x in seq(0,1,length.out=11)){
  lines(rep(x,2), c(0,1), lwd = 0.5, lty = 3, col = "gray")
}
for(y in seq(0,1,length.out=11)){
  lines(c(0,1), rep(y,2), lwd = 0.5, lty = 3, col = "gray")
}
points(sort(pvalue_vec[-c(1:10)]),
       seq(0,1,length.out = length(pvalue_vec[-c(1:10)])),
       pch = 16, col = "black")
lines(c(0,1), c(0,1), col = 2, lty = 2, lwd = 2)
axis(1); axis(2)
graphics.off()
