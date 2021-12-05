rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson4.RData")

mat[mat == min(mat)] <- 0
zero_prop <- apply(mat, 2, function(x){length(which(x == 0))/length(x)})
png("../../out/fig/writeup8e/ppt_library_size.png", height = 1100, width = 1700,
    res = 300, units = "px")
hist(zero_prop, breaks = seq(0, 1, by = 0.1),
     main = "Proportion of 0's for each gene (i.e., column)\n across all cells (i.e., rows)",
     xlab = "Proportion of 0's across all cells", ylab = "Number of genes")
graphics.off()

idx <- which.min(abs(zero_prop - 0.95))
png("../../out/fig/writeup8e/ppt_sparse_gene.png",
    height = 900, width = 1400,
    res = 300, units = "px")
par(mar = c(4,4,3,0.5))
hist(mat[,idx], breaks = seq(-.5, max(mat[,idx])+.5),
     main = "Example of sparse gene", ylab = "Number of cells", xlab = "Count")
rug(jitter(mat[,idx]), col = 2)
graphics.off()

idx <- which.min(abs(zero_prop - 0.1))
png("../../out/fig/writeup8e/ppt_dense_gene.png",
    height = 900, width = 1400,
    res = 300, units = "px")
par(mar = c(4,4,3,0.5))
hist(mat[,idx], breaks = seq(-.5, max(mat[,idx])+.5, by = 1),
     main = "Example of dense gene", ylab = "Number of cells", xlab = "Count")
rug(jitter(mat[,idx]), col = 2)
graphics.off()
