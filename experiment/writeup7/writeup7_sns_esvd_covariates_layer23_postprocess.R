rm(list=ls())

library(Seurat)

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_2.RData")

de_genes <- read.csv("../../../../data/sns_autism/de_genes.txt", header = F)
de_genes2 <- sort(unique(de_genes[,1]))
length(de_genes2)

idx <- which(colnames(mat) %in% de_genes2)
length(idx)

quantile(esvd_res$b_mat[idx,3]) # doesn't look particularly promising
quantile(esvd_res$b_mat[,3])

png("../../../../out/fig/writeup7/sns_esvd_covariates_layer23_asd_hist.png",
    width = 1800, height = 1500, units = "px", res = 300)
hist(esvd_res$b_mat[,3], main = "Human brain (SNS, with covariates)\neSVD, ASD coef. histogram",
     col = "gray", xlab = "ASD coefficient", breaks = 50)
rug(esvd_res$b_mat[idx,3], col = "red")
graphics.off()

###################

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates[,-2], esvd_res$b_mat[,-2])
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

idx <- which(metadata$diagnosis == "ASD")
pval_vec <- sapply(1:ncol(pred_mat), function(j){
  if(j %% 1000 == 0) print(paste0(j, " out of ", ncol(pred_mat)))

  stats::wilcox.test(x = pred_mat[idx,j],
                     y = pred_mat[-idx,j])$p.value
})

idx <- which(colnames(mat) %in% de_genes2)
png("../../../../out/fig/writeup7/sns_esvd_covariates_layer23_asd_wilcox.png",
    width = 3000, height = 1500, units = "px", res = 300)
par(mfrow = c(1,2))
hist(pval_vec, main = "Human brain (SNS, with covariates)\neSVD, Wilcoxon p-value",
     col = "gray", xlab = "ASD coefficient", breaks = 50)
rug(pval_vec[idx], col = "red")

p <- length(pval_vec)
ord <- order(pval_vec, decreasing = F)
col_vec <- rep(rgb(0.5, 0.5, 0.5), p)
col_vec[idx] <- rgb(1, 0, 0)
plot(pval_vec[ord], c(1:p)/p, pch = 16, col = col_vec[ord],
     main = "Human brain (SNS, with covariates)\neSVD, Wilcoxon p-value",
     xlab = "(Sorted) p-value", ylab = "rank", xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)
graphics.off()

################################

set.seed(10)
idx2 <- sample(1:n, length(idx))
pval_vec2 <- sapply(1:ncol(pred_mat), function(j){
  if(j %% 1000 == 0) print(paste0(j, " out of ", ncol(pred_mat)))

  set.seed(10)
  stats::wilcox.test(x = pred_mat[idx2,j],
                     y = pred_mat[-idx2,j])$p.value
})


idx <- which(colnames(mat) %in% de_genes2)
png("../../../../out/fig/writeup7/sns_esvd_covariates_layer23_asd_wilcox_reshuffled.png",
    width = 3000, height = 1500, units = "px", res = 300)
par(mfrow = c(1,2))
hist(pval_vec2, main = "Human brain (SNS, with covariates)\neSVD, Shuffled Wilcoxon p-value",
     col = "gray", xlab = "ASD coefficient", breaks = 50)
rug(pval_vec2[idx], col = "red")

p <- length(pval_vec)
ord <- order(pval_vec2, decreasing = F)
col_vec <- rep(rgb(0.5, 0.5, 0.5), p)
col_vec[idx] <- rgb(1, 0, 0)
plot(pval_vec2[ord], c(1:p)/p, pch = 16, col = col_vec[ord],
     main = "Human brain (SNS, with covariates)\neSVD, Shuffled Wilcoxon p-value",
     xlab = "(Sorted) p-value", ylab = "rank", xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1), col = "red", lwd = 2, lty = 2)
graphics.off()
