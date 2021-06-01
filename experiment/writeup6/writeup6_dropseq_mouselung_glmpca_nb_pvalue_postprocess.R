rm(list=ls())
load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb_pvalue.RData")
source("plotting_func.R")

pval_obj <- pval_res
uniq_val <- sort(unique(as.numeric(pval_obj$lookup_mat[,1])))
p <- ncol(pval_obj$pval_mat)
quantile(pval_obj$pval_mat)

# intersection-union
max_pval_mat <- t(sapply(uniq_val, function(i){
  idx <- which(pval_obj$lookup[,1] == i)
  apply(pval_obj$pval_mat[idx,], 2, max)
}))

max_pval_vec <- apply(max_pval_mat, 2, function(x){min(stats::p.adjust(x, method = "bonferroni"))})

adjust_pval_vec <- stats::p.adjust(max_pval_vec, method = "bonferroni")
pval_thres <- 5*10^(-8)
de_idx <- which(adjust_pval_vec <= pval_thres)
length(de_idx)

pval_vec <- apply(max_pval_mat, 2, min)

#########

png("../../../../out/fig/writeup6/dropseq_mouselung_glmpca_nb_volcano.png", height = 1500,
    width = 1500, res = 300, units = "px")
plot_volcano(pred_mat_unnorm, pval_vec, de_idx, main = paste0("P-value vs. standard dev.\n(GLM-PCA, NB, ",
                                                  round(100*length(de_idx)/ncol(pred_mat_unnorm)), "% significant)"))
graphics.off()

png("../../../../out/fig/writeup6/dropseq_mouselung_glmpca_nb_scatter.png", height = 1500,
    width = 1500, res = 300, units = "px")
plot_sd_scatter(pred_mat_unnorm, lung@meta.data["celltype"][,1], de_idx,
                main = "Within/Between standard dev.\n(GLM-PCA, NB)")
graphics.off()

png("../../../../out/fig/writeup6/dropseq_mouselung_glmpca_nb_pval_hist.png", height = 1500,
    width = 1500, res = 300, units = "px")
graphics::hist(as.numeric(pval_obj$pval_mat), col = "gray", xlab = "Individual pvalues",
               main = "Histogram of pvalues (GLM-PCA, NB)")
graphics.off()

png("../../../../out/fig/writeup6/dropseq_mouselung_glmpca_nb_prediction.png", height = 1500,
    width = 1500, res = 300, units = "px")
idx <- which(as.matrix(mat) != 0)
graphics::plot(y = as.matrix(mat)[idx], x = pred_mat[idx], pch = 16, asp = T,
               col = rgb(0.5, 0.5, 0.5, 0.5),
               xlab = "Predicted value", ylab = "Observed value",
               main = "Raw observed vs. predicted\n(GLM-PCA, NB)",
               xlim = c(1,2000), ylim = c(1,2000))
graphics.off()

png("../../../../out/fig/writeup6/dropseq_mouselung_glmpca_nb_prediction_rescaled.png", height = 1500,
    width = 1500, res = 300, units = "px")
idx <- which(as.matrix(mat) != 0)
vec <- as.matrix(mat)[idx]
med1 <- median(vec); med2 <- median(pred_mat[idx])
graphics::plot(y = as.matrix(mat)[idx]/med1, x = pred_mat[idx]/med2, pch = 16, asp = T,
               col = rgb(0.5, 0.5, 0.5, 0.5),
               xlab = "Rescaled predicted value", ylab = "Rescaled observed value",
               main = "Raw observed vs. predicted\n(GLM-PCA, NB, rescaled)",
               xlim = range(vec), ylim = range(vec))
graphics.off()

