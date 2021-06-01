rm(list = ls())
load("../../../../out/writeup6/writeup6_dropseq_mouselung_naive_pvalue.RData")
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

png("../../../../out/fig/writeup6/dropseq_mouselung_volcano.png", height = 1500,
    width = 1500, res = 300, units = "px")
plot_volcano(mat, pval_vec, de_idx, main = paste0("P-value vs. standard dev.\n(Naive, ",
                                                  round(100*length(de_idx)/ncol(mat)), "% significant)"),
             bool_iqr = F)
graphics.off()

png("../../../../out/fig/writeup6/dropseq_mouselung_scatter.png", height = 1500,
    width = 1500, res = 300, units = "px")
plot_sd_scatter(mat, lung@meta.data["celltype"][,1], de_idx,
                main = "Within/Between standard dev.\n(Naive)",
                bool_iqr = F)
graphics.off()

png("../../../../out/fig/writeup6/dropseq_mouselung_pval_hist.png", height = 1500,
    width = 1500, res = 300, units = "px")
graphics::hist(as.numeric(pval_obj$pval_mat), col = "gray", xlab = "Individual pvalues",
               main = "Histogram of pvalues (Naive)")
graphics.off()
