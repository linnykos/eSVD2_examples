rm(list=ls())
load("../../../../out/writeup6/writeup6_citeseq_bm_naive_pvalue.RData")

pval_obj <- pval_res

uniq_val <- sort(unique(as.numeric(pval_obj$lookup_mat[,1])))
p <- ncol(pval_obj$pval_mat)
quantile(pval_obj$pval_mat)

# intersection-union
max_pval_mat <- t(sapply(uniq_val, function(i){
  idx <- which(pval_obj$lookup[,1] == i)
  apply(pval_obj$pval_mat[idx,], 2, max)
}))

# max_pval_vec <- apply(max_pval_mat, 2, min)
# max_pval_vec <- sapply(max_pval_vec, function(x){stats::pbeta(x, length(uniq_val)-1,1)})
# max_pval_vec <- pmin(max_pval_vec*length(uniq_val), 1)
max_pval_vec <- pmin(apply(max_pval_mat*nrow(max_pval_mat), 2, min),1)

zz <- p.adjust(max_pval_vec, method = "bonferroni")
length(which(zz <= 5*10^(-8)))
