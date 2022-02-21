rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")
load("../../../../out/Writeup10/Writeup10_sns_invip_ideas.RData")

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
load("../../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "IN-VIP"),]
de_gene_specific <- tmp[,"Gene name"]
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list[-1], function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
de_genes <- de_genes[!de_genes %in% de_gene_specific]
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(colnames(mat) %in% c(hk_genes, cycling_genes))
de_idx <- which(colnames(mat) %in% de_gene_specific)
other_idx <- which(colnames(mat) %in% c(sfari_genes, de_genes))

#####################

pval_res[is.na(pval_res)] <- 0
pval_res[pval_res == 0] <- min(pval_res[pval_res != 0])/2
logpval <- -log10(pval_res)
fdr_idx <- which(stats::p.adjust(pval_res, method = "BH") <= 0.05)

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), ncol(mat))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

metadata <- sns@meta.data
indiv_list <- lapply(unique(metadata$individual), function(indiv){
  which(metadata$individual == indiv)
})
names(indiv_list) <- unique(metadata$individual)
case_individuals <- as.character(unique(metadata[which(metadata$diagnosis == "ASD"),"individual"]))
control_individuals <- as.character(unique(metadata[which(metadata$diagnosis == "Control"),"individual"]))
mat_avg <- t(sapply(indiv_list, function(idx_vec){
  matrixStats::colMeans2(mat[idx_vec,])
}))
rownames(mat_avg) <- names(indiv_list)
x_vec <- sapply(1:ncol(mat_avg), function(j){
  # log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))
  log2(mean(mat_avg[case_individuals,j])+1) - log2(mean(mat_avg[control_individuals,j])+1)
})
z_score <- sapply(1:length(pval_res), function(i){
  tmp <- abs(stats::qnorm(pval_res[i]/2, mean = 0, sd = 1))
  tmp
})


sparisity_vec <- sapply(1:ncol(mat), function(j){
  length(which(mat[,j] == 0))/nrow(mat)
})

quantile(logpval)
y_max <- max(logpval)
x_max <- ceiling(max(abs(x_vec)))
png("../../../../out/fig/Writeup10/sns_invip_ideas.png",
    height = 1000, width = 2500,
    units = "px", res = 300)
par(mfrow = c(1,3))
# z-scores
max_val <- max(abs(z_score))
break_vec <- seq(-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -0.05; break_vec[length(break_vec)] <- max_val+0.05
hist(z_score, breaks = break_vec,
     xlim = c(-0.05, max_val),
     main = "Histogram of test statistic\n(IDEAS)",
     xlab = "One-sided Z-score", ylab = "Frequency", freq = T)
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(z_score[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Other interest gene", "Housekeeping gene"),
       fill = c(2,4,3))

# volcano plot
plot(NA, xlim = c(-x_max, x_max), ylim = range(0, y_max), bty = "n",
     main = "Volcano plot for IN-VIP",
     xlab = "Log2 fold change (i.e., Log2 mean difference)", ylab = "-Log10(P value) via IDEAS")
for(x in seq(-x_max, x_max,by=0.5)){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
lines(rep(0,2), c(-1e5,1e5), col = "gray")
for(y in seq(0,y_max,by=2)){
  lines(c(-1e5,1e5), rep(y,2), lty = 2, col = "gray", lwd = 0.5)
}
points(x = x_vec[-unique(c(hk_idx,de_idx))],
       y = logpval[-unique(c(hk_idx,de_idx))],
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = x_vec[shuf_idx],
       y = logpval[shuf_idx],
       pch = 16, col = "white", cex = 1.5)
points(x = x_vec[shuf_idx],
       y = logpval[shuf_idx],
       pch = 16, col = col_vec[shuf_idx])
lines(x = c(-2*x_max, 2*x_max), y = rep(min(logpval[fdr_idx]),2),
      col = 2, lwd = 2, lty = 2)

#sparisty
xlim <- range(sparisity_vec)
ylim <- range(logpval)
plot(NA, xlim = xlim, ylim = ylim,
     xlab = "% of observed 0's in gene",
     ylab = "-Log10(P value)",
     main = "P-value vs. sparsity")
points(x = sparisity_vec[-unique(c(hk_idx,de_idx))],
       y = logpval[-unique(c(hk_idx,de_idx))],
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = sparisity_vec[shuf_idx],
       y = logpval[shuf_idx],
       pch = 16, col = "white", cex = 1.5)
points(x = sparisity_vec[shuf_idx],
       y = logpval[shuf_idx],
       pch = 16, col = col_vec[shuf_idx])
lines(x = c(-2*max(abs(sparisity_vec)), 2*max(abs(sparisity_vec))),
      y = rep(min(logpval[fdr_idx]),2),
      col = 2, lwd = 2, lty = 2)
graphics.off()


