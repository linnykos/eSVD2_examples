rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer56cc_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
teststat_vec <- eSVD_obj$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
fdr_vec <- locfdr_res$fdr
names(fdr_vec) <- names(gaussian_teststat)
null_mean <- locfdr_res$fp0["mlest", "delta"]
null_sd <- locfdr_res$fp0["mlest", "sigma"]
logpvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec <- -(logpvalue_vec/log(10) + log10(2))
logpvalue_vec <- pmin(logpvalue_vec, 15)

selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
length(selected_genes)

############

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L5/6-CC"),]
de_gene_specific <- tmp[,"Gene name"]
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list[-1], function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
de_genes <- de_genes[!de_genes %in% de_gene_specific]
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.005),"external_gene_name"]

gene_lists <- list(sfari_genes = sfari_genes,
                   bulk_de_genes = bulk_de_genes,
                   de_gene_specific = de_gene_specific,
                   de_genes = de_genes,
                   hk_genes = hk_genes)

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
sapply(gene_lists, function(vec){
  m <- length(intersect(vec, names(gaussian_teststat)))
  n <- length(gaussian_teststat) - m
  k <- length(selected_genes)
  x <- length(intersect(selected_genes, vec))
  # m;k;x
  expected <- k*(m/ length(gaussian_teststat))
  fisher <- sum(sapply(x:k, function(i){
    stats::dhyper(x = i, m = m, n = n, k = k, log = F)
  }))

  c(observed = x, expected = expected, pvalue = fisher)
})

#######

idx <- which(names(eSVD_obj$teststat_vec) %in% selected_genes)
sfari_idx <- which(names(eSVD_obj$teststat_vec) %in% sfari_genes)
hk_idx <- which(names(eSVD_obj$teststat_vec) %in% hk_genes)

x_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
xlim <- quantile(x_vec, probs = c(0.01, 1))
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- logpvalue_vec
ylim <- range(y_vec)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

# adjustment of selected genes for visual clarity
min_pthres <- min(y_vec[selected_genes])
selected_genes2 <- names(y_vec)[y_vec >= min_pthres]
idx <- which(names(y_vec) %in% selected_genes2)
sfari_idx <- which(names(y_vec) %in% sfari_genes)
bulk_idx <- which(names(y_vec) %in% bulk_de_genes)
hk_idx <- which(names(y_vec) %in% hk_genes)

png("../../../out/fig/main/sns_layer56cc_volcano.png",
    height = 3500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,1.5))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in seq(0,15,by = .5)){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}
lines(x = c(-1e4,1e4), y = rep(min_pthres, 2), col = orange_col, lty = 2, lwd = 2)

points(x = x_vec, y = y_vec,
       col = rgb(0.6, 0.6, 0.6, 0.3), pch = 16)
points(x = x_vec[idx], y = y_vec[idx],
       col = orange_col, pch = 16, cex = 1.5)

# plot non-overlapping genes
points(x = x_vec[setdiff(bulk_idx, idx)], y = y_vec[setdiff(bulk_idx, idx)],
       col = blue_col, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[setdiff(sfari_idx, idx)], y = y_vec[setdiff(sfari_idx, idx)],
       col = purple_col, pch = 16, cex = 1, lwd = 2)

# plot housekeeping
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = green_col_trans, pch = 16, cex = 1)

# circle overlapping gens
points(x = x_vec[intersect(idx, sfari_idx)], y = y_vec[intersect(idx, sfari_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx, sfari_idx)], y = y_vec[intersect(idx, sfari_idx)],
       col = purple_col, pch = 1, cex = 2, lwd = 2)
points(x = x_vec[intersect(idx, bulk_idx)], y = y_vec[intersect(idx, bulk_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx, bulk_idx)], y = y_vec[intersect(idx, bulk_idx)],
       col = blue_col, pch = 1, cex = 2, lwd = 2)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

###############

write_genes <- function(vec,
                        file){
  fileConn <- file(file)
  writeLines(vec, fileConn)
  close(fileConn)
}

write_genes(selected_genes, file = "../../../out/main/sns_layer56cc_esvd_de-genes.csv")
