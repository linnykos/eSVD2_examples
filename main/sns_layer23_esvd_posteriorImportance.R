rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L2/3"),]
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

hk_idx <- which(names(eSVD_obj$teststat_vec) %in% c(hk_genes, cycling_genes))
de_idx <- which(names(eSVD_obj$teststat_vec) %in% de_gene_specific)
other_idx <- which(names(eSVD_obj$teststat_vec) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj$teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

###############3

cc_var <- eSVD_obj$param$init_case_control_variable
nuisance_vec <- eSVD_obj[["fit_Second"]]$nuisance_vec
nat_mat1 <- tcrossprod(eSVD_obj[["fit_Second"]]$x_mat, eSVD_obj[["fit_Second"]]$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,cc_var], eSVD_obj[["fit_Second"]]$z_mat[,cc_var])
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
case_idx <- which(eSVD_obj$case_control == 1)
control_idx <- which(eSVD_obj$case_control == 0)
teststat_vec <- apply(mean_mat, 2, function(y){
  mean(y[case_idx]) - mean(y[control_idx])
})
teststat_vec <- pmax(pmin(teststat_vec, 5), -5)

# var_mat <- sweep(x = mean_mat, MARGIN = 2, STATS = nuisance_vec, FUN = "/")
#
# case_individuals <- unique(eSVD_obj$individual[which(eSVD_obj$case_control == 1)])
# control_individuals <- unique(eSVD_obj$individual[which(eSVD_obj$case_control == 0)])
# individual_vec <- eSVD_obj$individual
#
# res <- eSVD2:::compute_test_statistic.default(
#   input_obj = mean_mat,
#   posterior_var_mat = var_mat,
#   case_individuals = case_individuals,
#   control_individuals = control_individuals,
#   individual_vec = individual_vec
# )
# teststat_vec <- res$teststat_vec

#########

range_vec <- range(c(eSVD_obj$teststat_vec, teststat_vec))

den <- stats::density(-teststat_vec)
png("../../../out/fig/main/sns_layer23_posteriorImportance_raw-density.png",
    height = 250, width = 1250,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = -1*range_vec[c(2,1)],
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = 2, lwd = 1, lty = 2)
for(i in shuf_idx){
  rug(-teststat_vec[i], col = col_vec[i], lwd = 3)
}
graphics.off()


den <- stats::density(eSVD_obj$teststat_vec)
png("../../../out/fig/main/sns_layer23_posteriorImportance_posterior-density.png",
    height = 250, width = 1250,
    units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0.1))
plot(den, main = "", xlab = "", ylab = "", xlim = range_vec,
     xaxt = "n", yaxt = "n", bty ="n")
polygon(den, col="gray")
lines(x = rep(0,2), y = c(-1e4,1e4), col = 2, lwd = 1.5, lty = 2)
for(i in shuf_idx){
  rug(-teststat_vec[i], col = col_vec[i], lwd = 3)
}
graphics.off()

# de_idx <- which(abs(true_teststat_vec) >= 2.25)
# length(de_idx)

y_vec <- teststat_vec
x_vec <- eSVD_obj$teststat_vec
png("../../../out/fig/main/sns_layer23_posteriorImportance_teststatistic_scatterplot.png",
    height = 1750, width = 1250,
    units = "px", res = 500)
par(mar = c(3,3,0.5,0.5))
plot(NA, asp = T,
     xlim = range(x_vec), ylim = range(y_vec),
     xaxt = "n", yaxt = "n", bty = "n",
     cex.lab = 1.25, type = "n",
     xlab = "", ylab = "")
lines(c(-1e5,1e5), c(-1e5,1e5), col = 1, lty = 2, lwd = 2)
lines(rep(0,2), c(-1e5,1e5), col = 2, lty = 2, lwd = 2)
lines(c(-1e5,1e5), rep(0,2), col = 2, lty = 2, lwd = 2)

points(x_vec, y_vec, col = col_vec, pch = 16)
points(x_vec[shuf_idx], y_vec[shuf_idx], col = "black", pch = 16, cex = 1.5)
points(x_vec[shuf_idx], y_vec[shuf_idx], col = "white", pch = 16, cex = 1.25)
points(x_vec[shuf_idx], y_vec[shuf_idx], col = col_vec[shuf_idx], pch = 16)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
