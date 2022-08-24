rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

load("../../../../data/sns_autism/velmeshev_genes.RData")
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
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
deg_df <- readxl::read_xlsx(
  path = "../../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.005),"external_gene_name"]

hk_idx <- which(names(eSVD_obj$teststat_vec) %in% c(hk_genes, cycling_genes))
de_idx <- which(names(eSVD_obj$teststat_vec) %in% de_gene_specific)
other_idx <- which(names(eSVD_obj$teststat_vec) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj$teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

###########################

x_mat <- eSVD_obj$fit_Second$x_mat
y_mat <- eSVD_obj$fit_Second$y_mat
z_mat <-  eSVD_obj$fit_Second$z_mat
covariates <- eSVD_obj$covariates

input_obj <- eSVD_obj
metadata <- sns@meta.data
covariate_individual <- "individual"
case_control_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
covariates <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
cc_vec <- covariates[,case_control_variable]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])
individual_vec <- metadata[,covariate_individual]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])
tmp <- eSVD2:::.determine_individual_indices(case_individuals = case_individuals,
                                             control_individuals = control_individuals,
                                             covariate_individual = covariate_individual,
                                             metadata = metadata)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- eSVD2:::.construct_averaging_matrix(idx_list = all_indiv_idx,
                                               n = nrow(x_mat))
case_row_idx <- 1:length(case_individuals)
control_row_idx <- (length(case_individuals)+1):nrow(avg_mat)

# attempt 1: only residuals
nat_mat1 <- exp(tcrossprod(x_mat, y_mat))
nat_mat1 <- avg_mat %*% nat_mat1
case_gaussian_mean <- Matrix::colMeans(nat_mat1[case_row_idx,,drop = F])
control_gaussian_mean <- Matrix::colMeans(nat_mat1[control_row_idx,,drop = F])
test_vec1 <- case_gaussian_mean - control_gaussian_mean
quantile(test_vec1)
length(which(test_vec1 < 0))/length(test_vec1)

teststat_vec <- pmax(pmin(test_vec1, 30), -30)
max_val <- max(abs(teststat_vec))
png("../../../../out/fig/Writeup12/Writeup12_sns_layer23_teststat1.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
main_vec <- c("Published DE", "Other interest\n(SFARI, DE other region)", "Housekeeping+Cell-cycle")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency", freq = T)
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()

###############################################

# attempt 2: all individual-residuals
nat_mat1 <- exp(tcrossprod(x_mat, y_mat))
col_idx1 <- which(colnames(covariates) %in% c("age", "Intercept", "sex_F", "diagnosis_ASD"))
col_idx2 <- grep("^individual", colnames(covariates))
nat_mat2 <- exp(tcrossprod(covariates[,c(col_idx1, col_idx2)], z_mat[,c(col_idx1, col_idx2)]))
nat_mat <- nat_mat1 + nat_mat2
nat_mat <- avg_mat %*% nat_mat
case_gaussian_mean2 <- Matrix::colMeans(nat_mat[case_row_idx,,drop = F])
control_gaussian_mean2 <- Matrix::colMeans(nat_mat[control_row_idx,,drop = F])
test_vec2 <- case_gaussian_mean2 - control_gaussian_mean2
quantile(test_vec2)
length(which(test_vec2 < 0))/length(test_vec2)

teststat_vec <- pmax(pmin(test_vec2, 5), -5)
max_val <- max(abs(teststat_vec))
png("../../../../out/fig/Writeup12/Writeup12_sns_layer23_teststat2.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
main_vec <- c("Published DE", "Other interest\n(SFARI, DE other region)", "Housekeeping+Cell-cycle")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency", freq = T)
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()

###############################################

# attempt 3: just individuals
nat_mat1 <- exp(tcrossprod(x_mat, y_mat))
col_idx1 <- which(colnames(covariates) %in% c("diagnosis_ASD"))
col_idx2 <- grep("^individual", colnames(covariates))
nat_mat2 <- exp(tcrossprod(covariates[,c(col_idx1, col_idx2)], z_mat[,c(col_idx1, col_idx2)]))
nat_mat <- nat_mat1 + nat_mat2
nat_mat <- avg_mat %*% nat_mat
case_gaussian_mean <- Matrix::colMeans(nat_mat[case_row_idx,,drop = F])
control_gaussian_mean <- Matrix::colMeans(nat_mat[control_row_idx,,drop = F])
test_vec1 <- case_gaussian_mean - control_gaussian_mean
quantile(test_vec1)
length(which(test_vec1 < 0))/length(test_vec1)

teststat_vec <- pmax(pmin(test_vec1, 5), -5)
max_val <- max(abs(teststat_vec))
png("../../../../out/fig/Writeup12/Writeup12_sns_layer23_teststat3.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
main_vec <- c("Published DE", "Other interest\n(SFARI, DE other region)", "Housekeeping+Cell-cycle")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency", freq = T)
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()

###############################################

# attempt 4: all individual covariates aside from individual indicators
nat_mat1 <- exp(tcrossprod(x_mat, y_mat))
col_idx1 <- which(colnames(covariates) %in% c("age", "Intercept", "sex_F", "diagnosis_ASD", "post.mortem.hours", "region_ACC"))
nat_mat2 <- exp(tcrossprod(covariates[,c(col_idx1)], z_mat[,c(col_idx1)]))
nat_mat <- nat_mat1 + nat_mat2
nat_mat <- avg_mat %*% nat_mat
case_gaussian_mean <- Matrix::colMeans(nat_mat[case_row_idx,,drop = F])
control_gaussian_mean <- Matrix::colMeans(nat_mat[control_row_idx,,drop = F])
test_vec1 <- case_gaussian_mean - control_gaussian_mean
quantile(test_vec1)
length(which(test_vec1 < 0))/length(test_vec1)

teststat_vec <- pmax(pmin(test_vec1, 5), -5)
max_val <- max(abs(teststat_vec))
png("../../../../out/fig/Writeup12/Writeup12_sns_layer23_teststat4.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
main_vec <- c("Published DE", "Other interest\n(SFARI, DE other region)", "Housekeeping+Cell-cycle")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency", freq = T)
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()


###############################################

# attempt 5:just the data
library_size_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, which_fit = "param", what_obj = "init_library_size_variable")
library_size_variables <- library_size_variable
library_size_variables <- unique(c(library_size_variables, setdiff(colnames(covariates), c("Intercept", case_control_variable))))
library_size_variables <-  unique(c("Intercept", library_size_variables))

library_idx <- which(colnames(covariates) %in% library_size_variables)
library_mat <- exp(tcrossprod(
  covariates[,library_idx], z_mat[,library_idx]
))
library_mat <- pmax(library_mat, 1)
nat_mat <- eSVD_obj$dat / library_mat
nat_mat <- avg_mat %*% nat_mat
case_gaussian_mean5 <- Matrix::colMeans(nat_mat[case_row_idx,,drop = F])
control_gaussian_mean5 <- Matrix::colMeans(nat_mat[control_row_idx,,drop = F])
test_vec5 <- case_gaussian_mean5 - control_gaussian_mean5
quantile(test_vec5)
length(which(test_vec5 < 0))/length(test_vec5)

####################################

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj,
                             metadata = sns@meta.data,
                             covariate_individual = "individual")
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

selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]

############

idx <- which(names(eSVD_obj$teststat_vec) %in% selected_genes)
sfari_idx <- which(names(eSVD_obj$teststat_vec) %in% sfari_genes)
hk_idx <- which(names(eSVD_obj$teststat_vec) %in% hk_genes)

x_vec <- log2(case_gaussian_mean5) - log2(control_gaussian_mean5)
xlim <- quantile(x_vec, probs = c(0.01, 1))
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- logpvalue_vec
ylim <- range(y_vec)

# adjust x_vec
# median_idx <- which(y_vec <= stats::quantile(y_vec, probs = 0.05))
# x_vec <- x_vec - stats::median(x_vec[median_idx])

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

png("../../../../out/fig/Writeup12/Writeup12_sns_layer23_volcano_version5.png",
    height = 3500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in seq(0,7,by = .5)){
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

