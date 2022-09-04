rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)

load("../../../../out/Writeup12/adams_T_esvd.RData")
date_of_run

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

################################

# the hypothetical: what if we include the individual covariates inside the calculation of the "signal"?

eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj$teststat_vec <- NULL
eSVD_obj$case_mean <- NULL
eSVD_obj$control_mean <- NULL

input_obj = eSVD_obj
alpha_max = NULL
bool_adjust_covariates = F
bool_covariates_as_library = T
bool_return_components = F
bool_stabilize_underdispersion = F
library_min = 1
nuisance_lower_quantile = 0.01
pseudocount = 1

dat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "dat", which_fit = NULL)
covariates <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
latest_Fit <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
esvd_res <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = NULL, which_fit = latest_Fit)
nuisance_vec <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "nuisance", which_fit = latest_Fit)
case_control_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, which_fit = "param", what_obj = "init_case_control_variable")
library_size_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, which_fit = "param", what_obj = "init_library_size_variable")
bool_library_includes_interept <- eSVD2:::.get_object(eSVD_obj = input_obj, which_fit = "param", what_obj = "nuisance_bool_library_includes_interept")

input_obj = as.matrix(eSVD_obj$dat)
case_control_idx <- which(colnames(covariates) == case_control_variable)

library_size_variables <- library_size_variable
if(bool_covariates_as_library) library_size_variables <- unique(c(library_size_variables, setdiff(colnames(covariates),
                                                                                                  c("Intercept", case_control_variable,
                                                                                                  colnames(covariates)[grep("Subject_Identity", colnames(covariates))]))))
if(bool_library_includes_interept) library_size_variables <-  unique(c("Intercept", library_size_variables))

library_idx <- which(colnames(covariates) %in% library_size_variables)

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(covariates[,-library_idx], esvd_res$z_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  covariates[,library_idx], esvd_res$z_mat[,library_idx]
))
if(!is.null(library_min)) library_mat <- pmax(library_mat, library_min)

nuisance_vec <- pmax(nuisance_vec,
                     stats::quantile(nuisance_vec, probs = nuisance_lower_quantile))
Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
AplusAlpha <- input_obj + Alpha

if(!is.null(alpha_max)) AplusAlpha <- pmin(AplusAlpha, alpha_max)

SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2

rownames(posterior_mean_mat) <- rownames(input_obj)
rownames(posterior_var_mat) <- rownames(input_obj)
colnames(posterior_mean_mat) <- colnames(input_obj)
colnames(posterior_var_mat) <- colnames(input_obj)

eSVD_obj$fit_Second$posterior_mean_mat <- posterior_mean_mat
eSVD_obj$fit_Second$posterior_var_mat <- posterior_var_mat

eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           verbose = 1)

###########################

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("B", "Macrophage", "Macrophage Alveolar", "NK"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("Macrophages_disease_vs_control_.csv", "Monocytes_disease_vs_control_.csv",
              "B_Cells_disease_vs_control_.csv", "NK_Cells_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
habermann_df_genes <- setdiff(habermann_df_genes, adams_df_genes)
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
de_genes_others <- unique(c(adams_df_genes_others, habermann_df_genes_others))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
de_genes_others <- setdiff(de_genes_others, de_genes)
cycling_genes <- setdiff(cycling_genes, c(de_genes_others, de_genes))
hk_genes <- setdiff(hk_genes, c(cycling_genes, de_genes_others, de_genes))

gene_names <- names(eSVD_obj$teststat_vec)
adam_idx <- which(gene_names %in% adams_df_genes)
habermann_idx <- which(gene_names %in% habermann_df_genes)
de_other_idx <- which(gene_names %in% de_genes_others)
cycling_idx <- which(gene_names %in% cycling_genes)
hk_idx <- which(gene_names %in% hk_genes)

##########################################

# make volcano plot

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
idx <- order(logpvalue_vec, decreasing = T)[1:length(adams_df_genes)]
length(idx)
length(intersect(idx, c(adam_idx)))
length(intersect(idx, c(adam_idx, habermann_idx)))

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(adam_idx)
n <- length(gaussian_teststat) - m
k <- length(idx)
x <- length(intersect(idx, c(adam_idx)))
fisher <- stats::dhyper(x = x, m = m, n = n, k = k, log = F)
fisher

x_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
xlim <- range(x_vec)
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- abs(eSVD_obj$teststat_vec)
ylim <- range(y_vec)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

png("../../../../out/fig/Writeup12/adams_T_volcano_esvd_postprocess2.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in seq(0,3,by = 0.25)){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}
points(x = x_vec, y = y_vec,
       col = rgb(0.6, 0.6, 0.6, 0.3), pch = 16)
points(x = x_vec[idx], y = y_vec[idx],
       col = orange_col, pch = 16, cex = 1.5)
points(x = x_vec[setdiff(adam_idx, idx)], y = y_vec[setdiff(adam_idx, idx)],
       col = purple_col, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = green_col_trans, pch = 16, cex = 1)
points(x = x_vec[intersect(idx, adam_idx)], y = y_vec[intersect(idx, adam_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx, adam_idx)], y = y_vec[intersect(idx, adam_idx)],
       col = purple_col, pch = 1, cex = 2, lwd = 2)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

png("../../../../out/fig/Writeup12/adams_T_volcano_esvd_postprocess2_stats.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
plot(x = 1:10, y = 1:10, type = "n",
     main = paste0("Fisher = ", fisher, "\nTotal: ",
                   length(adam_idx), ", inter: ", length(intersect(idx, adam_idx))))
graphics.off()


