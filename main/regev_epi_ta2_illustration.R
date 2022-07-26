rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta2-inflamed_esvd.RData")

####################

set.seed(10)
umap_res <- Seurat::RunUMAP(eSVD_obj$fit_Second$x_mat)
umap_mat <- umap_res@cell.embeddings

case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(3)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(3)

tab <- table(regevEpi$Subject, regevEpi$Subject_Disease)
n <- ncol(regevEpi)
color_vec <- rep(NA, n)
case_subj <- rownames(tab)[which(tab[,"Colitis"] != 0)]
case_vec <- sample(1:3, replace = T, size = length(case_subj))
control_subj <- rownames(tab)[which(tab[,"Colitis"] == 0)]
control_vec <- sample(1:3, replace = T, size = length(case_subj))

for(i in 1:length(case_vec)){
  color_vec[regevEpi$Subject == case_subj[i]] <- case_color_palette[case_vec[i]]
}
for(i in 1:length(case_vec)){
  color_vec[regevEpi$Subject == control_subj[i]] <- control_color_palette[control_vec[i]]
}

png("../../../out/fig/main/regevEpi_ta2_illustration_umap.png",
    height = 1500, width = 2000,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
plot(umap_mat[,1], umap_mat[,2], col = color_vec,
     xaxt = "n", yaxt = "n", bty = "n", pch = 16)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

#######################3

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj,
                     metadata = regevEpi@meta.data,
                     covariate_individual = "Subject")
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
pvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    stats::pnorm(x, mean = null_mean, sd = null_sd)*2
  } else {
    (1-stats::pnorm(x, mean = null_mean, sd = null_sd))*2
  }
})
logpvalue_vec <- -log10(pvalue_vec)
idx <- which(fdr_vec < 0.1)

max(logpvalue_vec[-idx])
min(logpvalue_vec[idx])
gaussian_teststat[which(logpvalue_vec == max(logpvalue_vec[-idx]))]
fdr_vec[which(logpvalue_vec == max(logpvalue_vec[-idx]))]

tab <- table(regevEpi$Subject, regevEpi$Subject_Disease)
indiv_cases <- rownames(tab)[which(tab[,"Colitis"] != 0)]
indiv_controls <- rownames(tab)[which(tab[,"Colitis"] == 0)]
indiv_vec <- factor(as.character(regevEpi$Subject))

p <- length(gaussian_teststat)
posterior_mat <- eSVD_obj$fit_Second$posterior_mean_mat
x_vec <- sapply(1:p, function(j){
  log2(mean(posterior_mat[which(indiv_vec %in% indiv_cases),j])) - log2(mean(posterior_mat[which(indiv_vec %in% indiv_controls),j]))
})

png("../../../out/fig/main/regevEpi_ta2_illustration_volcano.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = logpvalue_vec, col = "gray",
     xaxt = "n", yaxt = "n", bty = "n", pch = 16,
     cex.lab = 1.25)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
y_val <- (max(logpvalue_vec[-idx]) + min(logpvalue_vec[idx]))/2
lines(x = c(-100,100), y = rep(y_val, 2), lwd = 2, lty = 2, col = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

#######################################

input_obj <- eSVD_obj
metadata <- regevEpi@meta.data
covariate_individual <- "Subject"

case_control_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
covariates <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
cc_vec <- covariates[,case_control_variable]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])

latest_Fit <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
posterior_mean_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_mean_mat", which_fit = latest_Fit)
posterior_var_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_var_mat", which_fit = latest_Fit)

individual_vec <- metadata[,covariate_individual]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])
tmp <- eSVD2:::.determine_individual_indices(case_individuals = case_individuals,
                                             control_individuals = control_individuals,
                                             covariate_individual = covariate_individual,
                                             metadata = metadata)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- eSVD2:::.construct_averaging_matrix(idx_list = all_indiv_idx,
                                               n = nrow(posterior_mean_mat))
avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

library_size_variables <- "Log_UMI"
library_size_variables <- unique(c(library_size_variables, setdiff(colnames(eSVD_obj$covariates), c("Intercept", case_control_variable))))
library_size_variables <- unique(c("Intercept", library_size_variables))

library_idx <- which(colnames(eSVD_obj$covariates) %in% library_size_variables)
library_mat <- exp(tcrossprod(
  eSVD_obj$covariates[,library_idx], eSVD_obj$fit_Second$z_mat[,library_idx]
))

nat_mat1 <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates[,-library_idx], eSVD_obj$fit_Second$z_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)

nuisance_vec <- eSVD_obj$fit_Second$nuisance_vec

eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      alpha_max = 10000,
                                      bool_adjust_covariates = F,
                                      bool_covariates_as_library = T)
p <- length(nuisance_vec)
tmp1 <- sapply(1:p, function(j){
  stats::cor(mean_mat_nolib[,j], eSVD_obj$fit_Second$posterior_mean_mat[,j])
})
dat <- as.matrix(eSVD_obj$dat)/library_mat
tmp2 <- sapply(1:p, function(j){
  abs(diff(range(mean_mat_nolib[,j])) - diff(range(dat[,j])))
})

# gene <- names(logpvalue_vec)[which.min(abs(logpvalue_vec - 4))]
# gene <- names(nuisance_vec)[which.min(abs(nuisance_vec - 100))]
set.seed(0); idx <- which(tmp2 <= 0.5); gene <- names(nuisance_vec)[idx[which.min(tmp1[idx])]]
# avg_posterior_mean_mat[,gene]
# sqrt(avg_posterior_var_mat[,gene])

library_vec <- library_mat[,gene]
obs_vec <- as.numeric(eSVD_obj$dat[,gene])/library_vec
fit_vec <- mean_mat_nolib[,gene]
posterior_vec <- eSVD_obj$fit_Second$posterior_mean_mat[,gene]
posterior_vec2 <- (as.numeric(eSVD_obj$dat[,gene]) + mean_mat_nolib[,gene]*nuisance_vec[gene])/(library_vec + nuisance_vec[gene])
n <- length(y_vec)

png("../../../out/fig/main/regevEpi_ta2_illustration_shrinkage.png",
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,3,0.3,0.3))
set.seed(10)
idx <- sample(n, size = 200)
plot(x = obs_vec[idx], y = fit_vec[idx], col = "gray",
     xaxt = "n", yaxt = "n", bty = "n",
     lwd = 2, cex = 2, asp = T)
lines(x = c(-100,100), y = c(-100,100), col = 2, lwd = 2, lty = 2)
points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, col = "white", cex = 3)
points(x = posterior_vec[idx], y = fit_vec[idx], pch = 16, cex = 2)
axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
axis(side = 2, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()

split_idx <- length(case_individuals); k <- nrow(avg_posterior_var_mat)
selected_case_subj <- order(avg_posterior_var_mat[1:split_idx,gene], decreasing = T)[1:3]
selected_control_subj <- order(avg_posterior_var_mat[(split_idx+1):k,gene], decreasing = T)[1:3] + split_idx
subj_idx <- c(selected_case_subj,selected_control_subj)

xlim <- c(0.25, 2.25)
avg_posterior_mean_vec <- avg_posterior_mean_mat[subj_idx,gene]
avg_posterior_mean_vec[4:6] <- avg_posterior_mean_vec[4:6]+.25
avg_posterior_mean_vec[4:6] <- avg_posterior_mean_vec[6:4]
avg_posterior_var_vec <- avg_posterior_var_mat[subj_idx,gene]
avg_posterior_var_vec[4:6] <- avg_posterior_var_vec[6:4]
gaussian_list <- lapply(1:6, function(i){
  mean_val <- avg_posterior_mean_vec[i]
  var_val <- avg_posterior_var_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sqrt(var_val))
  yseq <- yseq - min(yseq)
  yseq <- yseq/4*1.5
  cbind(xseq, yseq)
})

col_vec <- c(case_color_palette, control_color_palette)
png(paste0("../../../out/fig/main/regevEpi_ta2_illustration_gaussians.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,0,0,0), bg = NA)
plot(NA,
     xlim = c(0.25, 2.25),
     ylim = c(0, 7.5),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
for(i in 6:1){
  graphics::polygon(x = c(gaussian_list[[i]][,1], rev(gaussian_list[[i]][,1])),
                    y = 0.5+(i-1)+c(gaussian_list[[i]][,2], rep(0, nrow(gaussian_list[[i]]))),
                    col = col_vec[i])
}
axis(side = 1, at = c(0.5, 1, 1.5, 2), lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()

######################

mean_vec <- c(mean(avg_posterior_mean_vec[1:3]), mean(avg_posterior_mean_vec[4:6]))
var_vec <- c(mean(avg_posterior_var_vec[1:3]) + mean(avg_posterior_mean_vec[1:3]^2) - mean(avg_posterior_mean_vec[1:3])^2,
             mean(avg_posterior_var_vec[4:6]) + mean(avg_posterior_mean_vec[4:6]^2) - mean(avg_posterior_mean_vec[4:6])^2)

gaussian_list <- lapply(1:2, function(i){
  mean_val <- mean_vec[i]
  var_val <- var_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sqrt(var_val))
  yseq <- yseq - min(yseq)
  cbind(xseq, yseq)
})

ymax <- max(sapply(gaussian_list, function(x){x[,"yseq"]}))
png(paste0("../../../out/fig/main/regevEpi_ta2_illustration_gaussians-avg.png"),
    height = 1500, width = 1500,
    units = "px", res = 500)
par(mar = c(3,0,0,0), bg = NA)
plot(NA,
     xlim = c(0.25, 2.25),
     ylim = c(0, ymax+.5),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
col_vec <- c(rgb(244, 84, 84, 0.5*255, maxColorValue = 255),
             rgb(27, 198, 245, 0.5*255, maxColorValue = 255))
for(i in 1:2){
  graphics::polygon(x = c(gaussian_list[[i]][,1], rev(gaussian_list[[i]][,1])),
                    y = 0.5+c(gaussian_list[[i]][,2], rep(0, nrow(gaussian_list[[i]]))),
                    col = col_vec[i])
}
axis(side = 1, at = c(0.5, 1, 1.5, 2), lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()
