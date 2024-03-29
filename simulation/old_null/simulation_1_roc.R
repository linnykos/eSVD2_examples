rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

#########################################

load("../../out/simulation/simulation_1.RData")
load("../../out/simulation/simulation_1_esvd.RData")
load("../../out/simulation/simulation_1_deseq2.RData")
load("../../out/simulation/simulation_1_sctransform.RData")
load("../../out/simulation/simulation_1_mast.RData")

true_de_idx <- which(true_fdr_vec < 0.05)
length(true_de_idx)

#############

compute_roc <- function(estimated_teststat_vec,
                        true_de_idx){
  ordering_est <- order(abs(estimated_teststat_vec), decreasing = T)
  n <- length(ordering_est)
  tpr <- sapply(0:n, function(i){
    if(i == 0) return(0)
    est_de_idx <- ordering_est[1:i]
    length(intersect(est_de_idx, true_de_idx))/length(true_de_idx)
  })
  fpr <- sapply(0:n, function(i){
    if(i == 0) return(0)
    if(i == n) return(1)
    est_de_idx <- ordering_est[1:i]
    true_nonde_idx <- setdiff(1:n, true_de_idx)
    length(intersect(est_de_idx, true_nonde_idx))/length(true_nonde_idx)
  })

  names(tpr) <- paste0("selected-", 0:n)
  names(fpr) <- names(tpr)

  list(tpr = tpr, fpr = fpr)
}


smooth_roc <- function(tpr, fpr){
  # both vectors are assumed to start at 0 and end at 1, where tpr rises faster than fpr
  mat <- cbind(fpr, tpr)

  # first smooth to be unimodal
  radian <- 45*2*pi/360
  rotation_mat <- matrix(c(cos(radian), sin(radian), -sin(radian), cos(radian)), 2, 2)
  mat2 <- mat %*% rotation_mat
  mat2_old <- mat2
  ordering <- order(mat2[,1], decreasing = F)
  res <- UniIsoRegression::reg_1d(y_vec = mat2[ordering,2],
                                  w_vec = rep(1, nrow(mat2)),
                                  metric = 2,
                                  unimodal = T)
  mat2[ordering,2] <- res
  # plot(mat2[,1], mat2_old[,2]); points(mat2[,1], mat2[,2], col = 2)
  mat_new <- mat2 %*% t(rotation_mat)
  tpr <- mat_new[,2]; fpr <- mat_new[,1]

  # next, smooth to be monotonic
  res <- stats::isoreg(x = fpr, y = tpr)
  tpr_new <- res$yf

  # last, make sure the curve doesn't dip below the diagonal
  fpr <- pmax(pmin(fpr, 1), 0)
  tpr_new <- pmax(pmin(tpr_new, 1), 0)
  n <- length(tpr_new)
  for(i in 1:n){
    tpr_new[i] <- max(tpr_new[i], fpr[i])
  }

  list(tpr = tpr_new,
       fpr = fpr)
}

roc_fdr_point <- function(pvalue_vec,
                          true_de_idx,
                          fdr_threshold = 0.05){
  n <- length(pvalue_vec)
  fdr_vec <- stats::p.adjust(pvalue_vec, method = "BH")
  fdr_idx <- which(fdr_vec <= fdr_threshold)

  true_nonde_idx <- setdiff(1:n, true_de_idx)
  tpr <- length(intersect(true_de_idx,fdr_idx))/length(true_de_idx)
  fpr <- length(intersect(true_nonde_idx,fdr_idx))/length(true_nonde_idx)

  c(tpr = tpr, fpr = fpr,
    len = length(fdr_idx),
    intersection = length(intersect(true_de_idx,fdr_idx)))
}

search_for_point_on_curve <- function(fpr,
                                      fpr_curve,
                                      tpr,
                                      tpr_curve){
  n <- length(fpr_curve)
  dist_vec <- (fpr-fpr_curve)^2 + (tpr-tpr_curve)^2
  idx <- which.min(dist_vec)
  list(tpr = tpr_curve[idx],
       fpr = fpr_curve[idx])
}

##########################

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

esvd_pvalue <- logpvalue_vec[paste0("gene-", 1:p)]
esvd_roc <- compute_roc(estimated_teststat_vec = esvd_pvalue,
                        true_de_idx = true_de_idx)
esvd_roc <- smooth_roc(tpr = esvd_roc$tpr,
                       fpr = esvd_roc$fpr)
esvd_point <- roc_fdr_point(pvalue_vec = 10^(-esvd_pvalue),
                            true_de_idx = true_de_idx,
                            fdr_threshold = 0.1)
esvd_point <- search_for_point_on_curve(fpr = esvd_point["fpr"],
                                        fpr_curve = esvd_roc$fpr,
                                        tpr = esvd_point["tpr"],
                                        tpr_curve = esvd_roc$tpr)
esvd_point

deseq_pvalue <- deseq2_res[paste0("gene-", 1:p), "pvalue"]
deseq2_roc <- compute_roc(estimated_teststat_vec = -log10(deseq_pvalue),
                          true_de_idx = true_de_idx)
deseq2_roc <- smooth_roc(tpr = deseq2_roc$tpr,
                         fpr = deseq2_roc$fpr)
deseq2_point <- roc_fdr_point(pvalue_vec = deseq_pvalue,
                              true_de_idx = true_de_idx,
                              fdr_threshold = 0.1)
deseq2_point <- search_for_point_on_curve(fpr = deseq2_point["fpr"],
                                          fpr_curve = deseq2_roc$fpr,
                                          tpr = deseq2_point["tpr"],
                                          tpr_curve = deseq2_roc$tpr)
deseq2_point

sctransform_pvalue <- de_result[paste0("gene-", 1:p), "p_val"]
sctransform_roc <- compute_roc(estimated_teststat_vec = -log10(sctransform_pvalue),
                               true_de_idx = true_de_idx)
sctransform_roc <- smooth_roc(tpr = sctransform_roc$tpr,
                              fpr = sctransform_roc$fpr)
sctransform_point <- roc_fdr_point(pvalue_vec = sctransform_pvalue,
                                   true_de_idx = true_de_idx,
                                   fdr_threshold = 0.1)
sctransform_point <- search_for_point_on_curve(fpr = sctransform_point["fpr"],
                                               fpr_curve = sctransform_roc$fpr,
                                               tpr = sctransform_point["tpr"],
                                               tpr_curve = sctransform_roc$tpr)

mast_pvalue <- mast_pval_glmer[paste0("gene-", 1:p)]
mast_roc <- compute_roc(estimated_teststat_vec = -log10(mast_pval_glmer),
                        true_de_idx = true_de_idx)
mast_roc <- smooth_roc(tpr = mast_roc$tpr,
                       fpr = mast_roc$fpr)
mast_point <- roc_fdr_point(pvalue_vec = mast_pval_glmer,
                            true_de_idx = true_de_idx,
                            fdr_threshold = 0.1)
mast_point <- search_for_point_on_curve(fpr = mast_point["fpr"],
                                        fpr_curve = mast_roc$fpr,
                                        tpr = mast_point["tpr"],
                                        tpr_curve = mast_roc$tpr)

##################

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
yellow_col <- rgb(255, 205, 114, maxColorValue = 255)
blue_col <- rgb(48, 174, 255, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)

png("../../out/fig/simulation/simulation_1_roc.png",
    height = 2000, width = 2000,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
     xaxt = "n", yaxt = "n", bty = "n",
     cex.lab = 1.25, type = "n",
     xlab = "", ylab = "")
lines(sctransform_roc$fpr, sctransform_roc$tpr, col = blue_col, lwd = 4)
lines(mast_roc$fpr, mast_roc$tpr, col = purple_col, lwd = 4)
lines(deseq2_roc$fpr, deseq2_roc$tpr, col = yellow_col, lwd = 4)
lines(esvd_roc$fpr, esvd_roc$tpr, col = orange_col, lwd = 4)

lines(c(0,1), c(0,1), col = 1, lwd = 2, lty = 2)

points(sctransform_point$fpr, sctransform_point$tpr, col = "black", pch = 16, cex = 3)
points(sctransform_point$fpr, sctransform_point$tpr, col = "white", pch = 16, cex = 2.5)
points(sctransform_point$fpr, sctransform_point$tpr, col = blue_col, pch = 16, cex = 2)

points(mast_point$fpr, mast_point$tpr, col = "black", pch = 16, cex = 3)
points(mast_point$fpr, mast_point$tpr, col = "white", pch = 16, cex = 2.5)
points(mast_point$fpr, mast_point$tpr, col = purple_col, pch = 16, cex = 2)

points(deseq2_point$fpr, deseq2_point$tpr, col = "black", pch = 16, cex = 3)
points(deseq2_point$fpr, deseq2_point$tpr, col = "white", pch = 16, cex = 2.5)
points(deseq2_point$fpr, deseq2_point$tpr, col = yellow_col, pch = 16, cex = 2)

points(esvd_point$fpr, esvd_point$tpr, col = "black", pch = 16, cex = 3)
points(esvd_point$fpr, esvd_point$tpr, col = "white", pch = 16, cex = 2.5)
points(esvd_point$fpr, esvd_point$tpr, col = orange_col, pch = 16, cex = 2)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
