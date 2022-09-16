rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

#########################################

load("../eSVD2_examples/simulation/simulation_2.RData")

mean_mat <- exp(nat_mat)
case_idx <- which(covariates[,"cc"] == 1)
control_idx <- which(covariates[,"cc"] == 0)
case_mean <- colMeans(mean_mat[case_idx,])
control_mean <- colMeans(mean_mat[control_idx,])
diff_mean <- case_mean - control_mean
var_mat <- sweep(x = mean_mat, MARGIN = 2, STATS = nuisance_vec, FUN = "/")

res <- eSVD2:::compute_test_statistic.default(
  input_obj = mean_mat,
  posterior_var_mat = var_mat,
  case_individuals = case_individuals,
  control_individuals = control_individuals,
  individual_vec = individual_vec
)
true_teststat_vec_safe <- res$teststat_vec

load("../eSVD2_examples/simulation/simulation_2_esvd.RData")
load("../eSVD2_examples/simulation/simulation_2_deseq2.RData")
load("../eSVD2_examples/simulation/simulation_2_sctransform.RData")

true_teststat_vec <- true_teststat_vec_safe

#############

compute_roc <- function(estimated_teststat_vec,
                        true_teststat_vec){
  stopifnot(length(estimated_teststat_vec) == length(true_teststat_vec),
            length(names(true_teststat_vec)) > 0,
            all(names(estimated_teststat_vec) == names(true_teststat_vec)))

  len <- length(true_teststat_vec)
  ordering_est <- order(abs(estimated_teststat_vec), decreasing = T)
  ordering_true <- order(abs(true_teststat_vec), decreasing = T)
  tpr <- sapply(0:len, function(i){
    if(i == 0) return(0)
    est_de_idx <- ordering_est[1:i]
    true_de_idx <- ordering_true[1:i]
    length(intersect(est_de_idx, true_de_idx))/length(true_de_idx)
  })
  fpr <- sapply(0:len, function(i){
    if(i == 0) return(0)
    if(i == len) return(1)
    est_de_idx <- ordering_est[1:i]
    true_nonde_idx <- ordering_true[(i+1):len]
    length(intersect(est_de_idx, true_nonde_idx))/length(true_nonde_idx)
  })

  names(tpr) <- paste0("selected-", 0:len)
  names(fpr) <- names(tpr)

  # plot(fpr,tpr,asp=T)
  list(tpr = tpr,
       fpr = fpr)
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

  list(tpr = res$yf,
       fpr = fpr)
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

esvd_roc <- compute_roc(estimated_teststat_vec = logpvalue_vec,
                        true_teststat_vec = true_teststat_vec)
esvd_roc <- smooth_roc(tpr = esvd_roc$tpr,
                       fpr = esvd_roc$fpr)

deseq2_roc <- compute_roc(estimated_teststat_vec = -log10(deseq2_res$pvalue),
                          true_teststat_vec = true_teststat_vec)
deseq2_roc <- smooth_roc(tpr = deseq2_roc$tpr,
                         fpr = deseq2_roc$fpr)

sctransform_roc <- compute_roc(estimated_teststat_vec = -log10(de_result$p_val),
                               true_teststat_vec = true_teststat_vec)
sctransform_roc <- smooth_roc(tpr = sctransform_roc$tpr,
                              fpr = sctransform_roc$fpr)

##################

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
yellow_col <- rgb(255, 205, 114, maxColorValue = 255)
blue_col <- rgb(48, 174, 255, maxColorValue = 255)

png("../../out/fig/simulation/simulation_2_roc.png",
    height = 2000, width = 2000,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
     xaxt = "n", yaxt = "n", bty = "n",
     cex.lab = 1.25, type = "n",
     xlab = "", ylab = "")
points(deseq2_roc$fpr, deseq2_roc$tpr, col = yellow_col, pch = 16)
points(sctransform_roc$fpr, sctransform_roc$tpr, col = blue_col, pch = 16)
points(esvd_roc$fpr, esvd_roc$tpr, col = orange_col, pch = 16)
lines(c(0,1), c(0,1), col = 1, lwd = 2, lty = 2)
axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
