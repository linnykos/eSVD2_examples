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
