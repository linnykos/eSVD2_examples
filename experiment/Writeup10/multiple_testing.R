multttest_calibrate <- function(teststat_vec,
                                null_dens,
                                null_x,
                                fdr_cutoff = 0.05,
                                two_sided = T){
  normalizing_val <- sum(null_dens)
  cumsum_vec <- cumsum(null_dens)
  null_median <- null_x[which.min(abs(cumsum_vec - normalizing_val/2))]
  p_val_vec <- sapply(teststat_vec, function(teststat){
    if(teststat < null_median){
      tail_prob <- cumsum_vec[which.min(abs(null_x - teststat))]
    } else {
      tail_prob <- normalizing_val - cumsum_vec[which.min(abs(null_x - teststat))]
    }
    tail_prob <- tail_prob/normalizing_val

    if(two_sided) 2*tail_prob else tail_prob
  })

  fdr_vec <- stats::p.adjust(p_val_vec, method = "BH")
  fdr_idx <- which(fdr_vec <= fdr_cutoff)
  neglogp_val_vec <- -log10(p_val_vec)
  if(any(is.infinite(neglogp_val_vec))){
    inf_idx <- which(is.infinite(neglogp_val_vec))
    max_val <- max(neglogp_val_vec[-inf_idx])
    neglogp_val_vec[inf_idx] <- 1.5*max_val
  }

  list(p_val = p_val_vec,
       neglog_p_val = neglogp_val_vec,
       fdr = fdr_vec,
       idx = names(teststat_vec)[fdr_idx])
}

multtest_lfdr <- function(){

}

multtest_predictionrecursion <- function(){

}

write_genes <- function(vec,
                        file){
  fileConn <- file(file)
  writeLines(vec, fileConn)
  close(fileConn)
}

#############################################
