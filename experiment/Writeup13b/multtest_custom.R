multtest_custom <- function(eSVD_obj){
  gaussian_teststat <- eSVD_obj$pvalue_list$gaussian_teststat
  # median_val <- stats::median(gaussian_teststat)
  # mad_val <- stats::median(abs(gaussian_teststat - median_val))
  #
  # z_vec <- (gaussian_teststat - median_val)/mad_val
  # names(z_vec) <- names(gaussian_teststat)
  #
  # null_mean <- stats::median(z_vec)
  # null_sd <- stats::median(abs(z_vec - null_mean))

  locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
  null_mean <- locfdr_res$fp0["mlest", "delta"]
  null_sd <- locfdr_res$fp0["mlest", "sigma"]

  log10pvalue_vec <- sapply(gaussian_teststat, function(x){
    if(x < null_mean) {
      Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
    } else {
      Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
    }
  })
  log10pvalue_vec <- -(log10pvalue_vec/log(10) + log10(2))
  names(log10pvalue_vec) <- names(gaussian_teststat)

  pvalue_vec <- sapply(gaussian_teststat, function(x){
    if(x < null_mean) {
      2*Rmpfr::pnorm(x, mean = null_mean, sd = null_sd)
    } else {
      2*Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd)
    }
  })
  names(pvalue_vec) <- names(gaussian_teststat)

  fdr_vec <- stats::p.adjust(pvalue_vec, method = "BH")
  names(fdr_vec) <- names(gaussian_teststat)

  pvalue_list <- list(
    df_vec = eSVD_obj$pvalue_list$df_vec,
    fdr_vec = fdr_vec,
    gaussian_teststat = gaussian_teststat,
    log10pvalue = log10pvalue_vec,
    null_mean = null_mean,
    null_sd = null_sd,
    pvalue_vec = pvalue_vec
  )

  eSVD_obj[["pvalue_list"]] <- pvalue_list
  eSVD_obj
}
