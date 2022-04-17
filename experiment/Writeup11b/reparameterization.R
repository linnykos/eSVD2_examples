final_reparameterization <- function(esvd_res,
                                     exclude_vars,
                                     verbose = T){
  k <- ncol(esvd_res$x_mat)
  var_names <- colnames(esvd_res$covariates)
  var_names <- setdiff(var_names, exclude_vars)
  covariate_mat <- esvd_res$covariates[,var_names,drop = F]

  for(j in 1:k){
    if(verbose) print(paste0("Finished latent variable ", j))
    df_tmp <- as.data.frame(cbind(esvd_res$x_mat[,j], covariate_mat))
    colnames(df_tmp) <- c("y", paste0("x", 1:ncol(covariate_mat)))
    lm_res <- stats::lm(y ~ . - 1, data = df_tmp)
    coef_vec <- stats::coef(lm_res)
    residual_vec <- esvd_res$x_mat[,j] - stats::predict(lm_res)

    esvd_res$x_mat[,j] <- residual_vec
    esvd_res$b_mat[,var_names] <- esvd_res$b_mat[,var_names] +
      tcrossprod(esvd_res$y_mat[,j], coef_vec)
  }

  tmp <- eSVD2:::.reparameterize(esvd_res$x_mat, esvd_res$y_mat, equal_covariance = T)
  esvd_res$x_mat <- tmp$x_mat
  esvd_res$y_mat <- tmp$y_mat

  esvd_res
}
