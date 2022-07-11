reparameterization2 <- function(eSVD_obj,
                                fit_name,
                                omitted_variables){

  x_mat <- eSVD_obj[[fit_name]]$x_mat
  y_mat <- eSVD_obj[[fit_name]]$y_mat
  z_mat <- eSVD_obj[[fit_name]]$z_mat
  k <- ncol(x_mat)
  p <- nrow(y_mat)

  covariate_mat <- eSVD_obj$covariates
  covariate_mat <- covariate_mat[,which(!colnames(covariate_mat) %in% omitted_variables),drop=F]
  covariate_mat2 <- covariate_mat[,which(colnames(covariate_mat) != "Intercept")]

  for(ell in 1:k){
    tmp_df <- cbind(x_mat[,ell], covariate_mat2)
    colnames(tmp_df)[1] <- "x"
    tmp_df <- as.data.frame(tmp_df)

    lm_res <- stats::lm(x ~ . , data = tmp_df)
    coef_vec <- stats::coef(lm_res)
    names(coef_vec)[1] <- "Intercept"
    print(paste0(ell, ": R2 of ", round(summary(lm_res)$r.squared, 2)))

    for(j in 1:p){
      z_mat[j,names(coef_vec)] <- z_mat[j,names(coef_vec)] + coef_vec*y_mat[j,ell]
    }

    x_mat[,ell] <- stats::residuals(lm_res)
  }

  res <- eSVD2:::.reparameterize(x_mat, y_mat, equal_covariance = T)
  x_mat <- res$x_mat; y_mat <- res$y_mat

  eSVD_obj[[fit_name]]$x_mat <- x_mat
  eSVD_obj[[fit_name]]$y_mat <- y_mat
  eSVD_obj[[fit_name]]$z_mat <- z_mat
}
