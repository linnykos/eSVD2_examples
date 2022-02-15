initialize_esvd2 <- function(dat,
                             k,
                             family = "poisson",
                             covariates = NULL,
                             offset_vec = rep(0, nrow(dat)),
                             column_set_to_one = NULL,
                             tol = 1e-3,
                             verbose = 0){
  stopifnot(is.character(family),
            family == "poisson",
            all(is.null(covariates)) || is.matrix(covariates),
            length(offset_vec) == nrow(dat),
            k <= ncol(dat), k > 0, k %% 1 == 0)

  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0

  ##########

  # step 1: fix the covariates
  ## [[note to self: remove this hard-code]]
  stopifnot(which(colnames(covariates) == "Intercept") == 1)
  indiv_idx <- grep("individual_", colnames(covariates))
  keep_idx <- setdiff(c(2:ncol(covariates)), indiv_idx)
  covariates[,keep_idx] <- scale(covariates[,keep_idx], center = F, scale = T)

  # regress indiv_idx to all other covariates
  for(j in indiv_idx){
    df_tmp <- data.frame(covariates[,j], covariates[,c(1,keep_idx)])
    colnames(df_tmp)[1] <- "tmp"
    lm_fit <- stats::lm(tmp ~ . - 1, data = df_tmp)
    vec_tmp <- stats::residuals(lm_fit)

    if(eSVD2:::.l2norm(vec_tmp) > tol) {
      covariates[,j] <- vec_tmp
    } else {
      covariates[,j] <- rep(NA, n)
    }
  }

  na_idx <- which(is.na(covariates[1,]))
  if(length(na_idx) > 0) covariates <- covariates[,-na_idx,drop=F]

  ##########

  # step 2: poisson GLMS
  # covariates_noint <- covariates[,-which(colnames(covariates) == "Intercept")]
  coef_mat <- t(sapply(1:p, function(j){
    if(verbose >= 1 && p > 10 && j %% floor(p/10) == 0) cat('*')

    df <- as.data.frame(cbind(y = mat[,j], covariates))
    colnames(df)[1] <- "tmp"
    glm_fit <- stats::glm(tmp ~ . - 1, data = df)
    ## [[note to self: can be improved using the stat helper functions]]
    stats::coef(glm_fit)
  }))
  # coef_mat <- cbind(0, coef_mat)
  colnames(coef_mat) <- colnames(covariates)
  rownames(coef_mat) <- colnames(dat)

  na_idx <- which(is.na(coef_mat[1,]))
  if(length(na_idx) > 0) {
    covariates <- covariates[,-na_idx,drop=F]
    coef_mat <- coef_mat[,-na_idx,drop=F]
  }

  #########

  # step 3: SVD
  dat_transform <- log1p(dat)
  nat_mat <- tcrossprod(covariates, coef_mat)
  residual_mat <- dat_transform - nat_mat

  svd_res <- eSVD2:::.svd_truncated(residual_mat,
                                    K = k,
                                    symmetric = F,
                                    rescale = F,
                                    mean_vec = NULL,
                                    sd_vec = NULL,
                                    K_full_rank = F)
  x_init <- eSVD2:::.mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_init <- eSVD2:::.mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  rownames(x_init) <- rownames(dat)
  rownames(y_init) <- colnames(dat)

  structure(list(x_mat = x_init, y_mat = y_init, b_mat = coef_mat,
                 covariates = covariates,
                 nuisance_param_vec = rep(1, ncol(dat)),
                 offset_vec = offset_vec),
            class = "eSVD")
}
