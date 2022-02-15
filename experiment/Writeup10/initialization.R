initialize_esvd2 <- function(dat,
                             k,
                             family,
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

  family <- .string_to_distr_funcs(family)
  if(family$name != "gaussian") stopifnot(all(dat[!is.na(dat)] >= 0))

  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0

  ##########

  # step 1: fix the covariates

  for(j in 1:ncol(covariates)){
    if(colnames(covariates)[j] == "Intercept") next()
    covariates[,j] <- scale(covariates[,j], center = T, scale = T)
  }

  regress_idx <- which(colnames(covariates) == "diagnosis_ASD")
  indiv_idx <- grep("individual_", colnames(covariates))
  keep_idx <- setdiff(c(1:ncol(covariates)), c(regress_idx, indiv_idx))

  # regress diagnosis to keep_idx
  df_tmp <- data.frame(covariates[,regress_idx], covariates[,keep_idx])
  colnames(df_tmp)[1] <- "tmp"
  lm_fit <- stats::lm("tmp ~ . - 1", data = df_tmp)
  vec_tmp <- stats::residuals(lm_fit)
  covariates[,regress_idx] <- vec_tmp

  # regress indiv_idx to all other covariates
  for(j in indiv_idx){
    df_tmp <- data.frame(covariates[,j], covariates[,c(regress_idx,keep_idx)])
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

  covariates[,c(regress_idx,indiv_idx)] <- scale(covariates[,c(regress_idx,indiv_idx)],
                                                 center = T,
                                                 scale = T)

  ##########

  # step 2: poisson GLMS

  coef_mat <- sapply(1:p, function(j){
    if(verbose >= 1 && p > 10 && j %% floor(p/10) == 0) cat('*')

    df <- as.data.frame(cbind(y = mat[,j], covariates))
    colnames(df)[1] <- "tmp"
    glm_fit <- stats::glm(tmp ~ . - 1, data = df)
    ## [[note to self: can be improved using the stat helper functions]]
    stats::coef(glm_fit)
  })

  #########

  # step 3: SVD
  dat_transform <- log1p(dat)
  nat_mat <- tcrossprod(covariates, coef_mat)
  residual_mat <- dat_transform - nat_mat

  svd_res <- .svd_truncated(residual_mat,
                            K = k,
                            symmetric = F,
                            rescale = F,
                            mean_vec = NULL,
                            sd_vec = NULL,
                            K_full_rank = F)
  x_init <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_init <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  rownames(x_init) <- rownames(dat)
  rownames(y_init) <- colnames(dat)

  structure(list(x_mat = x_init, y_mat = y_init, b_mat = b_init,
                 covariates = covariates,
                 nuisance_param_vec = nuisance_init,
                 offset_vec = offset_vec),
            class = "eSVD")
}
