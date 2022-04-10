compute_posterior <- function(alpha_max, # default: 50
                              case_control_variable,
                              esvd_res,
                              mat,
                              nuisance_lower_quantile, # default: 0.01
                              nuisance_vec){
  offset_var <- setdiff(colnames(esvd_res$covariates), case_control_variable)

  nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
  nat_mat2 <- tcrossprod(esvd_res$covariates[,case_control_variable,drop = F],
                         esvd_res$b_mat[,case_control_variable,drop = F])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)
  library_mat <- exp(tcrossprod(
    esvd_res$covariates[,offset_var],
    esvd_res$b_mat[,offset_var]
  ))

  nuisance_vec <- pmax(nuisance_vec, quantile(nuisance_vec, probs = nuisance_lower_quantile))
  Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
                 STATS = nuisance_vec, FUN = "*")
  Alpha <- pmin(Alpha, alpha_max)
  AplusAlpha <- mat + Alpha
  SplusBeta <- sweep(library_mat, MARGIN = 2,
                     STATS = nuisance_vec, FUN = "+")
  posterior_mean_mat <- AplusAlpha/SplusBeta
  posterior_var_mat <- AplusAlpha/SplusBeta^2

  rownames(posterior_mean_mat) <- rownames(mat)
  rownames(posterior_var_mat) <- rownames(mat)
  colnames(posterior_mean_mat) <- colnames(mat)
  colnames(posterior_var_mat) <- colnames(mat)

  list(posterior_mean_mat = posterior_mean_mat,
       posterior_var_mat = posterior_var_mat)
}

compute_test_statistics <- function(case_individuals,
                                    control_individuals,
                                    covariate_individual,
                                    metadata,
                                    posterior_mean_mat,
                                    posterior_var_mat,
                                    verbose = F){
  p <- ncol(posterior_mean_mat)
  individual_stats <- lapply(1:p, function(j){
    if(verbose && p > 10 && j %% floor(p/10) == 0) cat('*')

    # next find the cells, then compute one gaussian per individual
    case_gaussians <- sapply(case_individuals, function(indiv){
      cell_names <- rownames(metadata)[which(metadata[,covariate_individual] == indiv)]
      cell_idx <- which(rownames(mat) %in% cell_names)

      mean_val <- mean(posterior_mean_mat[cell_idx,j])
      var_val <- mean(posterior_var_mat[cell_idx,j])
      c(mean_val = mean_val, var_val = var_val)
    })

    control_gaussians <- sapply(control_individuals, function(indiv){
      cell_names <- rownames(metadata)[which(metadata[,covariate_individual] == indiv)]
      cell_idx <- which(rownames(mat) %in% cell_names)

      mean_val <- mean(posterior_mean_mat[cell_idx,j])
      var_val <- mean(posterior_var_mat[cell_idx,j])
      c(mean_val = mean_val, var_val = var_val)
    })

    list(case_gaussians = case_gaussians,
         control_gaussians = control_gaussians)
  })

  # see https://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
  group_stats <- lapply(1:p, function(j){
    case_gaussians <- individual_stats[[j]]$case_gaussians
    control_gaussians <- individual_stats[[j]]$control_gaussians

    case_gaussian <- list(mean_val = mean(case_gaussians[1,]),
                          var_val = mean(case_gaussians[2,]) + mean(case_gaussians[1,]^2) - (mean(case_gaussians[1,]))^2,
                          n = ncol(case_gaussians))
    control_gaussian <- list(mean_val = mean(control_gaussians[1,]),
                             var_val = mean(control_gaussians[2,]) + mean(control_gaussians[1,]^2) - (mean(control_gaussians[1,]))^2,
                             n = ncol(control_gaussians))

    list(case_gaussian = case_gaussian,
         control_gaussian = control_gaussian)
  })

  teststat_vec <- sapply(1:p, function(j){
    case_gaussian <- group_stats[[j]]$case_gaussian
    control_gaussian <- group_stats[[j]]$control_gaussian

    n1 <- control_gaussian$n; n2 <- case_gaussian$n
    mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
    cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$var_val

    combined_cov <- cov1/n1 + cov2/n2
    (mean2 - mean1)/sqrt(combined_cov)
  })
  names(teststat_vec) <- colnames(posterior_mean_mat)

  teststat_vec
}

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


