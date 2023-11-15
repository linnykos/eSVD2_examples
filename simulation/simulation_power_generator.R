data_generator_nat_mat <- function(
    cell_latent_gaussian_mean,
    cell_latent_gaussian_covariance,
    gene_library_repeating_vec,
    gene_covariate_coefficient_proportion_mat,
    gene_covariate_coefficient_size,
    gene_nuisance_values,
    gene_nuisance_proporition_mat,
    gene_null_casecontrol_name,
    gene_null_casecontrol_proportion,
    gene_null_casecontrol_size,
    gene_null_latent_gaussian_noise,
    gene_null_nuisance_proportion,
    gene_num_mixed_membership,
    gene_num_null,
    gene_num_per_topic,
    gene_topic_latent_gaussian_noise,
    gene_topic_simplex,
    gene_topic_casecontrol_name,
    gene_topic_casecontrol_proportion,
    gene_topic_casecontrol_size,
    individual_covariates,
    individual_case_control_variable,
    individual_num_cells,
    gene_intercept_global_shift = 0
){
  stopifnot(nrow(gene_topic_simplex) == length(cell_latent_gaussian_mean))

  num_indiv <- nrow(individual_covariates)
  n <- individual_num_cells*num_indiv

  res <- .form_covariate_mat(individual_case_control_variable = individual_case_control_variable,
                             individual_covariates = individual_covariates,
                             individual_num_cells = individual_num_cells)
  covariates <- res$covariates
  case_individuals <- res$case_individuals
  control_individuals <- res$control_individuals
  individual_vec <- res$individual_vec

  x_mat <- .form_x_mat(
    cell_latent_gaussian_covariance = cell_latent_gaussian_covariance,
    cell_latent_gaussian_mean = cell_latent_gaussian_mean,
    covariates = covariates,
    individual_case_control_variable = individual_case_control_variable,
    individual_vec = individual_vec,
    n = n
  )

  # now work on the genes
  res <- .form_y_mat(gene_null_latent_gaussian_noise = gene_null_latent_gaussian_noise,
                     gene_num_null = gene_num_null,
                     gene_num_mixed_membership = gene_num_mixed_membership,
                     gene_num_per_topic = gene_num_per_topic,
                     gene_topic_latent_gaussian_noise = gene_topic_latent_gaussian_noise,
                     gene_topic_simplex = gene_topic_simplex)
  y_mat <- res$y_mat
  gene_labeling <- res$gene_labeling

  res <- .form_z_mat(covariates_colnames = colnames(covariates),
                     gene_covariate_coefficient_proportion_mat = gene_covariate_coefficient_proportion_mat,
                     gene_covariate_coefficient_size = gene_covariate_coefficient_size,
                     gene_intercept_global_shift = gene_intercept_global_shift,
                     gene_null_casecontrol_name = gene_null_casecontrol_name,
                     gene_null_casecontrol_proportion = gene_null_casecontrol_proportion,
                     gene_null_casecontrol_size = gene_null_casecontrol_size,
                     gene_topic_casecontrol_name = gene_topic_casecontrol_name,
                     gene_topic_casecontrol_proportion = gene_topic_casecontrol_proportion,
                     gene_topic_casecontrol_size = gene_topic_casecontrol_size,
                     gene_labeling = gene_labeling,
                     p = nrow(y_mat))
  z_mat <- res$z_mat; gene_labeling2 <- res$gene_labeling2

  nuisance_vec <- .form_nuisance_vec(gene_labeling2 = gene_labeling2,
                                     gene_nuisance_proporition_mat = gene_nuisance_proporition_mat,
                                     gene_nuisance_values = gene_nuisance_values,
                                     gene_null_casecontrol_name = gene_null_casecontrol_name,
                                     gene_topic_casecontrol_name = gene_topic_casecontrol_name,
                                     p = nrow(y_mat))

  p <- nrow(y_mat)
  gene_library_vec <- rep(gene_library_repeating_vec, times = ceiling(p/length(gene_library_repeating_vec)))
  gene_library_vec <- gene_library_vec[1:p]

  # append all the names
  rownames(x_mat) <- paste0("cell_", 1:nrow(x_mat))
  rownames(covariates) <- rownames(x_mat)
  rownames(y_mat) <- paste0("gene_", 1:p)
  names(nuisance_vec) <- rownames(y_mat)
  names(gene_labeling) <- rownames(y_mat)
  names(gene_labeling2) <- rownames(y_mat)
  names(gene_library_vec) <- rownames(y_mat)

  list(case_individuals = case_individuals,
       control_individuals = control_individuals,
       covariates = covariates,
       gene_labeling = gene_labeling,
       gene_labeling2 = gene_labeling2,
       gene_library_repeating_vec = gene_library_repeating_vec,
       individual_case_control_variable = individual_case_control_variable,
       individual_vec = individual_vec,
       nuisance_vec = nuisance_vec,
       x_mat = x_mat,
       y_mat = y_mat,
       z_mat = z_mat)
}

data_signal_enhancer <- function(input_obj,
                                 gene_case_control_none_size_sd = 0.3,
                                 global_shift = 0){
  case_individuals <- input_obj$case_individuals
  control_individuals <- input_obj$control_individuals
  covariates <- input_obj$covariates
  gene_labeling <- input_obj$gene_labeling
  gene_labeling2 <- input_obj$gene_labeling2
  individual_case_control_variable <- input_obj$individual_case_control_variable
  individual_vec <- input_obj$individual_vec
  x_mat <- input_obj$x_mat
  y_mat <- input_obj$y_mat
  z_mat <- input_obj$z_mat

  idx <- intersect(which(gene_labeling2 == "none"), which(z_mat[,individual_case_control_variable] == 0))
  z_mat[idx,individual_case_control_variable] <- stats::rnorm(length(idx), mean = 0, sd = gene_case_control_none_size_sd)

  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)
  nat_mat <- nat_mat + global_shift

  subset_case_indivuals <- sample(case_individuals, round(length(case_individuals)/4))
  size_addition <- runif(length(subset_case_indivuals), min = 0, max = 2)
  size_multiplier1 <- runif(length(subset_case_indivuals), min = 0, max = 1.5)

  subset_control_indivuals <- sample(control_individuals, round(length(control_individuals)/4))
  size_subtraction <- runif(length(subset_control_indivuals), min = 0, max = 2)
  size_multiplier2 <- runif(length(subset_control_indivuals), min = 0, max = 1.5)

  # start by exaggerating the strong signals
  gene_idx <- which(gene_labeling2 == "strong-positive")
  for(j in gene_idx){
    for(kk in 1:length(subset_case_indivuals)){
      individual_name <- subset_case_indivuals[kk]
      cell_idx <- which(individual_vec %in% individual_name)
      tmp <- nat_mat[cell_idx,j]
      nat_mat[cell_idx,j] <- pmax(tmp+size_addition[kk], tmp*size_multiplier1[kk])
    }

    for(kk in 1:length(subset_control_indivuals)){
      individual_name <- subset_control_indivuals[kk]
      cell_idx <- which(individual_vec %in% individual_name)
      tmp <- nat_mat[cell_idx,j]
      nat_mat[cell_idx,j] <- pmin(tmp-size_subtraction[kk], tmp*size_multiplier2[kk])
    }
  }

  gene_idx <- which(gene_labeling2 == "strong-negative")
  for(j in gene_idx){
    for(kk in 1:length(subset_case_indivuals)){
      individual_name <- subset_case_indivuals[kk]
      cell_idx <- which(individual_vec %in% individual_name)
      tmp <- nat_mat[cell_idx,j]
      nat_mat[cell_idx,j] <- pmin(tmp-size_addition[kk], tmp*size_multiplier1[kk])
    }

    for(kk in 1:length(subset_control_indivuals)){
      individual_name <- subset_control_indivuals[kk]
      cell_idx <- which(individual_vec %in% individual_name)
      tmp <- nat_mat[cell_idx,j]
      nat_mat[cell_idx,j] <- pmax(tmp+size_subtraction[kk], tmp*size_multiplier2[kk])
    }
  }

  # shrink none and weak signals towards the strong ones
  tab_mat <- table(gene_labeling, gene_labeling2)
  topic_names <- rownames(tab_mat)[grep("topic", rownames(tab_mat))]
  for(topic_name in topic_names){
    strong_pos_idx <- intersect(which(gene_labeling2 == "strong-positive"), which(gene_labeling == topic_name))
    strong_neg_idx <- intersect(which(gene_labeling2 == "strong-negative"), which(gene_labeling == topic_name))
    weak_pos_idx <- intersect(which(gene_labeling2 == "weak-positive"), which(gene_labeling == topic_name))
    weak_neg_idx <- intersect(which(gene_labeling2 == "weak-negative"), which(gene_labeling == topic_name))
    none_idx <- intersect(which(gene_labeling2 == "none"), which(gene_labeling == topic_name))

    none_number <- ceiling(length(none_idx)*.4)
    none_idx_posShrink <- sample(none_idx, none_number)
    none_idx <- setdiff(none_idx, none_idx_posShrink)
    none_idx_negShrink <- sample(none_idx, none_number)
    weak_pos_idx_posShrink <- sample(weak_pos_idx, ceiling(length(weak_pos_idx)*.3))
    weak_neg_idx_posShrink <- sample(weak_neg_idx, ceiling(length(weak_neg_idx)*.3))

    # denoise the strong genes
    if(length(strong_pos_idx) > 0){
      pos_nat_mat <- apply(nat_mat[,strong_pos_idx], 2, function(x){
        df <- data.frame(cbind(x, covariates))
        colnames(df)[1] <- "y"
        lm_res <- stats::lm(y ~ . - 1, data = df)
        stats::residuals(lm_res)
      })

      nat_mat <- .natural_shrinkage(nat_mat = nat_mat,
                                    shrinkage_idx = none_idx_posShrink,
                                    shrinkage_val = 0.7,
                                    target_mat = pos_nat_mat)
      nat_mat <- .natural_shrinkage(nat_mat = nat_mat,
                                    shrinkage_idx = weak_pos_idx_posShrink,
                                    shrinkage_val = 0.95,
                                    target_mat = pos_nat_mat)
    }

    if(length(strong_neg_idx) > 0){
      neg_nat_mat <- apply(nat_mat[,strong_neg_idx,drop = F], 2, function(x){
        df <- data.frame(cbind(x, covariates))
        colnames(df)[1] <- "y"
        lm_res <- stats::lm(y ~ . - 1, data = df)
        stats::residuals(lm_res)
      })

      nat_mat <- .natural_shrinkage(nat_mat = nat_mat,
                                    shrinkage_idx = none_idx_negShrink,
                                    shrinkage_val = 0.7,
                                    target_mat = neg_nat_mat)
      nat_mat <- .natural_shrinkage(nat_mat = nat_mat,
                                    shrinkage_idx = weak_neg_idx_posShrink,
                                    shrinkage_val = 0.95,
                                    target_mat = neg_nat_mat)
    }
  }

  list(case_individuals = input_obj$case_individuals,
       control_individuals = input_obj$control_individuals,
       covariates = covariates,
       gene_labeling = gene_labeling,
       gene_labeling2 = gene_labeling2,
       gene_library_repeating_vec = input_obj$gene_library_repeating_vec,
       individual_case_control_variable = input_obj$individual_case_control_variable,
       individual_vec = input_obj$individual_vec,
       nat_mat = nat_mat,
       nuisance_vec = input_obj$nuisance_vec,
       x_mat = input_obj$x_mat,
       y_mat = input_obj$y_mat,
       z_mat = input_obj$z_mat)
}

data_generator_obs_mat <- function(input_obj){
  gene_library_repeating_vec <- input_obj$gene_library_repeating_vec
  nuisance_vec <- input_obj$nuisance_vec
  nat_mat <- input_obj$nat_mat
  covariates <- input_obj$covariates
  z_mat <- input_obj$z_mat

  n <- nrow(nat_mat); p <- ncol(nat_mat)

  gamma_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    gamma_mat[,j] <- stats::rgamma(
      n = n,
      shape = exp(nat_mat[,j])*nuisance_vec[j],
      rate = nuisance_vec[j])
  }
  gamma_mat <- pmin(gamma_mat, 1000)

  gene_library_vec <- rep(gene_library_repeating_vec, times = ceiling(p/length(gene_library_repeating_vec)))
  gene_library_vec <- gene_library_vec[1:p]
  obs_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    obs_mat[,j] <- stats::rpois(n = n,
                                lambda = gene_library_vec[j]*gamma_mat[,j])
  }

  covariates <- cbind(covariates, log1p(rowMeans(obs_mat)))
  colnames(covariates)[ncol(covariates)] <- "Log_UMI"
  z_mat <- cbind(z_mat, 0)
  colnames(z_mat)[ncol(z_mat)] <- "Log_UMI"

  rownames(gamma_mat) <- rownames(nat_mat)
  rownames(obs_mat) <- rownames(nat_mat)
  colnames(gamma_mat) <- colnames(nat_mat)
  colnames(obs_mat) <- colnames(nat_mat)

  list(case_individuals = input_obj$case_individuals,
       control_individuals = input_obj$control_individuals,
       covariates = covariates,
       gamma_mat = gamma_mat,
       gene_labeling = input_obj$gene_labeling,
       gene_labeling2 = input_obj$gene_labeling2,
       gene_library_vec = gene_library_vec,
       individual_case_control_variable = input_obj$individual_case_control_variable,
       individual_vec = input_obj$individual_vec,
       nat_mat = nat_mat,
       nuisance_vec = nuisance_vec,
       obs_mat = obs_mat,
       x_mat = input_obj$x_mat,
       y_mat = input_obj$y_mat,
       z_mat = input_obj$z_mat)
}


.compute_population_quantities <- function(input_obj){
  case_individuals <- input_obj$case_individuals
  control_individuals <- input_obj$control_individuals
  covariates <- input_obj$covariates
  individual_case_control_variable <- input_obj$individual_case_control_variable
  individual_vec <- input_obj$individual_vec
  nuisance_vec <- input_obj$nuisance_vec
  x_mat <- input_obj$x_mat
  y_mat <- input_obj$y_mat
  z_mat <- input_obj$z_mat

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates[,individual_case_control_variable], z_mat[,individual_case_control_variable])
  nat_mat <- nat_mat1 + nat_mat2
  mean_mat <- exp(nat_mat)

  case_idx <- which(covariates[,individual_case_control_variable] == 1)
  control_idx <- which(covariates[,individual_case_control_variable] == 0)
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
  true_teststat_vec <- res$teststat_vec

  tmp <- eSVD2:::.determine_individual_indices(case_individuals = case_individuals,
                                               control_individuals = control_individuals,
                                               individual_vec = individual_vec)
  all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
  avg_mat <- eSVD2:::.construct_averaging_matrix(idx_list = all_indiv_idx,
                                                 n = nrow(mean_mat))
  avg_posterior_mean_mat <- as.matrix(avg_mat %*% mean_mat)
  avg_posterior_var_mat <- as.matrix(avg_mat %*% var_mat)

  case_row_idx <- 1:length(case_individuals)
  control_row_idx <- (length(case_individuals)+1):nrow(avg_posterior_mean_mat)
  case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_row_idx,,drop = F])
  control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_row_idx,,drop = F])
  case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[case_row_idx,,drop = F],
    avg_posterior_var_mat = avg_posterior_var_mat[case_row_idx,,drop = F]
  )
  control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[control_row_idx,,drop = F],
    avg_posterior_var_mat = avg_posterior_var_mat[control_row_idx,,drop = F]
  )
  n1 <- length(case_individuals); n2 <- length(control_individuals)
  numerator_vec <- (case_gaussian_var/n1 + control_gaussian_var/n2)^2
  denominator_vec <- (case_gaussian_var/n1)^2/(n1-1) + (control_gaussian_var/n2)^2/(n2-1)
  df_vec <- numerator_vec/denominator_vec
  names(df_vec) <- names(case_gaussian_var)
  p <- length(true_teststat_vec)
  gaussian_teststat <- sapply(1:p, function(j){
    qnorm(pt(true_teststat_vec[j], df = df_vec[j]))
  })

  fdr_res <- eSVD2:::multtest(gaussian_teststat)
  true_fdr_vec <- fdr_res$fdr_vec
  names(true_fdr_vec) <- names(gaussian_teststat)
  true_null_mean <- fdr_res$null_mean
  true_null_sd <- fdr_res$null_sd
  true_logpvalue_vec <- sapply(gaussian_teststat, function(x){
    if(x < true_null_mean) {
      Rmpfr::pnorm(x, mean = true_null_mean, sd = true_null_sd, log.p = T)
    } else {
      Rmpfr::pnorm(true_null_mean - (x-true_null_mean), mean = true_null_mean, sd = true_null_sd, log.p = T)
    }
  })
  true_logpvalue_vec <- -(true_logpvalue_vec/log(10) + log10(2))
  names(true_logpvalue_vec) <- names(gaussian_teststat)

  list(true_fdr_vec = true_fdr_vec,
       true_logpvalue_vec = true_logpvalue_vec,
       true_null_mean = true_null_mean,
       true_null_sd = true_null_sd,
       true_teststat_vec = true_teststat_vec)
}

####################

.form_covariate_mat <- function(individual_case_control_variable,
                                individual_covariates,
                                individual_num_cells){
  num_indiv <- nrow(individual_covariates)
  covariates <- do.call(rbind, lapply(1:num_indiv, function(i){
    matrix(rep(as.numeric(individual_covariates[i,]), each = individual_num_cells),
           nrow = individual_num_cells, ncol = ncol(individual_covariates))
  }))
  covariates <- cbind(1, covariates)
  colnames(covariates) <- c("Intercept", colnames(individual_covariates))

  individual_vec <- paste0("indiv_", 1:num_indiv)
  case_individuals <- individual_vec[which(individual_covariates[,individual_case_control_variable] == 1)]
  control_individuals <- individual_vec[which(individual_covariates[,individual_case_control_variable] == 0)]
  individual_vec_full <- rep(paste0("indiv_", 1:num_indiv), each = individual_num_cells)

  list(case_individuals = case_individuals,
       control_individuals = control_individuals,
       covariates = covariates,
       individual_vec = individual_vec_full)
}

.form_x_mat <- function(cell_latent_gaussian_covariance,
                        cell_latent_gaussian_mean,
                        covariates,
                        individual_case_control_variable,
                        individual_vec,
                        n){
  x_mat <- MASS::mvrnorm(n,
                         mu = cell_latent_gaussian_mean,
                         Sigma = cell_latent_gaussian_covariance)

  x_mat
}

.form_y_mat <- function(gene_null_latent_gaussian_noise,
                        gene_num_null,
                        gene_num_mixed_membership,
                        gene_num_per_topic,
                        gene_topic_latent_gaussian_noise,
                        gene_topic_simplex){
  r <- ncol(gene_topic_simplex)
  k <- nrow(gene_topic_simplex)

  gene_latent_embedding_list <- vector("list", length = r + 2)
  for(i in 1:r){
    mat <- matrix(rep(gene_topic_simplex[,i], each = gene_num_per_topic),
                  nrow = gene_num_per_topic, ncol = k)
    mat <- mat + MASS::mvrnorm(n = gene_num_per_topic,
                               mu = rep(0, k),
                               Sigma = gene_topic_latent_gaussian_noise)
    gene_latent_embedding_list[[i]] <- mat
  }

  mat <- t(sapply(1:gene_num_mixed_membership, function(i){
    proportion_vec <- runif(r)
    proportion_vec <- proportion_vec/sum(proportion_vec)
    gene_topic_simplex %*% proportion_vec
  }))
  mat <- mat + MASS::mvrnorm(n = gene_num_per_topic,
                             mu = rep(0, k),
                             Sigma = gene_topic_latent_gaussian_noise)
  gene_latent_embedding_list[[r+1]] <- mat
  gene_latent_embedding_list[[r+2]] <- MASS::mvrnorm(n = gene_num_null,
                                                     mu = rep(0, k),
                                                     Sigma = gene_null_latent_gaussian_noise)
  y_mat <- do.call(rbind, gene_latent_embedding_list)
  gene_labeling <- unlist(lapply(1:length(gene_latent_embedding_list), function(i){
    if(i <= r) return(rep(paste0("topic_", i), nrow(gene_latent_embedding_list[[i]])))
    if(i == r+1) return(rep("topic_mixture", nrow(gene_latent_embedding_list[[i]])))
    if(i == r+2) return(rep("null", nrow(gene_latent_embedding_list[[i]])))
  }))

  list(gene_labeling = gene_labeling,
       y_mat = y_mat)
}

.exact_sampler <- function(target_length,
                           value_proportion,
                           value_vec){
  if(target_length == 0) return(numeric(0))

  value_proportion <- value_proportion/sum(value_proportion)
  tmp <- unlist(lapply(1:length(value_vec), function(i){
    rep(value_vec[i], ceiling(target_length*value_proportion[i]))
  }))
  if(length(tmp) >= target_length) {
    tmp <- tmp[sample(1:length(tmp), size = target_length)]
  }
  sample(tmp)
}

.form_z_mat <- function(covariates_colnames,
                        gene_covariate_coefficient_proportion_mat,
                        gene_covariate_coefficient_size,
                        gene_intercept_global_shift,
                        gene_null_casecontrol_name,
                        gene_null_casecontrol_proportion,
                        gene_null_casecontrol_size,
                        gene_topic_casecontrol_name,
                        gene_topic_casecontrol_proportion,
                        gene_topic_casecontrol_size,
                        gene_labeling,
                        p){
  r <- length(covariates_colnames)
  z_mat <- matrix(0, nrow = p, ncol = r)
  colnames(z_mat) <- covariates_colnames
  z_mat[,"Intercept"] <- gene_intercept_global_shift

  # now for case-control
  gene_labeling2 <- rep("NA", p)
  gene_in_topics <- grep("topic", gene_labeling)
  z_mat[gene_in_topics, individual_case_control_variable] <- .exact_sampler(
    target_length = length(gene_in_topics),
    value_proportion = gene_topic_casecontrol_proportion,
    value_vec = gene_topic_casecontrol_size)
  gene_labeling2[gene_in_topics] <- plyr::mapvalues(
    x = as.character(z_mat[gene_in_topics, individual_case_control_variable]),
    from = as.character(gene_topic_casecontrol_size),
    to = gene_topic_casecontrol_name
  )

  gene_in_null <- grep("null", gene_labeling)
  z_mat[gene_in_null, individual_case_control_variable] <- .exact_sampler(
    target_length = length(gene_in_null),
    value_proportion = gene_null_casecontrol_proportion,
    value_vec = gene_null_casecontrol_size)
  gene_labeling2[gene_in_null] <- plyr::mapvalues(
    x = as.character(z_mat[gene_in_null, individual_case_control_variable]),
    from = as.character(gene_null_casecontrol_size),
    to = gene_null_casecontrol_name
  )

  # now for the covariates
  gene_casecontrol_name <- sort(unique(gene_topic_casecontrol_name, gene_null_casecontrol_name))
  covariates_other <- colnames(z_mat)[which(!colnames(z_mat) %in% c("Intercept", individual_case_control_variable))]
  for(variable in covariates_other){
    r <- length(gene_casecontrol_name)
    for(j in 1:r){
      gene_in_size <- grep(gene_casecontrol_name[j], gene_labeling2)
      z_mat[gene_in_size,variable] <- .exact_sampler(
        target_length = length(gene_in_size),
        value_proportion = gene_covariate_coefficient_proportion_mat[,paste0("cc_size:", gene_casecontrol_name[j])],
        value_vec = gene_covariate_coefficient_size
      )
    }
  }

  list(gene_labeling2 = gene_labeling2,
       z_mat = z_mat)
}

.form_nuisance_vec <- function(gene_labeling2,
                               gene_nuisance_proporition_mat,
                               gene_nuisance_values,
                               gene_null_casecontrol_name,
                               gene_topic_casecontrol_name,
                               p){
  gene_casecontrol_name <- sort(unique(gene_topic_casecontrol_name, gene_null_casecontrol_name))
  nuisance_vec <- rep(0, p)
  stopifnot(all(unique(gene_labeling2) %in% gene_casecontrol_name))

  r <- length(gene_casecontrol_name)
  for(j in 1:r){
    gene_in_size <- grep(gene_casecontrol_name[j], gene_labeling2)
    nuisance_vec[gene_in_size] <- .exact_sampler(
      target_length = length(gene_in_size),
      value_proportion = gene_nuisance_proporition_mat[,paste0("cc_size:", gene_casecontrol_name[j])],
      value_vec = gene_nuisance_values
    )
  }

  nuisance_vec
}

.natural_shrinkage <- function(nat_mat,
                               shrinkage_idx,
                               shrinkage_val,
                               target_mat){
  stopifnot(shrinkage_val >= 0, shrinkage_val <= 1)

  cor_mat <- stats::cor(nat_mat[,shrinkage_idx,drop=F], target_mat)
  for(j in 1:length(shrinkage_idx)){
    k <- which.max(cor_mat[j,])

    y_vec <- nat_mat[,shrinkage_idx[j]]
    df <- data.frame(cbind(y_vec, target_mat[,k]))
    colnames(df)[1] <- "y"
    lm_res <- stats::lm(y ~ . - 1, data = df)
    resid_vec <- stats::residuals(lm_res)

    y_vec <- (1-shrinkage_val)*y_vec + shrinkage_val*resid_vec
    nat_mat[,shrinkage_idx[j]] <- y_vec
  }

  nat_mat
}
