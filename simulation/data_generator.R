data_generator <- function(
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
    natural_param_max_quant,
    bool_include_extra_signal = T,
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
    bool_include_extra_signal = bool_include_extra_signal,
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

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates, z_mat)
  nat_mat <- nat_mat1 + nat_mat2
  nat_colQuant <- apply(nat_mat, 2, stats::quantile, probs = 0.9)
  idx <- which(nat_colQuant >= natural_param_max_quant)
  if(length(idx) > 0){
    for(j in idx){
      z_mat[j,"Intercept"] <- natural_param_max_quant - nat_colQuant[j] + z_mat[j,"Intercept"]
    }
  }

  nuisance_vec <- .form_nuisance_vec(gene_labeling2 = gene_labeling2,
                                     gene_nuisance_proporition_mat = gene_nuisance_proporition_mat,
                                     gene_nuisance_values = gene_nuisance_values,
                                     gene_null_casecontrol_name = gene_null_casecontrol_name,
                                     gene_topic_casecontrol_name = gene_topic_casecontrol_name,
                                     p = ncol(nat_mat))

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates, z_mat)
  nat_mat <- nat_mat1 + nat_mat2
  nat_mat[nat_mat >= natural_param_max_quant] <- natural_param_max_quant

  n <- nrow(nat_mat); p <- ncol(nat_mat)
  gamma_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    gamma_mat[,j] <- stats::rgamma(
      n = n,
      shape = exp(nat_mat[,j])*nuisance_vec[j],
      rate = nuisance_vec[j])
  }

  gene_library_vec <- rep(gene_library_repeating_vec, times = ceiling(p/length(gene_library_repeating_vec)))
  gene_library_vec <- gene_library_vec[1:p]
  obs_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    obs_mat[,j] <- stats::rpois(n = n,
                                lambda = gene_library_vec[j]*gamma_mat[,j])
  }

  ## sparsity_downsampling has no effect right now
  covariates <- cbind(covariates, log1p(rowMeans(obs_mat)))
  colnames(covariates)[ncol(covariates)] <- "Log_UMI"
  z_mat <- cbind(z_mat, 0)
  colnames(z_mat)[ncol(z_mat)] <- "Log_UMI"

  # append all the names
  rownames(x_mat) <- paste0("cell_", 1:n)
  rownames(gamma_mat) <-  rownames(x_mat)
  rownames(obs_mat) <- rownames(x_mat)
  rownames(covariates) <- rownames(x_mat)
  rownames(y_mat) <- paste0("gene_", 1:p)
  colnames(gamma_mat) <-  rownames(y_mat)
  colnames(obs_mat) <- rownames(y_mat)
  names(nuisance_vec) <- rownames(y_mat)
  names(gene_labeling) <- rownames(y_mat)
  names(gene_labeling2) <- rownames(y_mat)
  names(gene_library_vec) <- rownames(y_mat)

  list(case_individuals = case_individuals,
       control_individuals = control_individuals,
       covariates = covariates,
       gene_labeling = gene_labeling,
       gene_labeling2 = gene_labeling2,
       gene_library_vec = gene_library_vec,
       individual_vec = individual_vec,
       nuisance_vec = nuisance_vec,
       obs_mat = obs_mat,
       x_mat = x_mat,
       y_mat = y_mat,
       z_mat = z_mat)
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

.form_x_mat <- function(bool_include_extra_signal,
                        cell_latent_gaussian_covariance,
                        cell_latent_gaussian_mean,
                        covariates,
                        individual_case_control_variable,
                        individual_vec,
                        n){
  if(bool_include_extra_signal){
    tmp <- MASS::mvrnorm(n,
                         mu = cell_latent_gaussian_mean[-1],
                         Sigma = cell_latent_gaussian_covariance[-1,-1])
    tab <- table(individual_vec, covariates[,individual_case_control_variable])
    case_indiv <- rownames(tab)[which(tab[,"1"] != 0)]
    control_indiv <- rownames(tab)[which(tab[,"0"] != 0)]

    vec <- rep(0, n)
    vec[individual_vec %in% control_indiv] <- -1
    case_indiv_firsthalf <- case_indiv[1:ceiling(length(case_indiv)/2)]
    vec[individual_vec %in% case_indiv_firsthalf] <- 1/2
    case_indiv_secondhalf <- setdiff(case_indiv, case_indiv_firsthalf)
    for(indiv in case_indiv_secondhalf){
      idx <- which(individual_vec == indiv)
      vec[idx[1:ceiling(length(idx)/2)]] <- -1
    }
    x_mat <- cbind(vec, tmp)

  } else {
    x_mat <- MASS::mvrnorm(n,
                           mu = cell_latent_gaussian_mean,
                           Sigma = cell_latent_gaussian_covariance)
  }

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
