data_generator <- function(
    cell_latent_gaussian_mean,
    cell_latent_gaussian_covariance,
    gene_library_repeating_vec,
    gene_covariate_coefficient_proportion,
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
    sparsity_downsampling
){
  stopifnot(all(dim(gene_topic_nuisance_proporition_mat) ==
                  c(length(gene_nuisance_values), length(gene_topic_casecontrol_size))),
            nrow(gene_topic_simplex) == length(cell_latent_gaussian_mean),
            length(gene_covariate_coefficient_proportion) == length(gene_covariate_coefficient_size))

  num_indiv <- nrow(individual_covariates)
  n <- individual_num_cells*num_indiv

  x_mat <- .form_x_mat(cell_latent_gaussian_covariance = cell_latent_gaussian_covariance,
                       cell_latent_gaussian_mean = cell_latent_gaussian_mean,
                       n = n)
  covariates <- .form_covariate_mat(individual_covariates = individual_covariates,
                                    individual_num_cells = individual_num_cells)

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
                     gene_covariate_coefficient_proportion = gene_covariate_coefficient_proportion,
                     gene_covariate_coefficient_size = gene_covariate_coefficient_size,
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
      z_mat[j,"Intercept"] <- natural_param_max_quant - nat_colQuant[j]
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
      shape = exp(nat_mat[,j])/nuisance_vec[j],
      rate = 1/nuisance_vec[j])
  }

  gene_library_vec <- rep(gene_library_repeating_vec, times = ceiling(p/length(gene_library_repeating_vec)))
  gene_library_vec <- gene_library_vec[1:p]
  obs_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    obs_mat[,j] <- stats::rpois(n = n,
                                lambda = gene_library_vec[j]*gamma_mat[,j])
  }

  list(covariates = covariates,
       gene_labeling = gene_labeling,
       gene_labeling2 = gene_labeling2,
       nuisance_vec = nuisance_vec,
       obs_mat = obs_mat,
       x_mat = x_mat,
       y_mat = y_mat,
       z_mat = z_mat)
}

####################

.form_x_mat <- function(cell_latent_gaussian_covariance,
                        cell_latent_gaussian_mean,
                        n){
  MASS::mvrnorm(n,
                mu = cell_latent_gaussian_mean,
                Sigma = cell_latent_gaussian_covariance)
}

.form_covariate_mat <- function(individual_covariates,
                                individual_num_cells){
  num_indiv <- nrow(individual_covariates)
  covariates <- do.call(rbind, lapply(1:num_indiv, function(i){
    matrix(rep(as.numeric(individual_covariates[i,]), each = individual_num_cells),
           nrow = individual_num_cells, ncol = ncol(individual_covariates))
  }))
  covariates <- cbind(1, covariates)
  colnames(covariates) <- c("Intercept", colnames(individual_covariates))

  covariates
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
                        gene_covariate_coefficient_proportion,
                        gene_covariate_coefficient_size,
                        gene_null_casecontrol_name,
                        gene_null_casecontrol_proportion,
                        gene_null_casecontrol_size,
                        gene_topic_casecontrol_name,
                        gene_topic_casecontrol_proportion,
                        gene_topic_casecontrol_size,
                        gene_labeling,
                        p){
  p <- nrow(y_mat)
  r <- length(covariates_colnames)
  z_mat <- matrix(0, nrow = p, ncol = r)
  colnames(z_mat) <- covariates_colnames
  covariates_other <- colnames(z_mat)[which(!colnames(z_mat) %in% c("Intercept", individual_case_control_variable))]
  for(variable in covariates_other){
    z_mat[,variable] <- .exact_sampler(target_length = nrow(z_mat),
                                       value_proportion = gene_covariate_coefficient_proportion,
                                       value_vec = gene_covariate_coefficient_size)
  }

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
