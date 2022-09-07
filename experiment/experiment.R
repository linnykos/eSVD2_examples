rm(list=ls())
set.seed(10)

k <- 10
cell_latent_gaussian_mean = rep(0, k)
cov_mat <- matrix(0, k, k)
cov_mat[1:3,1:3] <- 0.9
cov_mat[4:7,4:7] <- 0.8
cov_mat[8:10,8:10] <- 0.6
diag(cov_mat) <- 1
cell_latent_gaussian_covariance = cov_mat

##############

gene_null_casecontrol_name <- c("none", "strong-negative", "strong-positive")
gene_null_casecontrol_proportion <- c(0.95, 0.03, 0.02)
gene_null_casecontrol_size <- c(0, -5, 2)
gene_null_latent_gaussian_noise <- 0.5*diag(k)

gene_num_mixed_membership <- 100
gene_num_null <- 200
gene_num_per_topic <- 100

##############

gene_library_repeating_vec <- c(0.8, rep(1, 8), 1.2)
gene_nuisance_values<- c(0.1, 1, 100)
gene_covariate_coefficient_proportion <- c(0.1, 0.1, 0.7, 0.1, 0.1, 0.1)
gene_covariate_coefficient_proportion <- gene_covariate_coefficient_proportion/sum(gene_covariate_coefficient_proportion)
gene_covariate_coefficient_size <- c(-1,-0.5,0,0.5,1,3)

num_topics <- 4
simplex <- matrix(0, nrow = k, ncol = num_topics)
for(j in 1:num_topics){
  idx <- sample(1:k, size = 3)
  simplex[idx,j] <- runif(length(idx))
}
for(i in 1:k){
  idx <- i %% num_topics + 1
  simplex[i,idx] <- simplex[i,idx] + runif(length(idx))
}
gene_topic_simplex <- simplex
gene_topic_latent_gaussian_noise <- 0.1*diag(k)
gene_topic_casecontrol_name <- c("none", "weak-negative", "strong-negative", "weak-positive", "strong-positive")
gene_topic_casecontrol_size <- c(0, -0.5, -5, 0.5,  2)
gene_topic_casecontrol_proportion <- c(0.74, 0.1, 0.03, 0.1, 0.03)
mat <- matrix(0, nrow = length(gene_nuisance_values), ncol = length(gene_topic_casecontrol_size))
mat[,1] <- c(0.8,0.15,0.05)
mat[,2] <- c(0.6,0.3,0.1)
mat[,3] <- c(0.1,0.2,0.7)
mat[,4] <- c(0.6,0.3,0.1)
mat[,5] <- c(0.1,0.2,0.7)
colnames(mat) <- paste0("cc_size:", gene_topic_casecontrol_name)
rownames(mat) <- paste0("nuisance:", gene_nuisance_values)
gene_nuisance_proporition_mat <- mat

##############

num_indiv <- 20
case_control_vec <- rep(c(0,1), each = num_indiv/2)
age_vec <- round(rnorm(num_indiv),2)
gender_vec <- rep(c(0,1), times = num_indiv/2)
tobacco_vec <- c(sample(c(0,1), size = num_indiv/2, replace = T, prob = c(0.7,0.3)),
                 sample(c(0,1), size = num_indiv/2, replace = T, prob = c(0.3,0.7)))
df <- data.frame(cc = as.numeric(case_control_vec),
                 age = as.numeric(age_vec),
                 gender = as.numeric(gender_vec),
                 tobacco = as.numeric(tobacco_vec))
individual_covariates = df
individual_case_control_variable = "cc"
individual_num_cells = 250

natural_param_max_quant <- 5
sparsity_downsampling <- 0


