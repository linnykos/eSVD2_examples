rm(list=ls())
set.seed(10)
library(Seurat)
source("../eSVD2_examples/simulation/data_generator.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

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
gene_null_casecontrol_proportion <- c(0.8, 0.15, 0.05)
gene_null_casecontrol_size <- c(0, -0.75, 0.75)
gene_null_latent_gaussian_noise <- 0.5*diag(k)

gene_num_mixed_membership <- 100
gene_num_null <- 200
gene_num_per_topic <- 100

##############

gene_library_repeating_vec <- c(0.8, rep(1, 8), 1.2)
gene_nuisance_values<- c(0.1, 1, 100)

set.seed(10)
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
gene_topic_casecontrol_size <- c(0, -0.5, -0.75, 0.5, .75)
gene_topic_casecontrol_proportion <- c(0.6, 0.15, 0.05, 0.15, 0.05)
mat <- matrix(0, nrow = length(gene_nuisance_values), ncol = length(gene_topic_casecontrol_size))
mat[,1] <- c(0.8,0.15,0.05)
mat[,2] <- c(0.6,0.35,0.05)
mat[,3] <- c(0.1,0.2,0.7)
mat[,4] <- c(0.6,0.35,0.05)
mat[,5] <- c(0.1,0.2,0.7)
colnames(mat) <- paste0("cc_size:", gene_topic_casecontrol_name)
rownames(mat) <- paste0("nuisance:", gene_nuisance_values)
gene_nuisance_proporition_mat <- mat

gene_covariate_coefficient_size <- c(-2,-0.5,0,0.5,1,2)
mat <- matrix(0, nrow = length(gene_covariate_coefficient_size),
              ncol = length(gene_topic_casecontrol_name))
colnames(mat) <- paste0("cc_size:", gene_topic_casecontrol_name)
rownames(mat) <- paste0("coef_size:", gene_nuisance_values)
mat[,1] <- c(0.1, 0.1, 0.6, 0.1, 0.1,   0)
mat[,2] <- c(0.6, 0.4,   0,   0,   0,   0)
mat[,3] <- c(  0, 0.2, 0.8,   0,   0,   0)
mat[,4] <- c(  0,   0,   0,   0, 0.4, 0.6)
mat[,5] <- c(  0,   0, 0.8, 0.2,   0,   0)
gene_covariate_coefficient_proportion_mat <- mat

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

#########################

set.seed(10)
input_obj <- data_generator_nat_mat(
  cell_latent_gaussian_mean = cell_latent_gaussian_mean,
  cell_latent_gaussian_covariance = cell_latent_gaussian_covariance,
  gene_library_repeating_vec = gene_library_repeating_vec,
  gene_covariate_coefficient_proportion_mat = gene_covariate_coefficient_proportion_mat,
  gene_covariate_coefficient_size = gene_covariate_coefficient_size,
  gene_nuisance_values = gene_nuisance_values,
  gene_nuisance_proporition_mat = gene_nuisance_proporition_mat,
  gene_null_casecontrol_name = gene_null_casecontrol_name,
  gene_null_casecontrol_proportion = gene_null_casecontrol_proportion,
  gene_null_casecontrol_size = gene_null_casecontrol_size,
  gene_null_latent_gaussian_noise = gene_null_latent_gaussian_noise,
  gene_null_nuisance_proportion = gene_null_nuisance_proportion,
  gene_num_mixed_membership = gene_num_mixed_membership,
  gene_num_null = gene_num_null,
  gene_num_per_topic = gene_num_per_topic,
  gene_topic_latent_gaussian_noise = gene_topic_latent_gaussian_noise,
  gene_topic_simplex = gene_topic_simplex,
  gene_topic_casecontrol_name = gene_topic_casecontrol_name,
  gene_topic_casecontrol_proportion = gene_topic_casecontrol_proportion,
  gene_topic_casecontrol_size = gene_topic_casecontrol_size,
  individual_covariates = individual_covariates,
  individual_case_control_variable = individual_case_control_variable,
  individual_num_cells = individual_num_cells
)

input_obj <- data_signal_enhancer(input_obj,
                                  global_shift = 0)
nat_mat <- input_obj$nat_mat
input_obj <- data_generator_obs_mat(input_obj)

#####################################3

case_individuals = input_obj$case_individuals
control_individuals = input_obj$control_individuals
covariates = input_obj$covariates
gene_labeling = input_obj$gene_labeling
gene_labeling2 = input_obj$gene_labeling2
gene_library_vec = input_obj$gene_library_vec
individual_vec = input_obj$individual_vec
nat_mat = input_obj$nat_mat
nuisance_vec = input_obj$nuisance_vec
obs_mat = input_obj$obs_mat
x_mat = input_obj$x_mat
y_mat = input_obj$y_mat
z_mat = input_obj$z_mat

# mean_mat <- exp(nat_mat)
# case_idx <- which(covariates[,"cc"] == 1)
# control_idx <- which(covariates[,"cc"] == 0)
# case_mean <- colMeans(mean_mat[case_idx,])
# control_mean <- colMeans(mean_mat[control_idx,])
# diff_mean <- case_mean - control_mean
# var_mat <- sweep(x = mean_mat, MARGIN = 2, STATS = nuisance_vec, FUN = "/")
#
# res <- eSVD2:::compute_test_statistic.default(
#   input_obj = mean_mat,
#   posterior_var_mat = var_mat,
#   case_individuals = case_individuals,
#   control_individuals = control_individuals,
#   individual_vec = individual_vec
# )
# teststat_vec <- res$teststat_vec
# true_teststat_vec <- teststat_vec

true_teststat_vec <- apply(nat_mat, 2, function(x){
  df <- as.data.frame(cbind(y = x, covariates[,1:6]))
  colnames(df)[1] <- "y"
  lm_res <- stats::lm(y ~ . - 1, data = df)
  stats::coef(lm_res)["cc"]
})

col_palette <- c("none" = rgb(0.5, 0.5, 0.5),
                 "strong-negative" = rgb(0.75, 0, 0),
                 "strong-positive" = rgb(0, 0.75, 0),
                 "weak-negative" = rgb(1, 0.5, 0.9),
                 "weak-positive" = rgb(0.5, 1, 0.9))
col_vec <- plyr::mapvalues(gene_labeling2, from = names(col_palette), to = col_palette)
plot(true_teststat_vec, col = col_vec, pch = 16)
hist(true_teststat_vec, breaks = 50)


quantile(apply(obs_mat, 2, function(x){length(which(x==0))/length(x)}))
table(gene_labeling2, z_mat[,"age"])

gene_casecontrol_name <- sort(unique(gene_topic_casecontrol_name, gene_null_casecontrol_name))
for(x in gene_casecontrol_name){
  idx <- which(gene_labeling2 == x)
  print(x)
  print(round(quantile(teststat_vec[idx], probs = seq(0,1,length.out=11)), 2))
  print("====")
}

####################

tmp <- data.frame(covariates[,2:5])
tmp$individual <- individual_vec
seurat_obj <- Seurat::CreateSeuratObject(counts = t(obs_mat),
                                         meta.data = tmp)
seurat_obj[["RNA"]]@var.features <- rownames(seurat_obj)
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = F)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:15)

save(seurat_obj,
     case_individuals,
     control_individuals,
     covariates,
     gene_labeling,
     gene_labeling2,
     gene_library_vec,
     individual_vec,
     nat_mat,
     nuisance_vec,
     obs_mat,
     true_teststat_vec,
     x_mat,
     y_mat,
     z_mat,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_2.RData")

Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "individual")
Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "cc")
Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "gender")
Seurat::FeaturePlot(seurat_obj, reduction = "umap", features = "age")

obs_mat <- t(as.matrix(seurat_obj[["RNA"]]@counts))
zero_prop <- apply(obs_mat, 2, function(x){
  length(which(x == 0))/length(x)
})
quantile(zero_prop)
