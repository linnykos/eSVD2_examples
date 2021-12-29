rm(list=ls())
library(Seurat)
load("../../../../out/writeup8f/writeup8f_sns_layer23_esvd.RData")

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

zero_prob <- apply(mat, 2, function(x){
  length(which(x == 0))/nrow(mat)
})
nuisance_param_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  val1 <- MASS::theta.ml(y = mat[,j], mu = mean_mat[,j])
  val2 <- MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
  val3 <- glmGamPoi::overdispersion_mle(y = mat[,j], mean = mean_mat[,j])$estimate

  vec <- c(val1, val2, val3)
  vec <- vec[!is.na(vec)]
  if(length(vec) == 1) return(vec[1])
  vec <- pmax(pmin(vec, 1e5), 1)
  vec <- c(vec, 1)

  obs_prob <- length(which(mat[,j] == 0))/nrow(mat)
  target_prob_vec <- sapply(vec, function(val){
    mean((1+mean_mat[,j]/val)^(-val))
  })
  return(vec[which.min(abs(target_prob_vec - obs_prob))])
})
quantile(nuisance_param_vec)

########

ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("esvd_res_full", "nuisance_param_vec", "metadata")]
rm(list = ls_vec)
true_esvd <- esvd_res_full

n <- nrow(true_esvd$x_mat)
p <- nrow(true_esvd$y_mat)
colnames(true_esvd$b_mat) <- colnames(true_esvd$covariates)
autism_idx <- which(colnames(true_esvd$covariates) == "diagnosis_ASD")
library_idx <- which(colnames(true_esvd$covariates) == "Log_UMI")

#################
# fix the simulation data
true_esvd$b_mat[,autism_idx] <- 0

y_max <- quantile(abs(true_esvd$y_mat), prob = 0.99)
true_esvd$y_mat[which(abs(true_esvd$y_mat) >= y_max)] <- y_max*sign(true_esvd$y_mat[which(abs(true_esvd$y_mat) >= y_max)])

x_max <- quantile(abs(true_esvd$x_mat), prob = 0.99)
true_esvd$x_mat[which(abs(true_esvd$x_mat) >= x_max)] <- x_max*sign(true_esvd$x_mat[which(abs(true_esvd$x_mat) >= x_max)])

b_max <- quantile(abs(true_esvd$b_mat), prob = 0.95)
true_esvd$b_mat[which(abs(true_esvd$b_mat) >= b_max)] <- b_max*sign(true_esvd$b_mat[which(abs(true_esvd$b_mat) >= b_max)])

################

nat_mat1 <- tcrossprod(true_esvd$x_mat, true_esvd$y_mat)
nat_mat2 <- tcrossprod(true_esvd$covariates[,-library_idx], true_esvd$b_mat[,-library_idx])
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

# add signal to the autistic genes
set.seed(10)
true_esvd$b_mat[,autism_idx] <- 0
autism_gene_idx <- sample(1:p, size = round(p/50))
multiplier_vec <- rep(1, p)
up_idx <- sample(autism_gene_idx, size = round(length(autism_gene_idx)/2))
down_idx <- setdiff(autism_gene_idx, up_idx)
multiplier_vec[up_idx] <- runif(n = length(up_idx), min = 1.1, max = 2.1)
multiplier_vec[down_idx] <- runif(n = length(down_idx), min = 1/2.1, max = 1/1.1)
autism_cell_idx <- which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5)
control_cell_idx <- which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5)
for(j in autism_gene_idx){
  vec <- nat_mat[,j]
  vec_autism <- vec[autism_cell_idx]
  vec_control <- vec[control_cell_idx]

  target <- multiplier_vec[j] * mean(exp(vec_control))
  current <- mean(exp(vec_autism))
  value <- log(target/current)
  true_esvd$b_mat[j,autism_idx] <- value
}
quantile(true_esvd$b_mat[up_idx, autism_idx])
quantile(true_esvd$b_mat[down_idx, autism_idx])

###########

# assess the relative effects
range_mat <- sapply(1:ncol(true_esvd$b_mat), function(j){
  print(j)
  quantile(tcrossprod(true_esvd$covariates[,j], true_esvd$b_mat[,j]))
})
colnames(range_mat) <- colnames(true_esvd$covariates)
round(range_mat,2)

indiv_columns <- grep("individual", colnames(true_esvd$b_mat))
indiv_vec <- unique(metadata$individual)[1:length(indiv_columns)]
diagnosis_vec <- sapply(indiv_vec, function(indiv_name){
  idx <- which(metadata$individual == indiv_name)
  unique(metadata$diagnosis[idx])
})
tmp_df <- data.frame(name = indiv_vec,
                     number = 1:length(indiv_columns),
                     diagnosis = diagnosis_vec)

# upweight and downweight certain genes
set.seed(10)
max_diff <- diff(range(range_mat))/2
tmp_gene_idx <- sample(1:nrow(true_esvd$b_mat), size = round(nrow(true_esvd$b_mat)/5))
for(j in tmp_gene_idx){
  diagnosis_val <- sample(c("Control", "ASD"), 1)
  number_vec <- tmp_df$number[which(tmp_df$diagnosis == diagnosis_val)]
  col_idx <- which(colnames(true_esvd$b_mat) %in% paste0("individual_", number_vec))
  true_esvd$b_mat[j,col_idx] <- runif(length(col_idx), min = 0, max = max_diff)
  true_esvd$b_mat[j,-col_idx] <- 0
}

###########

nat_mat1 <- tcrossprod(true_esvd$x_mat, true_esvd$y_mat)
nat_mat2 <- tcrossprod(true_esvd$covariates[,-library_idx], true_esvd$b_mat[,-library_idx])
nat_mat <- nat_mat1 + nat_mat2
nat_mat[nat_mat > log(500)] <- log(500)
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), up_idx])
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), up_idx])
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), down_idx])
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), down_idx])

#######################

# now simulate data
# first simulate the gamma
lambda_mat <- nat_mat
for(j in 1:ncol(lambda_mat)){
  set.seed(j)
  lambda_mat[,j] <- stats::rgamma(n,
                                  shape = nuisance_param_vec[j],
                                  scale = exp(nat_mat[,j])/nuisance_param_vec[j])
}
quantile(lambda_mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), up_idx])
quantile(lambda_mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), up_idx])
quantile(lambda_mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), down_idx])
quantile(lambda_mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), down_idx])

# scale up lambda according to the library size
library_mat <- sapply(1:p, function(j){
  exp(true_esvd$covariates[,"Log_UMI",drop = F]*true_esvd$b_mat[j,"Log_UMI"]/5)
})

mat <- lambda_mat
tol <- 1e-3
for(i in 1:n){
  set.seed(i)
  mat[i,] <- stats::rpois(p, lambda = library_mat[i,]*lambda_mat[i,] + tol)
}
length(which(mat == 0))/prod(dim(mat))

quantile(mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), up_idx])
quantile(mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), up_idx])
quantile(mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), down_idx])
quantile(mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), down_idx])

#########################
true_esvd$covariates[,library_idx] <- log(matrixStats::rowSums2(mat))

categorical_var <- c("individual", "diagnosis", "region", "sex", "Seqbatch") #, "individual")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
n <- nrow(mat)
covariates <- as.matrix(metadata[,numerical_var])
covariates <- cbind(1, log(matrixStats::rowSums2(mat)), covariates)
colnames(covariates)[1:2] <- c("Intercept", "Log_UMI")

for(variable in categorical_var){
  vec <- metadata[,variable]
  uniq_level <- unique(vec)
  for(i in uniq_level[-1]){
    tmp <- rep(0, n)
    tmp[which(vec == i)] <- 1

    var_name <- paste0(variable, "_", i)
    covariates <- cbind(covariates, tmp)
    colnames(covariates)[ncol(covariates)] <- var_name
  }
}

# regress all variables against intercept + Log_UMI + diagnosis_ASD
keep_idx <- which(colnames(covariates) %in% c("Intercept", "Log_UMI"))
other_idx <- which(colnames(covariates) %in% c("diagnosis_ASD", "age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt",
                                               "nFeature_RNA", "region_PFC", "sex_F", "Seqbatch_SB2", "Seqbatch_SB1"))
for(j in other_idx){
  df_tmp <- data.frame(covariates[,j], covariates[,keep_idx])
  colnames(df_tmp)[1] <- "tmp"
  lm_fit <- stats::lm("tmp ~ . ", data = df_tmp)
  vec_tmp <- stats::residuals(lm_fit)
  covariates[,j] <- vec_tmp
}

cols_regress_out <- grep("individual", colnames(covariates))
covariates_new <- covariates[,-cols_regress_out,drop = F]
for(i in 1:length(cols_regress_out)){
  df_tmp <- data.frame(covariates[,cols_regress_out[i]], covariates_new)
  colnames(df_tmp)[1] <- "tmp"
  lm_fit <- stats::lm("tmp ~ . - 1", data = df_tmp)
  vec_tmp <- stats::residuals(lm_fit)
  if(sum(abs(vec_tmp)) < 1e-6) break()
  covariates_new <- cbind(covariates_new, vec_tmp)
  colnames(covariates_new)[ncol(covariates_new)] <- paste0("individual_",i)
}
covariates <- covariates_new
stopifnot(all(dim(covariates) == dim(true_esvd$covariates)))

######################
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("true_esvd", "mat", "autism_gene_idx",
                                "up_idx", "down_idx",
                                "nat_mat", "lambda_mat" , "gamma_mat",
                                "library_mat", "covariates",
                                "nuisance_param_vec", "metadata")]
rm(list = ls_vec)
save.image("../../../../out/writeup8f/writeup8f_pseudoreal_data.RData")

######################

quantile(apply(mat, 2, max))
case_idx <- which(metadata[,"diagnosis"] == "ASD")
control_idx <- which(metadata[,"diagnosis"] == "Control")

x_vec <- sapply(1:ncol(mat), function(j){
  log2(mean(mat[case_idx,j])) - log2(mean(mat[control_idx,j]))
})
quantile(x_vec)
quantile(x_vec[autism_gene_idx])
quantile(x_vec[-autism_gene_idx])
