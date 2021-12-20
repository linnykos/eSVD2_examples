rm(list=ls())
library(Seurat)

load("../../../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson4.RData")
mat[which(mat == min(mat))] <- 0

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
nuisance_param_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  val1 <- MASS::theta.ml(y = mat[,j], mu = mean_mat[,j])
  val2 <- MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
  val3 <- glmGamPoi::overdispersion_mle(y = mat[,j], mean = mean_mat[,j])$estimate

  vec <- c(val1, val2, val3)
  vec <- vec[!is.na(vec)]
  if(length(vec) == 1) return(vec[1])
  vec <- pmax(pmin(vec, 1e5), 0.1)
  vec <- c(vec, c(0.1, 0.5, 1))

  obs_prob <- length(which(mat[,j] == 0))/nrow(mat)
  target_prob_vec <- sapply(vec, function(val){
    mean((1+mean_mat[,j]/val)^(-val))
  })
  return(vec[which.min(abs(target_prob_vec - obs_prob))])
})

########

true_esvd <- esvd_res_full
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("true_esvd", "nuisance_param_vec", "metadata")]
rm(list = ls_vec)

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

x_row_max <- matrixStats::rowMaxs(abs(true_esvd$x_mat))
idx <- which(x_row_max >= 2.5)
true_esvd$x_mat[idx,] <- true_esvd$x_mat[idx,]/2

y_row_max <- matrixStats::rowMaxs(abs(true_esvd$y_mat))
idx <- which(y_row_max >= 1)
true_esvd$y_mat[idx,] <- true_esvd$y_mat[idx,]/2

b_row_max <- matrixStats::rowMaxs(abs(true_esvd$b_mat))
idx <- which(b_row_max >= 3)
true_esvd$b_mat[idx,] <- true_esvd$b_mat[idx,]/5

true_esvd$b_mat[,c(3:6)] <- true_esvd$b_mat[,c(3:6)]/2
true_esvd$b_mat[,"age"] <- true_esvd$b_mat[,"age"]/10
true_esvd$b_mat[,"RNA.Integrity.Number"] <- true_esvd$b_mat[,"RNA.Integrity.Number"]/2
true_esvd$b_mat[,"post.mortem.hours"] <- true_esvd$b_mat[,"post.mortem.hours"]/10
true_esvd$b_mat[,"nFeature_RNA"] <- true_esvd$b_mat[,"nFeature_RNA"]/1000

################

nat_mat1 <- tcrossprod(true_esvd$x_mat, true_esvd$y_mat)
nat_mat2 <- tcrossprod(true_esvd$covariates[,-library_idx], true_esvd$b_mat[,-library_idx])
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
mean_mat[mean_mat >= 100] <- 100

# add signal to the autistic genes
set.seed(10)
autism_gene_idx <- sample(1:p, size = round(p/50))
multiplier_vec <- rep(1, p)
multiplier_vec[autism_gene_idx] <- rnorm(length(autism_gene_idx), mean = 5, sd = 1)
autism_cell_idx <- which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5)
control_cell_idx <- which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5)
for(j in 1:p){
  vec <- nat_mat[,j]
  vec_autism <- vec[autism_cell_idx]
  vec_control <- vec[control_cell_idx]

  if(j %in% autism_gene_idx){
    target <- multiplier_vec[j] * max(mean(exp(vec_control)), mean(exp(vec_autism)))
  } else {
    target <- multiplier_vec[j] * mean(exp(vec_control))
  }
  current <- mean(exp(vec_autism))
  value <- log(target/current)
  true_esvd$b_mat[j,autism_idx] <- value
}
quantile(true_esvd$b_mat[autism_gene_idx, autism_idx])
quantile(true_esvd$b_mat[-autism_gene_idx, autism_idx])

nat_mat1 <- tcrossprod(true_esvd$x_mat, true_esvd$y_mat)
nat_mat2 <- tcrossprod(true_esvd$covariates[,-library_idx], true_esvd$b_mat[,-library_idx])
nat_mat <- nat_mat1 + nat_mat2
nat_mat[nat_mat >= log(100)] <- log(100)
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx])
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx])

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

# scale up lambda according to the library size
library_mat <- sapply(1:p, function(j){
  exp(true_esvd$covariates[,"Log_UMI",drop = F]*true_esvd$b_mat[j,"Log_UMI"])
})
library_mat <- pmin(library_mat, 10)

mat <- lambda_mat
tol <- 1e-3
for(i in 1:n){
  set.seed(i)
  mat[i,] <- stats::rpois(p, lambda = library_mat[i,]*lambda_mat[i,] + tol)
}
true_esvd$covariates[,library_idx] <- log(matrixStats::rowSums2(mat))

length(which(mat == 0))/prod(dim(mat))
tmp <- mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx]
length(which(tmp == 0))/prod(dim(tmp))
quantile(tmp[tmp > 0])

tmp <- mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx]
length(which(tmp == 0))/prod(dim(tmp))
quantile(tmp[tmp > 0])


#########################
covariates <- true_esvd$covariates

# regress all variables against intercept + Log_UMI + diagnosis_ASD
keep_idx <- which(colnames(covariates) %in% c("Intercept", "Log_UMI", "diagnosis_ASD"))
other_idx <- c(1:ncol(covariates))[-keep_idx]

for(j in other_idx){
  df_tmp <- data.frame(covariates[,j], covariates[,keep_idx])
  colnames(df_tmp)[1] <- "tmp"
  lm_fit <- stats::lm("tmp ~ . ", data = df_tmp)
  vec_tmp <- stats::residuals(lm_fit)
  covariates[,j] <- vec_tmp
}


######################
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("true_esvd", "mat", "autism_gene_idx",
                                "nat_mat", "lambda_mat" , "gamma_mat",
                                "library_mat", "covariates",
                                "nuisance_param_vec", "metadata")]
rm(list = ls_vec)
save.image("../../../../out/writeup8e/writeup8e_sns_pseudoreal_data.RData")

