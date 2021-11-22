rm(list=ls())
load("../../../../out/writeup8e/writeup8e_sns_layer23_esvd_extended.RData")
load("../../../../data/sns_autism/sns_formatted.RData")

metadata <- sns@meta.data
cell_names <- rownames(esvd_res_nb2$x_mat)
cell_idx <- sapply(cell_names, function(name){
  which(rownames(metadata) == name)
})
metadata <- metadata[cell_idx,]

save(mat, esvd_res_nb2, metadata,
     file = "../../../../out/writeup8e/writeup8e_sns_layer23_esvd_extended_portable.RData")

###################
rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_extended_portable.RData")

nat_mat1 <- tcrossprod(esvd_res_nb2$x_mat, esvd_res_nb2$y_mat)
nat_mat2 <- tcrossprod(esvd_res_nb2$covariates, esvd_res_nb2$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)
quantile(angle_vec[gene_idx])

quantile(esvd_res_nb2$b_mat[,"Log_UMI"])
length(which(esvd_res_nb2$b_mat[,"Log_UMI"] < 0))

#########################

nat_mat1 <- tcrossprod(esvd_res_nb2$x_mat, esvd_res_nb2$y_mat)
nat_mat2 <- tcrossprod(esvd_res_nb2$covariates[,"Intercept",drop=F], esvd_res_nb2$b_mat[,"Intercept",drop=F])
lambda_mat <- exp(nat_mat1+nat_mat2)
library_vec <- matrixStats::rowSums2(mat)

AplusR <- sweep(mat, MARGIN = 2, STATS = esvd_res_nb2$nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(lambda_mat, MARGIN = 2, STATS = esvd_res_nb2$nuisance_param_vec, FUN = "/")
RoverMuplusS <- sweep(RoverMu, MARGIN = 1, STATS = library_vec, FUN = "+")
posterior_mean_mat <- AplusR/RoverMuplusS
posterior_var_mat <- AplusR/RoverMuplusS^2

posterior_mean_mat <- posterior_mean_mat*1e5
round(quantile(posterior_mean_mat),2)

###############################

categorical_var <- c("individual", "diagnosis", "region", "sex", "Seqbatch") #, "individual")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
n <- nrow(mat)
covariates <- as.matrix(metadata[,numerical_var])
covariates <- cbind(1, covariates)
colnames(covariates)[1] <- "Intercept"

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
covariates <- covariates[,-which(colnames(covariates) == "diagnosis_ASD")]

#############################

gene_idx <- which(colnames(mat) == "TTF2")
# plot(lambda_mat[,gene_idx], mat[,gene_idx], pch = 16)
col_vec <- rep(1, nrow(lambda_mat))
col_vec[which(metadata$diagnosis == "ASD")] <- 2
par(mfrow = c(1,2))
plot(posterior_mean_mat[,gene_idx]*1e5, mat[,gene_idx]/library_vec*1e5,
     pch = 16, col = col_vec)
vioplot::vioplot(posterior_mean_mat[which(metadata$diagnosis == "ASD"),gene_idx],
                 posterior_mean_mat[which(metadata$diagnosis == "Control"),gene_idx])

mu <- lambda_mat[,gene_idx]
r <- esvd_res_nb2$nuisance_param_vec[gene_idx]
A <- mat[,gene_idx]
s <- library_vec
zz <- (r+A)/(r/mu+1)

#########################
