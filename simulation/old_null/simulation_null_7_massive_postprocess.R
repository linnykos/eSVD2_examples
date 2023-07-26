rm(list=ls())

# load in data

result_list <- result_list[which(sapply(result_list, length) >= 1)]

# populate the gaussian for eSVD
num_result <- length(result_list)
p <- length(result_list[[1]]$gaussian_teststat)
teststat_mat <- matrix(NA, nrow = num_result, ncol = p)
for(i in 1:num_result){
  teststat_mat[i,] <- result_list[[i]]$gaussian_teststat
}
cor_mat <- stats::cor(teststat_mat)
round(cor_mat[1:10,11:20],1)

################################
n_each <- 100; s <- 20; p <- 1000; k <- 2
n <- n_each*s
df <- cbind(1,
            0,
            rep(c(0,1), each = s/2),
            rep(c(0,1), times = s/2),
            scale(round(rnorm(s, mean = 30, sd = 5)), center = F, scale = T),
            1:s)
colnames(df) <- c("Intercept", "Log_UMI", "CC", "Sex", "Age", "Individual")
# expand to covariate matrix
covariate <- do.call(rbind, lapply(1:nrow(df), function(i){
  matrix(rep(df[i,], each = n_each), nrow = n_each, ncol = ncol(df))
}))
colnames(covariate) <- colnames(df)[1:ncol(df)]

.compute_denoised_teststat <- function(covariate, df, eSVD_obj){
  denoised_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat,
                             eSVD_obj$fit_Second$y_mat) +
    tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
  denoised_mat <- exp(denoised_mat)

  cell_individual_list <- lapply(1:max(covariate[,"Individual"]), function(i){which(covariate[,"Individual"] == i)})
  mean_denoised <- t(sapply(cell_individual_list, function(idx){
    Matrix::colMeans(denoised_mat[idx,])
  }))
  var_denoised <- t(sapply(cell_individual_list, function(idx){
    matrixStats::colVars(denoised_mat[idx,])
  }))

  case_df_idx <- which(df[,"CC"] == 1)
  control_df_idx <- which(df[,"CC"] == 0)
  n1 <- length(case_df_idx)
  n2 <- length(control_df_idx)
  case_gaussian_mean <- Matrix::colMeans(mean_denoised[case_df_idx,,drop = F])
  control_gaussian_mean <- Matrix::colMeans(mean_denoised[control_df_idx,,drop = F])

  case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = mean_denoised[case_df_idx,,drop = F],
    avg_posterior_var_mat = var_denoised[case_df_idx,,drop = F]
  )
  control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = mean_denoised[control_df_idx,,drop = F],
    avg_posterior_var_mat = var_denoised[control_df_idx,,drop = F]
  )
  teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
    (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))
  plot(teststat_vec)

  numerator_vec <- (case_gaussian_var/n1 + control_gaussian_var/n2)^2
  denominator_vec <- (case_gaussian_var/n1)^2/(n1-1) + (control_gaussian_var/n2)^2/(n2-1)
  df_vec <- numerator_vec/denominator_vec
  p <- ncol(mean_denoised)
  gaussian_teststat <- sapply(1:p, function(j){
    qnorm(pt(teststat_vec[j], df = df_vec[j]))
  })

  gaussian_teststat
}

teststat_mat2 <- matrix(NA, nrow = num_result, ncol = p)
for(i in 1:num_result){
  print(i)
  teststat_mat2[i,] <- .compute_denoised_teststat(
    covariate = covariate,
    df = df,
    eSVD_obj = result_list[[i]]$eSVD_obj
  )
}

round(teststat_mat2[1:10,1:20],1)
cor_mat2 <- stats::cor(teststat_mat2)
round(cor_mat2[1:10,11:20],1)

diag(cor_mat2) <- NA
cor_mat2[1:10,1:10] <- NA
quantile(cor_mat2, probs = seq(0,1,length.out=11), na.rm = T)

save(teststat_mat, teststat_mat2,
     file = "tmp2.RData")


#######################
mean_vec <- sapply(1:nrow(teststat_mat), function(i){
  mean(teststat_mat[i,-c(1:10)])
})

mean_vec2 <- sapply(1:nrow(teststat_mat2), function(i){
  mean(teststat_mat2[i,-c(1:10)])
})
teststat_mat2 <- pmin(teststat_mat2, 10)

p <- 1000
y_block_assignment <- rep(c(1:3), times = ceiling(p/3))[1:p]
signal_idx <- 1:10
teststat_mat_condense <- matrix(NA, ncol = 6, nrow = nrow(teststat_mat))
teststat2_mat_condense <- matrix(NA, ncol = 6, nrow = nrow(teststat_mat2))

for(j in 1:3){
  teststat_mat_condense[,j] <- rowMeans(teststat_mat[,intersect(which(y_block_assignment == j),
                                                               signal_idx)])
}
for(j in 1:3){
  teststat_mat_condense[,j+3] <- rowMeans(teststat_mat[,intersect(which(y_block_assignment == j),
                                                                c(1:p)[-signal_idx])])
}

zz <- cor(teststat_mat_condense, method = "kendall"); zz[1,4]; zz[2,5]; zz[3,6]

for(j in 1:3){
  teststat2_mat_condense[,j] <- rowMeans(teststat_mat2[,intersect(which(y_block_assignment == j),
                                                                signal_idx)])
}
for(j in 1:3){
  teststat2_mat_condense[,j+3] <- rowMeans(teststat_mat2[,intersect(which(y_block_assignment == j),
                                                                  c(1:p)[-signal_idx])])
}
zz <- cor(teststat2_mat_condense, method = "kendall"); zz[1,4]; zz[2,5]; zz[3,6]



