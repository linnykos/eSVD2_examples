rm(list=ls())
load("../../git/eSVD2_examples/error_reports/2021-12-05_nuisance.RData")

dim(mat); mat[1:5,1:5]
dim(mean_mat); mean_mat[1:5,1:5]

zero_prop <- apply(mat, 2, function(x){length(which(x == 0))/length(x)})
quantile(zero_prop)

idx <- which.min(zero_prop)
plot(mean_mat[,idx], mat[,idx], asp = T, pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.1),
     xlab = "Predicted value", ylab = "Observed value",
     main = paste0("Dense gene: 0-percentage (", round(zero_prop[idx], 3), ")"))
lines(c(0, 1e5), c(0, 1e5), col = 2, lwd = 2, lty = 2)
# dense genes track the y=x line

idx <- which.max(zero_prop)
plot(mean_mat[,idx], mat[,idx], asp = T, pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.1),
     xlab = "Predicted value", ylab = "Observed value",
     main = paste0("Sparse gene: 0-percentage (", round(zero_prop[idx], 3), ")"))
lines(c(0, 1e5), c(0, 1e5), col = 2, lwd = 2, lty = 2)
# dense genes have downward bias in non-zero observed values

nuisance_param_vec <- lapply(1:ncol(mat), function(j){
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
  val_final <- vec[which.min(abs(target_prob_vec - obs_prob))]

  estimates <- c(val1, val2, val3)
  names(estimates) <- c("theta.ml", "theta.mm", "glamGamPoi.mle")
  return(list(estimates = estimates, value = val_final))
})

range_vec <- sapply(nuisance_param_vec, function(x){
  diff(range(x$estimates))
})
quantile(range_vec)
idx <- which.max(range_vec)
plot(mean_mat[,idx], mat[,idx], asp = T, pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.1),
     xlab = "Predicted value", ylab = "Observed value",
     main = paste0("0-percentage (", round(zero_prop[idx], 3), ")\n",
                   paste0(round(nuisance_param_vec[[idx]]$estimates, 2), collapse = ", ")))
lines(c(0, 1e5), c(0, 1e5), col = 2, lwd = 2, lty = 2)

plot(zero_prop, pmin(range_vec, 1e100),
     pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.1),
     xlab = "0-percentage", ylab = "Range of different nuisance parameters")
