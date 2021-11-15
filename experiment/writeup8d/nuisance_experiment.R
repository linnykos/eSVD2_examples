rm(list=ls())

set.seed(10)
n <- 1000
mu <- 10
x_vec <- 0:50
size_vec <- c(0.1, 1, 10, 100)
y_mat <- sapply(size_vec, function(size_val){
  stats::dnbinom(x_vec, size = size_val, mu = mu)
})

y_max <- max(y_mat)
plot(NA, xlim = range(x_vec), ylim = c(0, y_max))
for(i in 1:ncol(y_mat)){
  lines(x_vec, y_mat[,i], col = i, lwd = 2)
}

y_max <- max(y_mat[-1,])
max_idx <- 90
plot(NA, xlim = range(x_vec[2:max_idx]), ylim = c(0, y_max))
for(i in 1:ncol(y_mat)){
  lines(x_vec[2:max_idx], y_mat[2:max_idx,i], col = i, lwd = 2)
}
