rm(list=ls())
library(eSVD2)

set.seed(123)
n_each <- 100
num_cluster <- 3
n <- n_each*num_cluster
p <- 150
k <- 5
mu_vec <- c(-5, 0, 5)
x_mat <- do.call(rbind, lapply(1:num_cluster, function(i){
  MASS::mvrnorm(n = n_each, mu = rep(mu_vec[i], k), Sigma = diag(k)*.01)
}))
y_mat <- abs(MASS::mvrnorm(n = p,
                           mu = rep(0, k),
                           Sigma = 0.5*diag(k)+matrix(0.5, nrow = k, ncol = k))*.1)

# visualize initial data
set.seed(10)
umap_res <- Seurat::RunUMAP(x_mat, metric = "euclidean")
umap_res <- umap_res@cell.embeddings
plot(umap_res[,1], umap_res[,2], pch = 16,
     col = rep(1:3, each = n_each), asp = T,
     main = "True embedding")

# Simulate data
canonical_mat <- tcrossprod(x_mat, y_mat)
prob_mat <- matrix(NA, nrow = n, ncol = p)
dat <- matrix(NA, nrow = n, ncol = p)
set.seed(10)
for(i in 1:n){
  for(j in 1:p){
    canonical <- canonical_mat[i,j]
    prob_mat[i,j] <- exp(canonical)/(1+exp(canonical))
    dat[i,j] <- stats::rbinom(n = 1, size = 1, prob = prob_mat[i,j])
  }
}
image(prob_mat, main = "Probability matrix")
image(dat, main = "Observed binary matrix")
dat <- Matrix::Matrix(dat, sparse = T)
rownames(dat) <- paste0("c", 1:n)
colnames(dat) <- paste0("g", 1:p)

########
## The tricky part of using eSVD is actually finding a good way to initialize the factorization.
## Below is probably the "smarter" way to initialize, but this does too good of a job that
## it makes it hard to appreciate how well the alternating-minimization works.
# tol <- 1e-3
# svd_res <- irlba::irlba(dat, nv = 50)
# smoothed_mat <- tcrossprod(svd_res$u %*% diag(svd_res$d), svd_res$v)
# smoothed_mat <- pmax(pmin(smoothed_mat, 1-tol), tol)
# canonical_mat <- log(smoothed_mat / (1-smoothed_mat))
# svd_res2 <- irlba::irlba(canonical_mat, nv = 5)
# x_init <- svd_res2$u %*% diag(sqrt(svd_res2$d))
# y_init <- svd_res2$v %*% diag(sqrt(svd_res2$d))

## We do a random initialization just for this demonstration. (In practice, we wouldn't want to use a
## random initialization strategy.)
set.seed(10)
x_init <- MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = diag(k))
y_init <- MASS::mvrnorm(n = p, mu = rep(0, k), Sigma = diag(k))

# visualize results after initialization
set.seed(10)
umap_res <- Seurat::RunUMAP(x_init, metric = "euclidean")
umap_res <- umap_res@cell.embeddings
plot(umap_res[,1], umap_res[,2], pch = 16,
     col = rep(1:3, each = n_each), asp = T,
     main = "Initialization", )

# fit eSVD
eSVD_obj <- eSVD2::opt_esvd(input_obj = dat,
                            x_init = x_init,
                            y_init = y_init,
                            family = "bernoulli",
                            l2pen = 0.01)
plot(eSVD_obj$loss)

# visualize results after eSVD
set.seed(10)
umap_res <- Seurat::RunUMAP(eSVD_obj$x_mat, metric = "euclidean")
umap_res <- umap_res@cell.embeddings
plot(umap_res[,1], umap_res[,2], pch = 16,
     col = rep(1:3, each = n_each), asp = T,
     main = "After eSVD")

# plot results
estimated_canonical <- tcrossprod(eSVD_obj$x_mat, eSVD_obj$y_mat)
estimated_probability <- exp(estimated_canonical)/(1+exp(estimated_canonical))
image(estimated_probability, main = "Estimated probability")
plot(x = as.numeric(prob_mat),
     y = as.numeric(estimated_probability),
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1),
     asp = T,
     xlab = "True probability",
     ylab = "Estimated probability")

