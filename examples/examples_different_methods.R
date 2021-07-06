rm(list=ls())

# from https://cran.r-project.org/web/packages/glmpca/vignettes/glmpca.html
set.seed(202)
ngenes <- 5000 #must be divisible by 10
ngenes_informative <- ngenes*.1
ncells <- 50 #number of cells per cluster, must be divisible by 2
nclust <- 3
# simulate two batches with different depths
batch <- rep(1:2, each = nclust*ncells/2)
ncounts <- stats::rpois(ncells*nclust, lambda = 1000*batch)
# generate profiles for 3 clusters
profiles_informative <- replicate(nclust, exp(stats::rnorm(ngenes_informative)))
profiles_const <- matrix(ncol=nclust, rep(exp(stats::rnorm(ngenes-ngenes_informative)), nclust))
profiles <- rbind(profiles_informative,profiles_const)
# generate cluster labels
clust <- sample(rep(1:3, each = ncells))
# generate single-cell transcriptomes
counts <- sapply(seq_along(clust), function(i){
  stats::rmultinom(1, ncounts[i], prob = profiles[,clust[i]])
})
rownames(counts) <- paste("gene", seq(nrow(counts)), sep = "_")
colnames(counts) <- paste("cell", seq(ncol(counts)), sep = "_")
# clean up rows
dat <- t(counts[rowSums(counts) > 0, ])

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

########################

# try just PCA
pca_res <- stats::prcomp(dat)
plot(pca_res$x[,1], pca_res$x[,2], asp = T, col = clust, pch = 16)

# try UMAP
dat2 <- scale(dat, center = T, scale = T)
svd_res <- RSpectra::svds(dat2, k = 10)
dim_red <- .mult_mat_vec(svd_res$u, svd_res$d)
set.seed(10)
umap_res <- Seurat::RunUMAP(dim_red)@cell.embeddings
plot(umap_res[,1], umap_res[,2], asp = T, col = clust, pch = 16)

# try glm-pca
set.seed(10)
glm_res <- glmpca::glmpca(t(dat), 2, fam="poi")
plot(glm_res$factors[,1], glm_res$factors[,2], asp = T, col = clust, pch = 16)
plot(glm_res$loadings[,1], glm_res$loadings[,2], asp = T, pch = 16)
plot(glm_res$dev)
plot(glm_res$coefX[,1])
image(.rotate(glm_res$factors[order(clust),]))

# try eSVD: attempt 1
set.seed(10)
library_size_vec <- rowSums(dat)
n <- nrow(dat)
covariates <- matrix(1, nrow = n, ncol = 1)
init <- eSVD2::initialize_esvd(dat, k = 2, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = NA,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
# image(.rotate(init$x_mat[order(clust),]))
esvd_res <- opt_esvd(init$x_mat, init$y_mat, dat, family = "poisson",
                     nuisance_param_vec = NA, library_size_vec = NA,
                     b_init = init$b_mat, covariates = covariates,
                     max_iter = 50, verbose = 1)
plot(esvd_res$x_mat[,1], esvd_res$x_mat[,2], asp = T, col = clust, pch = 16)


# try eSVD: attempt 2
set.seed(10)
n <- nrow(dat)
covariates <- cbind(matrix(1, nrow = n, ncol = 1), log(library_size_vec))
init <- eSVD2::initialize_esvd(dat, k = 2, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = 1,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
# image(.rotate(init$x_mat[order(clust),]))
esvd_res <- opt_esvd(init$x_mat, init$y_mat, dat, family = "poisson",
                nuisance_param_vec = NA, library_size_vec = 1,
                b_init = init$b_mat, covariates = covariates,
                max_iter = 50, verbose = 1)
plot(esvd_res$x_mat[,1], esvd_res$x_mat[,2], asp = T, col = clust, pch = 16)

nat_mat <- x_mat %*% t(y_mat) + covariates %*% t(esvd_res$b_mat)
quantile(exp(nat_mat))
mean_mat <- .mult_vec_mat(library_size_vec, exp(nat_mat))
image(.rotate(nat_mat))
image(.rotate(esvd_res$x_mat[order(clust),]))
plot(esvd_res$loss)
plot(x_mat[,1], x_mat[,2], asp = T, col = clust, pch = 16)
plot(y_mat[,1], y_mat[,2], asp = T, pch = 16)
set.seed(10)
umap_res <- Seurat::RunUMAP(esvd_res$x)@cell.embeddings
plot(umap_res[,1], umap_res[,2], asp = T, col = clust, pch = 16)


round(crossprod(x_mat), 2)*sqrt(p/n)
round(crossprod(y_mat), 2)*sqrt(n/p)

