rm(list=ls())
set.seed(10)

source("../eSVD2_examples/simulation/simulator_generator.R")
source("../eSVD2_examples/simulation/simulator_methods.R") # you don't need to load this is you want to only want to generate data

# I've found the following parameters to work well (to yield matrices that contain values between 0 and 1000)
# - curved_gaussian: modifiter = 1/250, nuisance_param = 2
# - neg_binom: modifier = 1/250, nuisance_param = 50
# - pcmf: modifier = 1/2, dropout_prob = 0.1
# - poisson: modifier = 1/300
# - zinbwave: modifier = 1/250, nuisance_param is a vector of length p (number of genes) sampled from values (80, 120, 600) with equal probabilities

modifier <- 1/250
nuisance_param <- 2

n_each <- 50; p_each <- 200
cell_pop <- matrix(c(4,10, 25,100, 60,80, 25,100,
                     40,10, 60,80, 60,80, 100, 25),
                   nrow = 4, ncol = 4, byrow = T)
gene_pop <- matrix(c(20,90, 25,100,
                     90,20, 100,25)/20, nrow = 2, ncol = 4, byrow = T)
sigma <- 5

# this generates cells according to 4 cell types
res <- generate_natural_mat(cell_pop, gene_pop, n_each, p_each, sigma, modifier)
nat_mat <- res$nat_mat

dat <- generator_curved_gaussian(nat_mat, nuisance_param)

# visualize the generated n by p data
clockwise90 <- function(a) { t(a[nrow(a):1,]) }
image(clockwise90(dat))
image(clockwise90(log(dat+1)))

# visualize the true cell embedding and true gene embedding
par(mfrow = c(1,2))
plot(res$cell_mat[,1], res$cell_mat[,2], asp = T, pch = 16, col = rep(1:4, each = n_each),
     main = "True cell embedding", xlab = "Latent dimension 1", ylab = "Latent dimension 2")
plot(res$gene_mat[,1], res$gene_mat[,2], asp = T, pch = 16, col = rep(1:2, each = p_each),
     main = "True gene embedding", xlab = "Latent dimension 1", ylab = "Latent dimension 2")

# try estimating this via the SVD (or anything in simulator_methods.R)
# -- these only return the cell embeddings for consistency across methods, but we can plot the estimated gene embeddings as well in general for eSVD
est_embedding <- method_svd(dat)
plot(est_embedding$fit[,1], est_embedding$fit[,2], asp = T, pch = 16, col = rep(1:4, each = n_each),
     main = "Est. cell embedding", xlab = "Latent dimension 1", ylab = "Latent dimension 2")

