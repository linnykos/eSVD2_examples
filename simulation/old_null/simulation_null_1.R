rm(list=ls())
set.seed(10)
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

n_each <- 100
s <- 20
p <- 1000
k <- 2

n <- n_each*s

# form latent variables
x_mat <- cbind(runif(n, min = -0.5, max = 2),
               runif(n, min = -0.5, max = 2))
y_mat <- cbind(runif(p, min = -0.5, max = 2),
               runif(p, min = -0.5, max = 2))

# form covariates
# first form the table
df <- cbind(1,
            0,
            rep(c(0,1), each = s/2),
            rep(c(0,1), times = s/2),
            1:s)
colnames(df) <- c("Intercept", "Log_UMI", "CC", "Sex", "Individual")
# expand to covariate matrix
covariate <- do.call(rbind, lapply(1:nrow(df), function(i){
  matrix(rep(df[i,], each = n_each), nrow = n_each, ncol = ncol(df))
}))
colnames(covariate) <- colnames(df)[1:ncol(df)]
z_mat <- cbind(rep(0, p),
               rep(0, p),
               rep(0, p),
               rep(0, p))
colnames(z_mat) <- colnames(df)[1:(ncol(df)-1)]
# form nuisance
dispersion_vec <- rep(1, p)

# generate data
nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariate[,"CC"], z_mat[,"CC"])
gamma_mat <- matrix(0, nrow = n, ncol = p)
for(j in 1:p){
  gamma_mat[,j] <- stats::rgamma(
    n = n,
    shape = exp(nat_mat[,j])*dispersion_vec[j],
    rate = dispersion_vec[j])
}
gamma_mat <- pmin(gamma_mat, 50)

lib_mat <- tcrossprod(covariate[,c("Intercept", "Log_UMI", "Sex")], z_mat[,c("Intercept", "Log_UMI", "Sex")])
lib_mat <- exp(lib_mat)
obs_mat <- matrix(0, nrow = n, ncol = p)
for(j in 1:p){
  obs_mat[,j] <- stats::rpois(n = n,
                              lambda = lib_mat[,j]*gamma_mat[,j])
}
covariate[,"Log_UMI"] <- log1p(Matrix::rowSums(obs_mat))

rownames(obs_mat) <- paste0("c", 1:nrow(obs_mat))
colnames(obs_mat) <- paste0("g", 1:ncol(obs_mat))
rownames(covariate) <- rownames(obs_mat)

seurat_obj <- Seurat::CreateSeuratObject(counts = t(obs_mat),
                                         meta.data = as.data.frame(covariate[,c("CC", "Sex", "Individual")]))
seurat_obj[["RNA"]]@var.features <- rownames(seurat_obj)
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = F)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:5)

save(seurat_obj,
     covariate,
     df,
     obs_mat,
     x_mat,
     y_mat,
     z_mat,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_null_1.RData")

Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "CC")
Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "Sex")
