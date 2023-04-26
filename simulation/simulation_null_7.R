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
x_mat <- cbind(runif(n, min = -0.1, max = 1),
               runif(n, min = -0.1, max = 1),
               runif(n, min = -0.1, max = 1))
# have a block structure for the y's
y_block_assignment <- rep(c(1:3), times = ceiling(p/3))[1:p]
y_centers <- matrix(c(1.5, 0.1, 0.1,
                      0.1, 1.5, 0.1,
                      0.1, 0.1, 1.5), ncol = 3, nrow = 3, byrow = T)
y_mat <- y_centers[y_block_assignment,] + matrix(rnorm(p*3, mean = 0, sd = 0.1), ncol = 3, nrow = p)
# image(y_mat[c(which(y_block_assignment == 1), which(y_block_assignment == 2), which(y_block_assignment == 3)),])

# form covariates
# first form the table
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
z_mat <- cbind(rep(0, p), # intercept
               rep(0.1, p), # library
               c(rep(1,10),rep(0, p-10)), # cc
               rnorm(p, mean = 0, sd = 0.2), # sex
               rnorm(p, mean = 0, sd = 0.5)) # age
colnames(z_mat) <- colnames(df)[1:(ncol(df)-1)]
# form nuisance
set.seed(10)
dispersion_vec <- sample(rep(c(0.1, 1, 10), each = ceiling(p/3))[1:p])

# generate data
nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariate[,"CC"], z_mat[,"CC"])
# idx <- c(which(y_block_assignment == 1), which(y_block_assignment == 2), which(y_block_assignment == 3)); image(cor(nat_mat)[idx,idx])

# for each gene, shrink all the cells in an individual to its mean
cell_individual_list <- lapply(1:s, function(i){which(covariate[,"Individual"] == i)})
shrink_percentage <- 0.6 # higher means we shrink more
for(j in 1:p){
  for(idx in cell_individual_list){
    mean_val <- mean(nat_mat[idx,j])
    nat_mat[idx,j] <- mean_val*shrink_percentage + nat_mat[idx,j]*(1-shrink_percentage)
  }
}

# manually force more correlation
set.seed(10)
shrink_percentage <- 0.7 # higher means we shrink more
for(j in 11:p){
  target_idx <- sample(intersect(1:10, which(y_block_assignment == y_block_assignment[j])),1)
  tmp_df <- data.frame(x = nat_mat[,target_idx], y = nat_mat[,j])
  lm_res <- stats::lm(y ~ x - 1, data = tmp_df)
  pred_y <- lm_res$fitted.values
  nat_mat[,j] <- pred_y*shrink_percentage + nat_mat[,j]*(1-shrink_percentage)
}

# manually force the off-genes to have no DE
case_idx <- which(covariate[,"CC"] == 1)
for(j in 11:p){
  mean_val_case <- mean(nat_mat[case_idx,j])
  mean_val_control <- mean(nat_mat[-case_idx,j])

  if(mean_val_control < mean_val_case){
    nat_mat[case_idx,j] - mean_val_case + mean_val_control
  } else {
    nat_mat[-case_idx,j] - mean_val_control + mean_val_case
  }
}

gamma_mat <- matrix(0, nrow = n, ncol = p)
set.seed(10)
for(j in 1:p){
  gamma_mat[,j] <- stats::rgamma(
    n = n,
    shape = exp(nat_mat[,j])*dispersion_vec[j],
    rate = dispersion_vec[j])
}
gamma_mat <- pmin(gamma_mat, 50)
round(quantile(gamma_mat),2)

lib_mat <- tcrossprod(covariate[,c("Intercept", "Log_UMI", "Sex", "Age")], z_mat[,c("Intercept", "Log_UMI", "Sex", "Age")])
lib_mat <- exp(lib_mat)
round(quantile(lib_mat),2)
obs_mat <- matrix(0, nrow = n, ncol = p)
set.seed(10)
for(j in 1:p){
  obs_mat[,j] <- stats::rpois(n = n,
                              lambda = lib_mat[,j]*gamma_mat[,j])
}
covariate[,"Log_UMI"] <- log1p(Matrix::rowSums(obs_mat))
length(which(obs_mat == 0))/prod(dim(obs_mat))
# idx <- c(which(y_block_assignment == 1), which(y_block_assignment == 2), which(y_block_assignment == 3)); image(cor(obs_mat)[idx,idx])

rownames(obs_mat) <- paste0("c", 1:nrow(obs_mat))
colnames(obs_mat) <- paste0("g", 1:ncol(obs_mat))
rownames(covariate) <- rownames(obs_mat)
round(quantile(obs_mat),2)

seurat_obj <- Seurat::CreateSeuratObject(counts = t(obs_mat),
                                         meta.data = as.data.frame(covariate[,c("CC", "Sex", "Age", "Individual")]))
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
     y_block_assignment,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_null_7.RData")

Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "CC")
Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "Sex")
Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "Individual")
Seurat::FeaturePlot(seurat_obj, features = "Age")

