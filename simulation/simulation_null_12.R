# new idea: let's draw the natrual matrix (aside from covariates) from a gaussian-gaussian model
rm(list=ls())
set.seed(10)
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

n_each <- 100 # number of cells per person
s <- 20 # number of people
p <- 1000 # number of genes

n <- n_each*s # total number of cells

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
               c(rep(0.8,10),rep(0, p-10)), # cc
               rnorm(p, mean = 0, sd = 0.2), # sex
               rnorm(p, mean = 0, sd = 0.5)) # age
# z_mat <- cbind(rep(0, p), # intercept
#                rep(0, p), # library
#                rep(0, p), # cc
#                rep(0, p), # sex
#                rep(0, p)) # age
colnames(z_mat) <- colnames(df)[1:(ncol(df)-1)]

# form nuisance
set.seed(10)
dispersion_vec <- sample(rep(c(10, 1, 0.1), each = ceiling(p/3))[1:p])
dispersion_vec[1:10] <- 1

# construct natural matrix
# first label what type of gene it'll be
gene_type <- rep(NA, p)
gene_type[1:10] <- "true_DE"
gene_type[seq(11, p, by=2)] <- "null_large_var"
gene_type[seq(12, p, by=2)] <- "null_interleaved"
gene_type <- as.factor(gene_type)

nat_mat <- matrix(NA, nrow = n, ncol = p)
rownames(nat_mat) <- 1:n # dummy
colnames(nat_mat) <- paste0("g", 1:p)
ind_idx <- lapply(sort(unique(covariate[,"Individual"])), function(indiv){
  which(covariate[,"Individual"] == indiv)
})
for(i in 1:length(ind_idx)){
  rownames(nat_mat)[ind_idx[[i]]] <- paste0("c", i, "_", ind_idx[[i]])
}
case_indiv <- df[which(df[,"CC"] == 1), "Individual"]
control_indiv <- df[which(df[,"CC"] == 0), "Individual"]

# set true DE
set.seed(10)
gene_idx <- which(gene_type == "true_DE")
cc_diff <- 1
indiv_sd <- 0.1
group_sd <- 0.1
min_mean <- 0.5
max_mean <- 3
for(j in gene_idx){
  case_mean <- stats::runif(1, min = min_mean, max = max_mean)
  control_mean <- case_mean + sample(c(-1,1), size = 1)*cc_diff

  for(indiv in case_indiv){
    mean_val <- stats::rnorm(1, mean = case_mean, sd = group_sd)
    nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                mean = mean_val,
                                                sd = indiv_sd)
  }

  for(indiv in control_indiv){
    mean_val <- stats::rnorm(1, mean = control_mean, sd = group_sd)
    nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                mean = mean_val,
                                                sd = indiv_sd)
  }
}

# set null large variance
set.seed(10)
gene_idx <- which(gene_type == "null_large_var")
cc_diff <- 0.2
indiv_sd <- 1.5
group_sd <- 0.01
min_mean <- 0
max_mean <- 2
for(j in gene_idx){
  print(paste0("===\nGene ", j))
  case_mean <- stats::runif(1, min = min_mean, max = max_mean)
  sign_val <- sample(c(-1,1), size = 1)
  control_mean <- case_mean + sign_val*cc_diff

  for(indiv in case_indiv){
    mean_val <- stats::rnorm(1, mean = case_mean, sd = group_sd)
    # tmp <- abs(stats::rnorm(1, mean = 0, sd = group_sd))
    # mean_val <- case_mean - tmp
    # print(paste0("Case mean: ", round(mean_val, 2)))
    nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                mean = mean_val,
                                                sd = indiv_sd)
  }

  for(indiv in control_indiv){
    mean_val <- stats::rnorm(1, mean = control_mean, sd = group_sd)
    # tmp <- abs(stats::rnorm(1, mean = 0, sd = group_sd))
    # mean_val <- control_mean + tmp
    # print(paste0("Control mean: ", round(mean_val, 2)))
    nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                mean = mean_val,
                                                sd = indiv_sd)
  }
}

# set null interleaved
set.seed(10)
gene_idx <- which(gene_type == "null_interleaved")
indiv_sd_min <- 0.1
indiv_sd_max <- 0.2
group_sd_min <- 0.1
group_sd_max <- 0.2
min_mean <- 0.5
max_mean <- 2
for(j in gene_idx){
  mean_val <- stats::runif(1, min = min_mean, max = max_mean)
  sd_val <- stats::runif(1, min = group_sd_min, max = group_sd_max)

  for(indiv in 1:length(ind_idx)){
    mean_indiv <- stats::rnorm(1, mean = mean_val, sd = sd_val)
    sd_indiv <- stats::runif(1, min = indiv_sd_min, max = indiv_sd_max)

    nat_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                                mean = mean_indiv,
                                                sd = sd_indiv)
  }
}

nat_mat <- pmin(nat_mat, log(100))

# manually force more correlation
# set.seed(10)
# shrink_percentage <- 0.95 # higher means we shrink more. 0.4 "works" but then there is no difference visible
# for(j in 11:round(p*.3)){
#   target_idx <- sample(1:10, 1)
#   tmp_df <- data.frame(x = nat_mat[,target_idx], y = nat_mat[,j])
#   lm_res <- stats::lm(y ~ x - 1, data = tmp_df)
#   new_y <- lm_res$fitted.values
#   nat_mat[,j] <- new_y*shrink_percentage + nat_mat[,j]*(1-shrink_percentage)
# }

# # manually force the off-genes to have no DE
# case_idx <- which(covariate[,"CC"] == 1)
# for(j in 11:p){
#   mean_val_case <- mean(nat_mat[case_idx,j])
#   mean_val_control <- mean(nat_mat[-case_idx,j])
#
#   if(mean_val_control < mean_val_case){
#     nat_mat[case_idx,j] <- nat_mat[case_idx,j] - mean_val_case + mean_val_control
#   } else {
#     nat_mat[-case_idx,j] <- nat_mat[-case_idx,j] - mean_val_control + mean_val_case
#   }
# }

# some global shifting and shrinking
# nat_mat[,(round(p*.3)+1):p] <- (nat_mat[,(round(p*.3)+1):p] - 0.5)*.6
# nat_mat <- nat_mat - 1

gamma_mat <- matrix(0, nrow = n, ncol = p)
set.seed(10)
for(j in 1:p){
  gamma_mat[,j] <- stats::rgamma(
    n = n,
    shape = exp(nat_mat[,j])*dispersion_vec[j],
    rate = dispersion_vec[j])
}
gamma_mat <- pmin(gamma_mat, 50)

lib_mat <- tcrossprod(covariate[,c("Intercept", "Log_UMI", "Sex", "Age")], z_mat[,c("Intercept", "Log_UMI", "Sex", "Age")])
lib_mat <- exp(lib_mat)
obs_mat <- matrix(0, nrow = n, ncol = p)
set.seed(10)
for(j in 1:p){
  obs_mat[,j] <- stats::rpois(n = n,
                              lambda = lib_mat[,j]*gamma_mat[,j])
}
covariate[,"Log_UMI"] <- log1p(Matrix::rowSums(obs_mat))
length(which(obs_mat == 0))/prod(dim(obs_mat))

rownames(obs_mat) <- rownames(nat_mat)
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
     nat_mat,
     obs_mat,
     z_mat,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_null_12.RData")

#########################

Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "CC")
Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "Sex")
Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "Individual")
Seurat::FeaturePlot(seurat_obj, features = "Age")

