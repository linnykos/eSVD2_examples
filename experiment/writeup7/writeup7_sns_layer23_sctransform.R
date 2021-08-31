rm(list=ls())

library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_36501genes.RData")
mat <- as.matrix(mat)

library_size_vec <- rowSums(mat)
covariates <- log(library_size_vec)

n <- nrow(mat)
uniq_diagnos <- unique(metadata$diagnosis)
uniq_sex <- unique(metadata$sex)
uniq_indiv <- unique(metadata$individual)
for(i in uniq_diagnos[-1]){
  tmp <- rep(0, n)
  tmp[which(metadata$diagnosis == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
for(i in uniq_sex[-1]){
  tmp <- rep(0, n)
  tmp[which(metadata$sex == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
for(i in uniq_indiv[-1]){
  tmp <- rep(0, n)
  tmp[which(metadata$individual == i)] <- 1
  covariates <- cbind(covariates, tmp)
}

cell_attr <- as.data.frame(covariates)
colnames(cell_attr) <- c("log_umi", uniq_diagnos[-1], uniq_sex[-1], paste0("i", uniq_indiv[-1]))
cell_attr<- as.matrix(cell_attr)

mat <- as.matrix(Matrix::t(mat))

set.seed(10)
res <- sctransform::vst(umi = mat,
                        cell_attr = cell_attr,
                        latent_var = colnames(cell_attr),
                        n_genes = NULL,
                        verbosity = 2)

# something is sus about this is.nan(cell_attr). see how they do it in the code itself..
