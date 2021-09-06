rm(list=ls())

library(eSVD2); library(glmGamPoi); library(scran)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

print("Loading in data")
dat <- anndata::read_h5ad("../../../../data/10x_mousepancreas/endocrinogenesis_day15.5.h5ad")
tmp <- Matrix::t(dat$X)
tmp <- Matrix::Matrix(as.matrix(tmp), sparse = T)
clusters <- dat$obs$clusters
rm(list = "dat")
gc()

print("Starting Seurat")
pancreas <- Seurat::CreateSeuratObject(counts = tmp)
pancreas[["celltype"]] <- clusters
pancreas <- Seurat::NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas <-  Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)

# extract data
mat <- pancreas[["RNA"]]@counts[Seurat::VariableFeatures(pancreas),]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

binary_mat <- mat
binary_mat[binary_mat > 0] <- 1
tmp <- matrixStats::colSums2(binary_mat)
idx <- which(tmp < nrow(binary_mat)/100)
if(length(idx) > 0){
  mat <- mat[,-idx]
}
rm(list = c("binary_mat", "tmp", "idx"))

set.seed(10)
res <- glmGamPoi::glm_gp(t(mat),
                         design = ~1,
                         col_data = NULL,
                         size_factors = "deconvolution",
                         overdispersion = T,
                         overdispersion_shrinkage = FALSE,
                         do_cox_reid_adjustment = TRUE,
                         subsample = FALSE,
                         on_disk = NULL,
                         verbose = TRUE)
names(res)
quantile(res$overdispersions)
quantile(res$size_factors)
quantile(res$Beta[,1])
dim(res$Mu)
dim(res$Offset)

