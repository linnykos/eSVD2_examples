rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##########

mat <- Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])

downsample_values <- seq(0.95, 0.6, by = -0.05)
for(downsample_value in downsample_values){
  print(downsample_value)
  x_vec <- mat@x

  set.seed(10)
  x_vec2 <- sapply(x_vec, function(x){
    stats::rbinom(1, size = x, prob = downsample_value)
  })

  mat@x <- as.numeric(x_vec2)
  mat <- Matrix::Matrix(as.matrix(mat), sparse = T)
  save(mat, date_of_run, session_info,
       file = paste0("../../../out/main/sns_layer23_processed_downsampled-", downsample_value, ".RData"))
}
