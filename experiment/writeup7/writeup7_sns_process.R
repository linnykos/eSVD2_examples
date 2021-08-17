rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

library(Seurat)

mat <- Matrix::readMM("../../../../data/sns_autism/matrix.mtx")
gene_mat <- read.csv("../../../../data/sns_autism/genes.tsv", sep = "\t", header = F)
cell_mat <- read.csv("../../../../data/sns_autism/barcodes.tsv", sep = "\t", header = F)
metadata <- read.csv("../../../../data/sns_autism/meta.txt", sep = "\t", header = T)

rownames(mat) <- gene_mat[,2]
colnames(mat) <- cell_mat[,1]
save(mat, file = "../../../../data/sns_autism/raw_counts_mat.RData")





