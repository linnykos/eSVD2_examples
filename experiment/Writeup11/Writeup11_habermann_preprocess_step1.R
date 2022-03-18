rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- Matrix::readMM("~/project/eSVD/data/GSE135893_habermann_lung/GSE135893_matrix.mtx.gz")
genes <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/GSE135893_genes.tsv.gz",
                  header = F, sep = "\t")
barcodes <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/GSE135893_barcodes.tsv.gz",
                     header = F, sep = "\t")
metadata <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/GSE135893_IPF_metadata.csv.gz",
                     header = T, sep = ",")
table(metadata$population, metadata$celltype)
table(metadata$Status, metadata$Sample_Name)
table(metadata$Status, metadata$Diagnosis)

rownames(mat) <- genes[,1]
colnames(mat) <- barcodes[,1]
all(metadata$X %in% colnames(mat))
mat2 <- mat[,metadata$X]
all(colnames(mat2) == metadata$X)

habermann <- Seurat::CreateSeuratObject(mat2)
habermann$Diagnosis <- metadata$Diagnosis
habermann$Sample_Name <- metadata$Sample_Name
habermann$Sample_Source <- metadata$Sample_Source
habermann$Status <- metadata$Status
habermann$population <- metadata$population
habermann$celltype <- metadata$celltype

habermann[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = habermann, pattern = "^MT-")
habermann[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = habermann, pattern = "^RPS")

set.seed(10)
habermann <- Seurat::SCTransform(habermann, variable.features.n = 3000)
save(habermann, date_of_run, session_info,
     file = "../../../../out/Writeup11/Writeup11_habermann_preprocessed.RData")


