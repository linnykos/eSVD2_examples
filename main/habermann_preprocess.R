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
rownames(mat) <- genes[,1]
colnames(mat) <- barcodes[,1]

mat2 <- mat[,metadata$X]
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
set.seed(10); habermann <- Seurat::RunPCA(habermann, verbose = F)
set.seed(10)
habermann <- Seurat::RunUMAP(habermann, dims = 1:50,
                             reduction.name = 'umap.rna',
                             reduction.key = 'rnaUMAP_')
habermann <- Seurat::CellCycleScoring(habermann, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

##############

demographic_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Demographics_Information-RawData.csv")
demographic_mat[which(demographic_mat$Sample_Name == "VUHD71"), "Sample_Name"] <- "VUHD071"

n <- ncol(habermann)
gender_vec <- rep(NA, n)
age_vec <- rep(NA, n)
ethnicity_vec <- rep(NA, n)
tobacco_vec <- rep(NA, n)
for(i in 1:length(demographic_mat$Sample_Name)){
  subj <- demographic_mat$Sample_Name[i]
  idx <- which(habermann$Sample_Name == subj)
  gender_vec[idx] <- demographic_mat$Gender[i]
  age_vec[idx] <- demographic_mat$Age[i]
  ethnicity_vec[idx] <- demographic_mat$Ethnicity[i]
  tobacco_vec[idx] <- demographic_mat$Tobacco[i]
}
gender_vec <- as.factor(gender_vec)
ethnicity_vec <- as.factor(ethnicity_vec)
tobacco_vec <- as.factor(tobacco_vec)

habermann$Gender <- gender_vec
habermann$Age <- age_vec
habermann$Ethnicity <- ethnicity_vec
habermann$Tobacco <- tobacco_vec

save(habermann, date_of_run, session_info,
     file = "../../../out/main/habermann_preprocessed.RData")

