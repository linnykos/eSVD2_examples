rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- Matrix::readMM("~/project/eSVD/data/GSE136831_adams_lung/GSE136831_RawCounts_Sparse.mtx.gz")
genes <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/GSE136831_AllCells.GeneIDs.txt.gz",
                  header = T, sep = "\t")
barcodes <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/GSE136831_AllCells.cellBarcodes.txt.gz",
                     header = F, sep = "\t")
metadata <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz",
                     header = T, sep = "\t")
rownames(mat) <- genes[,2]
colnames(mat) <- barcodes[,1]

adams <- Seurat::CreateSeuratObject(mat)
adams$CellType_Category <- metadata$CellType_Category
adams$Manuscript_Identity <- metadata$Manuscript_Identity
adams$Subclass_Cell_Identity <- metadata$Subclass_Cell_Identity
adams$Disease_Identity <- metadata$Disease_Identity
adams$Subject_Identity <- metadata$Subject_Identity

adams[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = adams, pattern = "^MT-")
adams[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = adams, pattern = "^RPS")

set.seed(10)
adams <- Seurat::NormalizeData(adams,
                               normalization.method = "LogNormalize", scale.factor = 10000)
adams <- Seurat::FindVariableFeatures(adams,
                                      selection.method = "vst", nfeatures = 5000)
adams <- Seurat::ScaleData(adams)
adams <- Seurat::CellCycleScoring(adams, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

set.seed(10)
adams <- Seurat::RunPCA(adams, verbose = F)
set.seed(10)
adams <- Seurat::RunUMAP(adams, dims = 1:50)

#########

demographic_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/covariates.txt")

n <- ncol(adams)
gender_vec <- rep(NA, n)
age_vec <- rep(NA, n)
ethnicity_vec <- rep(NA, n)
tobacco_vec <- rep(NA, n)
for(i in 1:length(demographic_mat$Subject_ID)){
  subj <- demographic_mat$Subject_ID[i]
  idx <- which(adams$Subject_Identity == subj)
  gender_vec[idx] <- demographic_mat$Sex[i]
  age_vec[idx] <- demographic_mat$Age[i]
  ethnicity_vec[idx] <- demographic_mat$Race[i]
  tobacco_vec[idx] <- demographic_mat$Ever_Smoker[i]
}
gender_vec <- as.factor(gender_vec)
ethnicity_vec <- as.factor(ethnicity_vec)
tobacco_vec <- as.factor(tobacco_vec)
age_vec <- as.numeric(age_vec)

adams$Gender <- gender_vec
adams$Age <- age_vec
adams$Ethnicity <- ethnicity_vec
adams$Tobacco <- tobacco_vec

save(adams, date_of_run, session_info,
     file = "../../../out/main/adams_preprocessed.RData")
