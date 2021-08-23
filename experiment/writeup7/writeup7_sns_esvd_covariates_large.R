rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

library(Seurat)

load("../../../../data/sns_autism/raw_counts_mat.RData")
metadata <- read.csv("../../../../data/sns_autism/meta.txt", sep = "\t", header = T)

##########

print("Preprocessing")
sns <- Seurat::CreateSeuratObject(counts = mat)
sns[["percent.mt"]] <- Seurat::PercentageFeatureSet(sns, pattern = "^MT-")

# the following set is already satisfied
# sns <- subset(sns, subset = nFeature_RNA > 500 & percent.mt < 5)

set.seed(10)
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::SCTransform(sns, verbose = T, variable.features.n = 10000)
sns <- Seurat::RunPCA(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

sns[["celltype"]] <- metadata$cluster
sns[["sample"]] <- metadata$sample
sns[["individual"]] <- metadata$individual
sns[["region"]] <- metadata$region
sns[["age"]] <- metadata$age
sns[["sex"]] <- metadata$sex
sns[["RNA.Integrity.Number"]] <- metadata$RNA.Integrity.Number
sns[["post.mortem.hours"]] <- metadata$post.mortem.interval..hours.
sns[["diagnosis"]] <- metadata$diagnosis
sns[["Capbatch"]] <- metadata$Capbatch
sns[["Seqbatch"]] <- metadata$Seqbatch
rm(list = "metadata")
save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates_large.RData")

######################

print("Forming covariates")
mat <- sns[["RNA"]]@counts[Seurat::VariableFeatures(sns),]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

set.seed(10)
n <- nrow(mat)
library_size_vec <- rowSums(mat)
covariates <- cbind(matrix(1, nrow = n, ncol = 1), log(library_size_vec))
uniq_diagnos <- unique(sns@meta.data$diagnosis)
# uniq_sex <- unique(sns@meta.data$sex)
for(i in uniq_diagnos[-1]){
  tmp <- rep(0, n)
  tmp[which(sns@meta.data$diagnosis == i)] <- 1
  covariates <- cbind(covariates, tmp)
}
# for(i in uniq_sex[-1]){
#   tmp <- rep(0, n)
#   tmp[which(sns@meta.data$sex == i)] <- 1
#   covariates <- cbind(covariates, tmp)
# }
# colnames(covariates) <- c("Intercept", "Log-library", uniq_diagnos[-1], uniq_sex[-1])
colnames(covariates) <- c("Intercept", "Log-library", uniq_diagnos[-1])

K <- 30
print("Starting initialization")
time_start1 <- Sys.time()
set.seed(10)
init <- eSVD2::initialize_esvd(mat, k = K, family = "poisson", nuisance_param_vec = NA,
                               library_size_vec = 1,
                               covariates = covariates,
                               config = eSVD2::initialization_options(), verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates_large.RData")

print("Starting estimation")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, mat, family = "poisson",
                            nuisance_param_vec = NA, library_size_vec = 1,
                            b_init = init$b_mat, covariates = covariates,
                            max_iter = 100, verbose = 1)
time_end2 <- Sys.time()

save.image("../../../../out/writeup7/writeup7_sns_esvd_covariates_large.RData")


