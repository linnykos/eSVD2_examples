rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

barcodes <- read.csv("~/nzhanglab/data/SCP259_regev_colitis/Epi.barcodes2.tsv",
                     sep = "\t", header = F)
genes <- read.csv("~/nzhanglab/data/SCP259_regev_colitis/Epi.genes.tsv",
                  sep = "\t", header = F)
mat <- Matrix::readMM("~/nzhanglab/data/SCP259_regev_colitis/gene_sorted-Epi.matrix.mtx")
colnames(mat) <- barcodes[,1]
rownames(mat) <- genes[,1]

regevEpi <- Seurat::CreateSeuratObject(counts = mat)

metadata <- read.csv("~/nzhanglab/data/SCP259_regev_colitis/all.meta2.txt",
                     sep = "\t", header = T)
metadata <- metadata[-1,]

n <- ncol(regevEpi)
mapping_idx <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  which(metadata[,"NAME"] == colnames(regevEpi)[i])
})
regevEpi$Celltype <- metadata[mapping_idx,"Cluster"]
regevEpi$Subject <- metadata[mapping_idx,"Subject"]
regevEpi$Sample_Health <- metadata[mapping_idx,"Health"]
regevEpi$Sample_Location <- metadata[mapping_idx,"Location"]
regevEpi$Sample <- metadata[mapping_idx,"Sample"]

indiv_metadata <- readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-8.xlsx",
                                    sheet = "Subject_Metadata")
indiv_metadata <- as.data.frame(indiv_metadata)
sample_metadata <- readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-8.xlsx",
                                     sheet = "Sample_Metadata")
sample_metadata <- as.data.frame(sample_metadata)

missing_subjects <- unique(regevEpi$Subject)[which(!unique(regevEpi$Subject) %in% indiv_metadata[,"Subject ID"])]
if(length(missing_subjects) > 0){
  keep_vec <- rep(1, n)
  for(subj in missing_subjects){
    keep_vec[which(regevEpi$Subject == subj)] <- 0
  }
  regevEpi$keep <- keep_vec
  regevEpi <- subset(regevEpi, keep == 1)
}

n <- ncol(regevEpi)
disease_vec <- rep(NA, n)
gender_vec <- rep(NA, n)
location_vec <- rep(NA, n)
smoking_vec <- rep(NA, n)
for(indiv in unique(regevEpi$Subject)){
  idx1 <- which(regevEpi$Subject == indiv)
  idx2 <- which(indiv_metadata[,"Subject ID"] == indiv)
  disease_vec[idx1] <- indiv_metadata[idx2,"Disease"]
  gender_vec[idx1] <- indiv_metadata[idx2,"Gender"]
  location_vec[idx1] <- indiv_metadata[idx2,"Location"]
  smoking_vec[idx1] <- indiv_metadata[idx2,"Smoking"]
}
regevEpi$Subject_Disease <- disease_vec
regevEpi$Subject_Gender <- gender_vec
regevEpi$Subject_Location <- location_vec
regevEpi$Subject_Smoking <- smoking_vec

regevEpi$Sample_ID <- sapply(regevEpi$Sample, function(str){
  strsplit(str, split = "\\.")[[1]][2]
})

########

set.seed(10)
regevEpi <- Seurat::NormalizeData(regevEpi,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)
regevEpi <- Seurat::FindVariableFeatures(regevEpi,
                                         selection.method = "vst",
                                         nfeatures = 5000)
regevEpi <- Seurat::ScaleData(regevEpi)

set.seed(10)
regevEpi <- Seurat::RunPCA(regevEpi, verbose = F)
set.seed(10)
regevEpi <- Seurat::RunUMAP(regevEpi, dims = 1:50)
regevEpi[["percent.mt"]] <- Seurat::PercentageFeatureSet(regevEpi, pattern = "^MT-")

group_vec <- c("Celltype", "Subject", "Sample_Health", "Sample_Location", "Subject_Disease", "Subject_Gender", "Subject_Smoking")
for(group in group_vec){
  plot1 <- Seurat::DimPlot(regevEpi, reduction = "umap",
                           group.by = group, label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           raster = FALSE)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Regev Epithelial Colitis\n", group))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_umap-", group, ".png"),
                  plot1, device = "png", width = 8, height = 5, units = "in")
}

save(regevEpi, date_of_run, session_info,
     file = "../../../out/main/regevEpi_preprocessed.RData")
