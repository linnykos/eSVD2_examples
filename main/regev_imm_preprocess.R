rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

barcodes <- read.csv("~/nzhanglab/data/SCP259_regev_colitis/Imm.barcodes2.tsv",
                     sep = "\t", header = F)
genes <- read.csv("~/nzhanglab/data/SCP259_regev_colitis/Imm.genes.tsv",
                  sep = "\t", header = F)
mat <- Matrix::readMM("~/nzhanglab/data/SCP259_regev_colitis/gene_sorted-Imm.matrix.mtx")
colnames(mat) <- barcodes[,1]
rownames(mat) <- genes[,1]

regevImm <- Seurat::CreateSeuratObject(counts = mat)

metadata <- read.csv("~/nzhanglab/data/SCP259_regev_colitis/all.meta2.txt",
                     sep = "\t", header = T)
metadata <- metadata[-1,]

n <- ncol(regevImm)
mapping_idx <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  which(metadata[,"NAME"] == colnames(regevImm)[i])
})
regevImm$Celltype <- metadata[mapping_idx,"Cluster"]
regevImm$Subject <- metadata[mapping_idx,"Subject"]
regevImm$Sample_Health <- metadata[mapping_idx,"Health"]
regevImm$Sample_Location <- metadata[mapping_idx,"Location"]
regevImm$Sample <- metadata[mapping_idx,"Sample"]

indiv_metadata <- readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-8.xlsx",
                                    sheet = "Subject_Metadata")
indiv_metadata <- as.data.frame(indiv_metadata)
sample_metadata <- readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-8.xlsx",
                                     sheet = "Sample_Metadata")
sample_metadata <- as.data.frame(sample_metadata)

missing_subjects <- unique(regevImm$Subject)[which(!unique(regevImm$Subject) %in% indiv_metadata[,"Subject ID"])]
if(length(missing_subjects) > 0){
  keep_vec <- rep(1, n)
  for(subj in missing_subjects){
    keep_vec[which(regevImm$Subject == subj)] <- 0
  }
  regevImm$keep <- keep_vec
  regevImm <- subset(regevImm, keep == 1)
}

n <- ncol(regevImm)
disease_vec <- rep(NA, n)
gender_vec <- rep(NA, n)
location_vec <- rep(NA, n)
smoking_vec <- rep(NA, n)
for(indiv in unique(regevImm$Subject)){
  idx1 <- which(regevImm$Subject == indiv)
  idx2 <- which(indiv_metadata[,"Subject ID"] == indiv)
  disease_vec[idx1] <- indiv_metadata[idx2,"Disease"]
  gender_vec[idx1] <- indiv_metadata[idx2,"Gender"]
  location_vec[idx1] <- indiv_metadata[idx2,"Location"]
  smoking_vec[idx1] <- indiv_metadata[idx2,"Smoking"]
}
regevImm$Subject_Disease <- disease_vec
regevImm$Subject_Gender <- gender_vec
regevImm$Subject_Location <- location_vec
regevImm$Subject_Smoking <- smoking_vec

regevImm$Sample_ID <- sapply(regevImm$Sample, function(str){
  strsplit(str, split = "\\.")[[1]][2]
})

########

set.seed(10)
regevImm <- Seurat::NormalizeData(regevImm,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)
regevImm <- Seurat::FindVariableFeatures(regevImm,
                                         selection.method = "vst",
                                         nfeatures = 5000)
regevImm <- Seurat::ScaleData(regevImm)

set.seed(10)
regevImm <- Seurat::RunPCA(regevImm, verbose = F)
set.seed(10)
regevImm <- Seurat::RunUMAP(regevImm, dims = 1:50)
regevImm[["percent.mt"]] <- Seurat::PercentageFeatureSet(regevImm, pattern = "^MT-")

group_vec <- c("Celltype", "Subject", "Sample_Health", "Sample_Location", "Subject_Disease", "Subject_Gender", "Subject_Smoking")
for(group in group_vec){
  plot1 <- Seurat::DimPlot(regevImm, reduction = "umap",
                           group.by = group, label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           raster = FALSE)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Regev Immune Colitis\n", group))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevImm_umap-", group, ".png"),
                  plot1, device = "png", width = 8, height = 5, units = "in")
}

save(regevImm, date_of_run, session_info,
     file = "../../../out/main/regevImm_preprocessed.RData")
