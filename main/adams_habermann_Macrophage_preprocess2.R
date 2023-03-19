rm(list=ls())
library(Seurat)
load("../../../out/main/adams_Macrophage_preprocessed.RData")
load("../../../out/main/habermann_Macrophage_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(adams))
keep_vec[which(adams$Disease_Identity %in% c("Control", "IPF"))] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
tab_vec <- table(adams$Subject_Identity)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(adams))
keep_vec[which(adams$Subject_Identity %in% subj_keep)] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
# table(adams$Subject_Identity, adams$Disease_Identity)

adams <- Seurat::NormalizeData(adams,
                               normalization.method = "LogNormalize", scale.factor = 10000)
adams <- Seurat::FindVariableFeatures(adams,
                                      selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "Macrophage")]

keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$Diagnosis %in% c("Control", "IPF"))] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)
tab_vec <- table(habermann$Sample_Name)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$Sample_Name %in% subj_keep)] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)
# table(habermann$Sample_Name, habermann$Diagnosis)

Seurat::DefaultAssay(habermann) <- "RNA"
habermann <- Seurat::NormalizeData(habermann,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
habermann <- Seurat::FindVariableFeatures(habermann,
                                          selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/Macrophages_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
# length(intersect(habermann_df_genes, adams_df_genes))
# length(unique(c(habermann_df_genes, adams_df_genes)))
# length(intersect(Seurat::VariableFeatures(habermann), Seurat::VariableFeatures(adams)))
# length(unique(c(Seurat::VariableFeatures(habermann), Seurat::VariableFeatures(adams))))

hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

all_genes <- unique(c(Seurat::VariableFeatures(habermann),
                      Seurat::VariableFeatures(adams),
                      habermann_df_genes,
                      adams_df_genes,
                      hk_genes,
                      cycling_genes))
all_available_genes <- intersect(rownames(adams), rownames(habermann))
all_genes <- intersect(all_available_genes, all_genes)

adams[["RNA"]]@var.features <- all_genes
habermann[["RNA"]]@var.features <- all_genes

all(adams[["RNA"]]@var.features %in% rownames(adams))
all(habermann[["RNA"]]@var.features %in% rownames(habermann))

save(adams, date_of_run, session_info,
     file = "../../../out/main/adams_Macrophage_preprocessed.RData")

save(habermann, date_of_run, session_info,
     file = "../../../out/main/habermann_Macrophage_preprocessed.RData")


