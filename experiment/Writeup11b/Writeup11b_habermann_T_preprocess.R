rm(list=ls())
load("../../../../out/Writeup11/Writeup11_habermann_preprocessed.RData")
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

table(habermann$celltype[which(habermann$celltype == "T Cells")])
keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$celltype == "T Cells")] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)

Seurat::DefaultAssay(habermann) <- "RNA"
habermann[["SCT"]] <- NULL
habermann[["umap.rna"]] <- NULL
set.seed(10)
habermann <- Seurat::SCTransform(habermann, variable.features.n = 3000)

set.seed(10); habermann <- Seurat::RunPCA(habermann, verbose = F)
set.seed(10)
habermann <- Seurat::RunUMAP(habermann, dims = 1:50)

demographic_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Demographics_Information-RawData.csv")
demographic_mat[which(demographic_mat$Sample_Name == "VUHD71"), "Sample_Name"] <- "VUHD071"
length(unique(demographic_mat$Sample_Name)) == length(unique(habermann$Sample_Name))
all(sort(demographic_mat$Sample_Name) == sort(unique(habermann$Sample_Name)))

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

covariates <- c("Status", "Sample_Name", "celltype", "Sample_Source", "Diagnosis",
                "Gender", "Ethnicity", "Tobacco")
for(covariate in covariates){
  plot1 <-Seurat::DimPlot(habermann, reduction = "umap",
                          group.by = covariate, label = TRUE,
                          repel = TRUE, label.size = 2.5,
                          raster = F)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Habermann (T-cells): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  width <- ifelse(covariate == "Sample_Name", 7, 5)
  ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/Writeup11b_habermann_exploration_umap_T_", covariate, ".png"),
                  plot1, device = "png", width = width, height = 5, units = "in")
}

covariates <- c("percent.mt", "Age")
for(covariate in covariates){
  plot1 <-Seurat::FeaturePlot(habermann, reduction = "umap",
                              features = covariate)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Habermann (T-cells): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  width <- 5
  ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/Writeup11b_habermann_exploration_umap_T_", covariate, ".png"),
                  plot1, device = "png", width = width, height = 5, units = "in")
}

save(habermann, date_of_run, session_info,
     file = "../../../../out/Writeup11b/Writeup11b_habermann_T_preprocessed.RData")
