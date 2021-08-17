rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

library(Seurat)

load("../../../../data/sns_autism/raw_counts_mat.RData")
metadata <- read.csv("../../../../data/sns_autism/meta.txt", sep = "\t", header = T)

##########

sns <- Seurat::CreateSeuratObject(counts = mat)
sns[["percent.mt"]] <- Seurat::PercentageFeatureSet(sns, pattern = "^MT-")

# the following set is already satisfied
# sns <- subset(sns, subset = nFeature_RNA > 500 & percent.mt < 5)

set.seed(10)
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::SCTransform(sns, verbose = T)
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


save(sns, file = "../../../../data/sns_autism/sns_formatted.RData")

plot1 <- Seurat::DimPlot(sns, reduction = "umap.rna", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

###########

plot1 <- Seurat::DimPlot(sns, reduction = "umap.rna", group.by = "Capbatch", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS)\nCapbatch")
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_umap_capbatch.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

#######

plot1 <- Seurat::DimPlot(sns, reduction = "umap.rna", group.by = "individual", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS)\nIndividual")
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_umap_individual.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


#######

plot1 <- Seurat::DimPlot(sns, reduction = "umap.rna", group.by = "diagnosis", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS)\nDiagnosis")
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_umap_diagnosis.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


