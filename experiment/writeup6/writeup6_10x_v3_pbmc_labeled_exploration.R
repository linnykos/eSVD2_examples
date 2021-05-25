rm(list=ls())
dat <- read.csv("../../../data/10x_v3_pbmc_labeled/10Xv3/10Xv3_pbmc1.csv",
                row.names = 1)
dat <- Matrix::Matrix(as.matrix(dat), sparse = T)
dim(dat)

celltype <- read.csv("../../../data/10x_v3_pbmc_labeled/10Xv3/10Xv3_pbmc1Labels.csv")
celltype <- as.factor(celltype[,1])
table(celltype)
length(celltype)
quantile(sparseMatrixStats::rowSums2(dat))

save.image("../../../data/10x_v3_pbmc_labeled/10x_v3_pbmc_formatted.RData")

############

rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
load("../../../data/10x_v3_pbmc_labeled/10x_v3_pbmc_formatted.RData")

quantile(sparseMatrixStats::rowSums2(dat))
pbmc <- Seurat::CreateSeuratObject(counts = Matrix::t(dat))
pbmc[["celltype"]] <- celltype
pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <-  Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <-  Seurat::ScaleData(pbmc)
pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc),
                       verbose = F)

diff(pbmc[["pca"]]@stdev)/pbmc[["pca"]]@stdev[-1]
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:30)

plot1 <- Seurat::DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("PBMC (10x v3)")
ggplot2::ggsave(filename = "../../../out/fig/writeup6/10x_v3_pbmc_umap_labeled.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

