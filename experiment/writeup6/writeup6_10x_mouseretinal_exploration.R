# inspiration from
# https://github.com/berenslab/umi-normalization/blob/main/06_preprocess_retina.ipynb

rm(list=ls())
dat <- read.csv("../../../data/10x_mouseretinal/GSE133382_AtlasRGCs_CountMatrix.csv", stringsAsFactors = FALSE,
                header = TRUE)
rownames(dat) <- dat$X
dat$X <- NULL
dat <- Matrix::Matrix(as.matrix(dat), sparse = T)
dat <- Matrix::t(dat)

metadata <- read.table("../../../data/10x_mouseretinal/RGC_Atlas_coordinates.txt", sep = "\t", header = TRUE, 
                       stringsAsFactors = F)
metadata <- metadata[-1,]

dim(dat)
dim(metadata)
dat[1:5,1:5]
head(metadata)

cell_vec <- rownames(dat)
cell_vec <- gsub(pattern = "\\.", "-", cell_vec)
idx <- which(cell_vec %in% metadata$NAME)
dat <- dat[idx,]
idx <- which(metadata$NAME %in% cell_vec)
metadata <- metadata[idx,]

dim(dat)
dim(metadata)
dat <- dat[order(rownames(dat)),]
rownames(dat) <- gsub(pattern = "\\.", "-", rownames(dat))
metadata <- metadata[order(metadata$NAME),]
all(rownames(dat) == metadata$NAME)

table(metadata$Cluster)
quantile(dat@x)
quantile(sparseMatrixStats::rowSums2(dat))

rm(list=c("idx", "cell_vec"))
save.image("../../../data/10x_mouseretinal/10x_mouseretinal_formatted.RData")

##########################

rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
load("../../../data/10x_mouseretinal/10x_mouseretinal_formatted.RData")

rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

retinal <- Seurat::CreateSeuratObject(counts = Matrix::t(dat), 
                                      meta.data = metadata, min.cells = 5)
largest_batch <- names(which.max(table(retinal@meta.data$BatchID)))
cell_keep <- rownames(retinal@meta.data[retinal@meta.data$BatchID == largest_batch,])
retinal <- retinal[,cell_keep]
retinal <- Seurat::NormalizeData(retinal, normalization.method = "LogNormalize", scale.factor = 10000)
retinal <-  Seurat::FindVariableFeatures(retinal, selection.method = "vst", nfeatures = 2000)
retinal <-  Seurat::ScaleData(retinal)
retinal <- Seurat::RunPCA(retinal, features = Seurat::VariableFeatures(retinal),
                        verbose = F)

diff(retinal[["pca"]]@stdev)/retinal[["pca"]]@stdev[-1]
set.seed(10)
retinal <- Seurat::RunUMAP(retinal, dims = 1:24)

plot1 <- Seurat::DimPlot(retinal, reduction = "umap", group.by = "Cluster", label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse retinal (10x)")
ggplot2::ggsave(filename = "../../../out/fig/writeup6/10x_mouseretinal_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")



