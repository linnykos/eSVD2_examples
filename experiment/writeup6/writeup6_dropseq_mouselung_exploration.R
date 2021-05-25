rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
dat <- anndata::read_h5ad("../../../../data/dropseq_mouselung/lung_regeneration_after_bleo")

dat
head(dat$obs)
head(dat$var)
dim(dat$X); dim(dat$var); dim(dat$obs)
dat$X[1:5,1:5]

lung <- Seurat::CreateSeuratObject(counts = Matrix::t(dat$X))
lung[["celltype"]] <- dat$obs$clusters
lung <- Seurat::NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung <-  Seurat::FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)
lung <-  Seurat::ScaleData(lung)
lung <- Seurat::RunPCA(lung, features = Seurat::VariableFeatures(lung),
                           verbose = F)


diff(lung[["pca"]]@stdev)/lung[["pca"]]@stdev[-1]
set.seed(10)
# however...  70 PCs.. suggested by https://github.com/theislab/cellrank_reproducibility/blob/main/notebooks/fig_6_lung/ML-2020-10-17_fig_6_lung.ipynb
lung <- Seurat::RunUMAP(lung, dims = 1:24)

plot1 <- Seurat::DimPlot(lung, reduction = "umap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/dropseq_mouselung_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
