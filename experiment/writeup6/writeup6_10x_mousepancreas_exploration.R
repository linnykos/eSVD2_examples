rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
dat <- anndata::read_h5ad("../../../../data/10x_mousepancreas/GSE132188_adata.h5ad.h5")

dat
head(dat$obs)
head(dat$var)
dim(dat$X); dim(dat$var); dim(dat$obs)
dat$X[1:5,1:5]
# seems to only include the preprocessed data, not the raw

#########

# let's look at the data at https://figshare.com/collections/CellRank_for_directed_single-cell_fate_mapping_-_datasets/5172299
dat <- anndata::read_h5ad("../../../../data/10x_mousepancreas/endocrinogenesis_day15.5.h5ad")
dat
head(dat$obs)
head(dat$var)
dim(dat$X); dim(dat$var); dim(dat$obs)
dat$X[1:5,1:5]
# it feels like some of the data isn't loaded in
# for example: dat$layers$spliced... Maybe try loading it in python to extract it?

pancreas <- Seurat::CreateSeuratObject(counts = Matrix::t(dat$X))
pancreas[["celltype"]] <- dat$obs$clusters
pancreas <- Seurat::NormalizeData(pancreas, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas <-  Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)
pancreas <-  Seurat::ScaleData(pancreas)
pancreas <- Seurat::RunPCA(pancreas, features = Seurat::VariableFeatures(pancreas),
                       verbose = F)


diff(pancreas[["pca"]]@stdev)/pancreas[["pca"]]@stdev[-1]
set.seed(10)
# however...  70 PCs.. suggested by https://github.com/theislab/cellrank_reproducibility/blob/main/notebooks/fig_2_pancreas_main/ML-2020-10-14_fig_2_and_3_pancreas_main.ipynb
pancreas <- Seurat::RunUMAP(pancreas, dims = 1:8)

plot1 <- Seurat::DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Pancreas (10x Chromium)")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/10x_mousepancreas_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")


