# from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
bm <- SeuratData::LoadData(ds = "bmcite")

Seurat::DefaultAssay(bm) <- "RNA"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm <-  Seurat::FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)
bm <-  Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, features = Seurat::VariableFeatures(bm),
                          verbose = F)

set.seed(10)
bm <- Seurat::RunUMAP(bm, dims = 1:30)

plot1 <- Seurat::DimPlot(bm, reduction = "umap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

#################
