rm(list=ls())

library(SeuratData)
AvailableData()
# seems like it doesn't have 10x v3
# also look at the DuoClustering2018 package

# data from https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3
# vignette from https://www.bioconductor.org/packages/release/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
rhdf5::h5ls("../../../data/10x_pbmc/pbmc_10k_v3_filtered_feature_bc_matrix.h5")
# dat <- rhdf5::h5read("../../../data/10x_pbmc/pbmc_10k_v3_filtered_feature_bc_matrix.h5", "matrix/data")
dat <- Seurat::Read10X_h5("../../../data/10x_pbmc/pbmc_10k_v3_filtered_feature_bc_matrix.h5")
dim(dat)
head(colnames(dat))
head(rownames(dat))
dat[1:5,1:5]
pbmc <- Seurat::CreateSeuratObject(counts = dat, min.cells = 3, min.features = 200)
pbmc

# https://satijalab.org/seurat/articles/sctransform_vignette.html
pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
head(pbmc@meta.data)
pbmc <- Seurat::SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
pbmc <- Seurat::RunPCA(pbmc, verbose = FALSE)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

# [note to self: we should use some cell annotator]
# [see https://bioconductor.org/books/release/OSCA/cell-type-annotation.html]
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- Seurat::FindClusters(pbmc, verbose = FALSE)

png("../../../out/fig/writeup6/10x_v3_pbmc_umap_unlabeled.png", height = 1500, width = 1500, units = "px", res = 300)
Seurat::DimPlot(pbmc, label = TRUE) + Seurat::NoLegend()
graphics.off()


############
# trying the annotation, see https://github.com/dviraran/SingleR/issues/115
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html
pbmc.sce <- Seurat::as.SingleCellExperiment(pbmc, assay = "RNA")
pbmc.sce <- scuttle::logNormCounts(pbmc.sce)
pbmc.sce@assays@data$logcounts[1:5,1:5]
range(pbmc.sce@assays@data$logcounts)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
hpca.se@assays@data$logcounts[1:5,1:5]
range(hpca.se@assays@data$logcounts)

pred.hesc <- SingleR::SingleR(test = pbmc.sce, ref = hpca.se,
                     labels = hpca.se$label.main)
pbmc$SingleR.calls <- pred.hesc$labels

png("../../../out/fig/writeup6/10x_v3_pbmc_umap_unlabeled_singler.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, label = TRUE, group.by = 'SingleR.calls',
                repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("10x PBMC (v3), SingleR HPCA labels")
graphics.off()

