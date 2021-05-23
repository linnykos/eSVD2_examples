rm(list=ls())
dat <- read.table(file = "../../../data/smartseq_mousebrain/mouse_VISp_2018-06-14_exon-matrix.csv",
           sep = ',', stringsAsFactors = FALSE, header = TRUE)
genes <- read.table(file = "../../../data/smartseq_mousebrain/mouse_VISp_2018-06-14_genes-rows.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
rownames(dat) <- make.unique(names = genes$gene_symbol)
dat$X <- NULL
dat <- Matrix::Matrix(as.matrix(dat), sparse = T)
dat <- Matrix::t(dat)

dat[1:5,1:5]
dim(dat)
quantile(dat@x)

metadata <- read.csv(file = "../../../data/smartseq_mousebrain/mouse_VISp_2018-06-14_samples-columns.csv",
                      row.names = 1, stringsAsFactors = FALSE)
table(metadata$subclass)
quantile(sparseMatrixStats::rowSums2(dat))

save.image("../../../data/smartseq_mousebrain/smartseq_mousebrain_formatted.RData")

###################

rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
load("../../../../data/smartseq_mousebrain/smartseq_mousebrain_formatted.RData")

quantile(sparseMatrixStats::rowSums2(dat))
brain <- Seurat::CreateSeuratObject(counts = Matrix::t(dat), meta.data = metadata, min.cells = 10)
low_q_cells <- rownames(brain@meta.data[brain@meta.data$class %in% c('Low Quality', 'No Class'), ])
ok_cells <- rownames(brain@meta.data)[!(rownames(x = brain@meta.data) %in% low_q_cells)]
brain <- brain[, ok_cells]
brain <- Seurat::NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
brain <-  Seurat::FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <-  Seurat::ScaleData(brain)
brain <- Seurat::RunPCA(brain, features = Seurat::VariableFeatures(brain),
                       verbose = F)

set.seed(10)
brain <- Seurat::RunUMAP(brain, dims = 1:50, nneighbors = 5)

plot1 <- Seurat::DimPlot(brain, reduction = "umap", group.by = "subclass", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Brain (Smartseq)")
ggplot2::ggsave(filename = "../../../out/fig/writeup6/smartseq_mousebrain_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")




