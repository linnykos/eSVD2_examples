rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

library(Seurat)

load("../../../../data/dropseq_humancortical/raw_counts_mat.rdata")
metadata <- read.csv("../../../../data/dropseq_humancortical/cell_metadata.csv")
dim(metadata)
table(metadata$Cluster)
dim(raw_counts_mat)

all(metadata$Cell %in% colnames(raw_counts_mat))
n <- nrow(metadata)

# rearrange the columns to match metadata
idx_mapping <- sapply(1:n, function(i){
  which(colnames(raw_counts_mat) == metadata$Cell[i])
})
mat <- raw_counts_mat[,idx_mapping]

# for preprocessing, see http://solo.bmap.ucla.edu/shiny/webapp/
# remove genes with no cells
idx <- which(sparseMatrixStats::rowSums2(mat) == 0)
mat <- mat[-idx,]

# remove cells with too few
mat <- as.matrix(mat)
num_uniq <- apply(mat, 2, function(x){length(which(x != 0))})
quantile(num_uniq)
idx <- which(num_uniq > mean(num_uniq)+3*stats::sd(num_uniq))
mat <- mat[,-idx]
dim(mat)

##############

cortical <- Seurat::CreateSeuratObject(counts = mat)
cortical[["percent.mt"]] <- Seurat::PercentageFeatureSet(cortical, pattern = "^MT-")
cortical <- subset(cortical, subset = percent.mt < 5)
mat <- cortical[["RNA"]]@counts
num_uniq <- apply(mat, 1, function(x){length(which(x != 0))})
idx <- which(num_uniq <= 3)
mat <- mat[-idx,]
cortical <- Seurat::CreateSeuratObject(counts = mat)

Seurat::DefaultAssay(cortical) <- "RNA"
cortical <- Seurat::SCTransform(cortical, verbose = T)
cortical <- Seurat::RunPCA(cortical, verbose = F)
set.seed(10)
cortical <- Seurat::RunUMAP(cortical, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

n <- ncol(cortical)
celltype <- sapply(1:n, function(i){
  idx <- which(metadata$Cell == rownames(cortical@meta.data)[i])[1]
  metadata$Cluster[idx]
})
which(is.na(celltype))
table(celltype)
# replace names, from Figure 1F in https://www.cell.com/action/showPdf?pii=S0896-6273%2819%2930561-6
celltype[celltype == "ExDp1"] <- "Excitatory deep layer 1"
celltype[celltype == "ExDp2"] <- "Excitatory deep layer 2"
celltype[celltype == "ExM"] <- "Maturing excitatory"
celltype[celltype == "ExM-U"] <- "Maturing excitatory upper enriched"
celltype[celltype == "ExN"] <- "Migrating excitatory"
celltype[celltype == "InCGE"] <- "Interneuron CGE"
celltype[celltype == "InMGE"] <- "Interneuron MGE"
celltype[celltype == "PgG2M"] <- "Cycling progenitors (G2/M phase)"
celltype[celltype == "PgS"] <- "Cycling progenitors (S phase)"
celltype[celltype == "Mic"] <- "Microglia"
celltype[celltype == "Per"] <- "Pericyte"
celltype[celltype == "End"] <- "Endothelial"


cortical[["celltype"]] <- celltype


plot1 <- Seurat::DimPlot(cortical, reduction = "umap.rna", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

save(cortical, file = "../../../../data/dropseq_humancortical/dropseq_humancortical_formatted.RData")
