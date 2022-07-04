rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../data/sns_autism/sns_formatted2.RData")
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "IN-VIP")] <- 1
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)
mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts))
print(dim(mat))

# remove genes with all 0's
gene_sum <- Matrix::colSums(mat)
if(any(gene_sum == 0)) mat <- mat[,-which(gene_sum == 0)]
print(dim(mat))

# remove cells with too high counts
n <- nrow(mat)
rank_cutoff <- 5
set.seed(10)
high_expression_mat <- sapply(1:ncol(mat), function(j){
  order(mat[,j] + stats::runif(n, min = 0, max = .1), decreasing = T)[1:rank_cutoff]
})
upper_limit <- 4*(rank_cutoff/nrow(mat)*ncol(mat)) # inflation factor of 4
cell_high_expression <- table(as.numeric(high_expression_mat))
remove_cell_idx <- names(cell_high_expression)[cell_high_expression >= upper_limit]
keep_vec <- rep(1, nrow(mat))
keep_vec[as.numeric(remove_cell_idx)] <- 0
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)
mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts))
print(dim(mat))

# remove genes that are too sparsely observed
case_idx <- which(sns$diagnosis == "ASD")
control_idx <- which(sns$diagnosis == "Control")
threshold <- 0.05
gene_bool <- sapply(1:ncol(mat), function(j){
  case_percentage <- length(which(mat[case_idx,j] > 0))/length(case_idx)
  control_percentage <- length(which(mat[control_idx,j] > 0))/length(control_idx)

  max(case_percentage, control_percentage) > threshold
})
mat <- mat[,which(gene_bool)]
print(dim(mat))

# remove genes with too high counts
gene_total <- matrixStats::colSums2(log1p(mat))
gene_threshold <- log1p(10)*nrow(mat)
gene_keep <- colnames(mat)[which(gene_total <= gene_threshold)]
mat <- mat[,gene_keep]
print(dim(mat))

sns <- subset(sns, features = colnames(mat))
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::FindVariableFeatures(sns,
                                    selection.method = "vst",
                                    nfeatures = 5000)

load("../../../data/sns_autism/velmeshev_genes.RData")
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list, function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_keep <- unique(c(de_genes,  hk_genes, sfari_genes, cycling_genes, Seurat::VariableFeatures(sns)))
gene_keep <- rownames(sns[["RNA"]]@counts)[which(rownames(sns[["RNA"]]@counts) %in% gene_keep)]
sns[["RNA"]]@var.features <- gene_keep

################

# keep cells with non-trivial variance
cell_variance <- matrixStats::rowSds(mat)
keep_vec <- rep(1, nrow(mat))
keep_vec[cell_variance < 0.1] <- 0
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)
print(dim(sns))

save(sns, date_of_run, session_info,
     file = "../../../out/main/sns_invip_processed.RData")

