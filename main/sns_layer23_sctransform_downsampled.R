rm(list=ls())

library(Seurat)
load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sns_clean <- sns

sns$region <- factor(sns$region)
sns$sex <- factor(sns$sex)
sns$Seqbatch <- factor(sns$Seqbatch)
set.seed(10)
sns <- Seurat::SCTransform(sns, method = "glmGamPoi",
                             residual.features = sns[["RNA"]]@var.features,
                             vars.to.regress = c("region", "sex", "Seqbatch", "percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours"),
                             verbose = T)
Seurat::Idents(sns) <- "diagnosis"
levels(sns)

Seurat::DefaultAssay(sns) <- "SCT"
de_result <- Seurat::FindMarkers(sns, ident.1 = "ASD", ident.2 = "Control",
                                 slot = "scale.data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)

save(sns, de_result,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_sctransform_downsampled.RData")

###################

downsample_values <- c(0.9, 0.8, 0.7, 0.6, 0.5)
de_result_downsampled <- vector("list", length = length(downsample_values))
names(de_result_downsampled) <- paste0("downsampled_", downsample_values)

for(kk in 1:length(downsample_values)){
  downsample_value <- downsample_values[kk]
  print(paste0("Working on downsample: ", downsample_value))
  if("mat" %in% ls()) rm(list = "mat")

  load(paste0("../../../out/main/sns_layer23_processed_downsampled-", downsample_value, ".RData"))

  metadata <- sns_clean@meta.data[,c("region", "sex", "Seqbatch", "percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours", "diagnosis")]
  sns2 <- Seurat::CreateSeuratObject(counts = Matrix::t(mat), meta.data = metadata)
  sns2[["RNA"]]@var.features <- colnames(mat)

  sns2$region <- factor(sns2$region)
  sns2$sex <- factor(sns2$sex)
  sns2$Seqbatch <- factor(sns2$Seqbatch)
  sns2$diagnosis <- factor(sns2$diagnosis)
  set.seed(10)
  sns2 <- Seurat::SCTransform(sns2, method = "glmGamPoi",
                             residual.features = sns2[["RNA"]]@var.features,
                             vars.to.regress = c("region", "sex", "Seqbatch", "percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours"),
                             verbose = T)
  Seurat::Idents(sns2) <- "diagnosis"
  levels(sns2)

  Seurat::DefaultAssay(sns2) <- "SCT"
  de_result2 <- Seurat::FindMarkers(sns2, ident.1 = "ASD", ident.2 = "Control",
                                   slot = "scale.data",
                                   test.use = "wilcox",
                                   logfc.threshold = 0,
                                   min.pct = 0,
                                   verbose = T)

  de_result_downsampled[[kk]] <- de_result2
}

save(sns, de_result, de_result_downsampled,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_sctransform_downsampled.RData")

##################

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L2/3"),]
de_gene_specific <- tmp[,"Gene name"]
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list[-1], function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
de_genes <- de_genes[!de_genes %in% de_gene_specific]
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.005),"external_gene_name"]

########################

original_selected_genes <- rownames(de_result)[order(de_result$p_val_adj, decreasing = F)[1:100]]
downsampled_selected_genes <- lapply(de_result_downsampled, function(x){
  rownames(x)[order(x$p_val_adj, decreasing = F)[1:100]]
})


sapply(downsampled_selected_genes, function(x){
  length(intersect(x, original_selected_genes))
})

length(intersect(original_selected_genes, sfari_genes))
sapply(downsampled_selected_genes, function(x){
  length(intersect(x, sfari_genes))
})

length(intersect(original_selected_genes, bulk_de_genes))
sapply(downsampled_selected_genes, function(x){
  length(intersect(x, bulk_de_genes))
})

length(intersect(original_selected_genes, hk_genes))
sapply(downsampled_selected_genes, function(x){
  length(intersect(x, hk_genes))
})




