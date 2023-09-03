rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)

load("../../../../out/Writeup13b/sns_layer23_esvd.RData")

fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
esvd_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
all_genes <- names(fdr_vec)

sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]
sfari_genes <- sfari_genes[sfari_genes %in% all_genes]
bulk_de_genes <- bulk_de_genes[bulk_de_genes %in% all_genes]

##########################

m <- length(bulk_de_genes)
n <- length(all_genes) - m
k <- length(esvd_genes)
x <- length(intersect(esvd_genes, bulk_de_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(all_genes)),1) )
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher
