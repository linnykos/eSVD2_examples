rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_esvd.RData")
# load("../../../out/Writeup12/Writeup12_sns_layer23_esvd3.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

source("../experiment/Writeup13b/multtest_custom.R")
eSVD_obj <- eSVD2:::compute_pvalue(eSVD_obj)
fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
gene_names <- names(fdr_vec)

length(selected_genes)
length(gene_names)

############

hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]

hk_genes <- hk_genes[hk_genes %in% gene_names]
sfari_genes <- sfari_genes[sfari_genes %in% gene_names]
bulk_de_genes <- bulk_de_genes[bulk_de_genes %in% gene_names]

length(sfari_genes)
length(bulk_de_genes)
length(intersect(selected_genes, bulk_de_genes))
length(intersect(selected_genes, sfari_genes))

###############

m <- length(bulk_de_genes)
n <- length(gene_names) - m
k <- length(selected_genes)
x <- length(intersect(selected_genes, bulk_de_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(gene_names)),1) )
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher


