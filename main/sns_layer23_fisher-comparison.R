rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)

load("../../../out/main/sns_layer23_esvd.RData")

eSVD_obj <- eSVD2:::compute_pvalue(input_obj = eSVD_obj)
fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
esvd_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
rm_idx <- grep("^MT", esvd_genes)
if(length(rm_idx) > 0) esvd_genes <- esvd_genes[-rm_idx]
# pvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
# esvd_genes <- names(pvalue_vec)[order(pvalue_vec, decreasing = T)[1:100]]
# all_genes_vec1 <- names(eSVD_obj$teststat_vec)

load("../../../out/main/sns_layer23_deseq2.RData")
fdr_vec <- stats::p.adjust(deseq2_res$pvalue, method = "BH")
deseq_genes <- rownames(deseq2_res)[which(fdr_vec <= 0.05)]
# pvalue_vec <- -log10(deseq2_res$pvalue)
# deseq_genes <- rownames(deseq2_res)[order(pvalue_vec, decreasing = T)[1:100]]
all_genes_vec2 <- rownames(deseq2_res)
all_genes <- sort(intersect(all_genes_vec1, all_genes_vec2))

sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]
sfari_genes <- sfari_genes[sfari_genes %in% all_genes]
bulk_de_genes <- bulk_de_genes[bulk_de_genes %in% all_genes]

hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
hk_genes <- hk_genes[hk_genes %in% all_genes]

#########################

length(intersect(deseq_genes, sfari_genes))
length(intersect(esvd_genes, sfari_genes))
length(intersect(deseq_genes, sfari_genes))/length(deseq_genes)
length(intersect(esvd_genes, sfari_genes))/length(esvd_genes)

length(intersect(deseq_genes, bulk_de_genes))
length(intersect(esvd_genes, bulk_de_genes))
length(intersect(deseq_genes, bulk_de_genes))/length(deseq_genes)
length(intersect(esvd_genes, bulk_de_genes))/length(esvd_genes)

length(intersect(deseq_genes, hk_genes))
length(intersect(esvd_genes, hk_genes))
length(intersect(deseq_genes, hk_genes))/length(deseq_genes)
length(intersect(esvd_genes, hk_genes))/length(esvd_genes)

##########################

length(intersect(esvd_genes, bulk_de_genes))/length(esvd_genes)
length(bulk_de_genes)/length(all_genes)

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

