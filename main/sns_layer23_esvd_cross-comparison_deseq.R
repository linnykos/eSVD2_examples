rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)


load("../../../out/main/sns_layer23_deseq2.RData")
load("../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]

eSVD_obj <- eSVD2:::compute_pvalue(eSVD_obj)

gene_names <- names(eSVD_obj$case_mean)
hk_genes2 <- hk_genes[hk_genes %in% gene_names]
sfari_genes2 <- sfari_genes[sfari_genes %in% gene_names]
bulk_de_genes2 <- bulk_de_genes[bulk_de_genes %in% gene_names]

fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
esvd_selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
esvd_logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
esvd_pthres <- min(esvd_logpvalue_vec[esvd_selected_genes])

deseq_fdr_val <- stats::p.adjust(deseq2_res$pvalue, method = "BH")
names(deseq_fdr_val) <- rownames(deseq2_res)
deseq_selected_genes <- names(deseq_fdr_val)[which(deseq_fdr_val <= 0.05)]
deseq_logpvalue_vec <- -log10(deseq2_res$pvalue)
names(deseq_logpvalue_vec) <- rownames(deseq2_res)
deseq_pthres <- min(deseq_logpvalue_vec[deseq_selected_genes])

# minor adjustments
deseq_logpvalue_vec <- deseq_logpvalue_vec[names(esvd_logpvalue_vec)]
deseq_selected_genes <- intersect(deseq_selected_genes, names(esvd_logpvalue_vec))
deseq_selected_genes <- names(deseq_logpvalue_vec)[which(deseq_logpvalue_vec >= sct_pthres)]
esvd_selected_genes <- names(esvd_logpvalue_vec)[which(esvd_logpvalue_vec >= esvd_pthres)]

length(deseq_selected_genes)
length(esvd_selected_genes)
length(intersect(deseq_selected_genes, esvd_selected_genes))
