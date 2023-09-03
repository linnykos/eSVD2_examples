rm(list=ls())
library(Seurat)
library(eSVD2)

file_vec <- c("../../../out/main/sns_astpp_esvd.RData",
              "../../../out/main/sns_endothelial_esvd.RData",
              "../../../out/main/sns_insst_esvd.RData",
              "../../../out/main/sns_invip_esvd.RData",
              "../../../out/main/sns_layer4_esvd.RData",
              "../../../out/main/sns_layer23_esvd.RData",
              "../../../out/main/sns_layer56_esvd.RData",
              "../../../out/main/sns_layer56cc_esvd.RData",
              "../../../out/main/sns_microglia_esvd.RData",
              "../../../out/main/sns_oligo_esvd.RData",
              "../../../out/main/sns_opc_esvd.RData")
names(file_vec) <- c("astpp", "endothelial", "insst", "invip", "layer4", "layer23",
                     "layer56", "layer56cc", "microglia", "oligo", "opc")

#################

de_gene_list <- lapply(file_vec, function(file){
  print(file)
  load(file)

  eSVD_obj <- eSVD2:::compute_pvalue(input_obj = eSVD_obj)

  fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
  selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]

  list(fdr = fdr_vec,
       genes = sort(selected_genes))
})

load("../../../out/main/sns_astpp_esvd.RData")
all_gene_vec <- names(eSVD_obj$teststat_vec)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.005),"external_gene_name"]

sfari_genes <- sfari_genes[sfari_genes %in% all_gene_vec]
bulk_de_genes <- bulk_de_genes[bulk_de_genes %in% all_gene_vec]

###########################

de_gene_list2 <- lapply(de_gene_list[c("layer4", "layer23", "layer56", "layer56cc")], function(gene_list){
  names(gene_list$fdr)[which(gene_list$fdr <= 0.01)]
})

de_gene_table <- table(unlist(de_gene_list2))
table(de_gene_table)
gene_vec <- names(de_gene_table)[de_gene_table >= 2]

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(sfari_genes)
n <- length(all_gene_vec) - m
k <- length(gene_vec)
x <- length(intersect(gene_vec, sfari_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(all_gene_vec)),1) )

fisher1 <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))

####
m <- length(bulk_de_genes)
n <- length(all_gene_vec) - m
k <- length(gene_vec)
x <- length(intersect(gene_vec, bulk_de_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(all_gene_vec)),1) )

fisher2 <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))

####
tmp <- intersect(bulk_de_genes,sfari_genes)
m <- length(tmp)
n <- length(all_gene_vec) - m
k <- length(gene_vec)
x <- length(intersect(gene_vec, tmp))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(all_gene_vec)),1) )

fisher3 <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))


fisher1; fisher2; fisher3
