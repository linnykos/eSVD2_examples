rm(list=ls())
library(Seurat)
library(eSVD2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

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

de_gene_list <- lapply(1:length(file_vec), function(kk){
  file <- file_vec[kk]
  print(file)
  load(file)

  fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
  selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
  gene_names <- names(fdr_vec)

  # see https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
  ego <- clusterProfiler::enrichGO(gene          = selected_genes,
                                   universe      = gene_names,
                                   OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                   keyType       = "SYMBOL",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
  head(ego)
  length(ego$ID)

  print(length(ego$ID))
  if(length(ego$ID) > 5) {

    ## see https://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
    simMatrix <- rrvgo::calculateSimMatrix(ego$ID,
                                           orgdb="org.Hs.eg.db",
                                           ont="BP",
                                           method="Rel")

    scores <- setNames(-log10(ego$qvalue), ego$ID)
    reducedTerms <- rrvgo::reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=0.7,
                                           orgdb="org.Hs.eg.db")

    png(paste0("../../../out/fig/main/sns_", names(file_vec)[kk], "_revigo-treemap.png"),
        height = 2000, width = 2000,
        units = "px", res = 500)
    rrvgo::treemapPlot(reducedTerms)
    graphics.off()
  }
})
