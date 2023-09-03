rm(list=ls())
library(Seurat)
library(eSVD2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

load("../../../../out/Writeup13b/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

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

png("../../../../out/fig/Writeup13b/sns_layer23_revigo-treemap.png",
    height = 1500, width = 3000,
    units = "px", res = 500)
rrvgo::treemapPlot(reducedTerms)
graphics.off()
