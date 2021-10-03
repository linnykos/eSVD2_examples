rm(list=ls())

library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/sns_formatted.RData")
head(sns@meta.data)
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "L2/3")] <- 1
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)

set.seed(10)
sns <- Seurat::SCTransform(sns, verbose = T)
sns <- Seurat::RunPCA(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns, dims = 1:50,
                       reduction.name = 'umap.rna',
                       reduction.key = 'rnaUMAP_')

################

covariates <- c("diagnosis", "sex", "individual", "region", "Capbatch", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "umap.rna",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_default_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "umap.rna")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_default_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

plot1 <- Seurat::FeaturePlot(sns,
                             features = c("TTF2", "MX2", "ASCC1",
                                          "GLRA3", "CIRBP", "SAT2",
                                          "QTRT1", "CDH2", "LUC7L"),
                             reduction = "umap.rna")
ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_default_umap_genes.png"),
                plot1, device = "png", width = 11, height = 8, units = "in")

###################

set.seed(10)
sns_de <- Seurat::FindMarkers(sns,
                              assay = "SCT",
                              slot = "data",
                              ident.1 = "Control",
                              ident.2 = "ASD",
                              group.by = "diagnosis",
                              logfc.threshold = 0,
                              min.pct = 0,
                              min.cells.feature = 0,
                              verbose = T)


