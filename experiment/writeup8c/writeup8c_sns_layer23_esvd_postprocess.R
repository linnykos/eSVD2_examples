rm(list=ls())
library(Seurat)
load("../../../../out/writeup8c/writeup8c_sns_layer23_esvd.RData")

set.seed(10)
umap_res <- Seurat::RunUMAP(esvd_res$x_mat)
umap_res <- umap_res@cell.embeddings
sns[["esvdumap"]] <- Seurat::CreateDimReducObject(embedding = umap_res,
                                                  key = "esvdumap_",
                                                  assay = "RNA")

pred_mat <- exp(tcrossprod(esvd_res$x_mat, esvd_res$y_mat))
pred_mat <- pmin(pred_mat, 100)
colnames(pred_mat) <- paste0("pred-", colnames(pred_mat))
pred_mat <- t(pred_mat)
sns[["pred"]] <- Seurat::CreateAssayObject(counts = pred_mat)

############

covariates <- c("diagnosis", "sex", "individual", "region")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "esvdumap",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_layer23_esvd_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

plot1 <- Seurat::FeaturePlot(sns,
                             features = "nCount_RNA",
                             reduction = "esvdumap")
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): nCount_RNA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_layer23_esvd_umap_nCount_RNA.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

gene_vec <- c("TTF2", "MX2", "ASCC1",
              "GLRA3", "CIRBP", "SAT2",
              "QTRT1", "CDH2", "LUC7L")
gene_vec <- paste0("pred-", gene_vec)
plot1 <- Seurat::FeaturePlot(sns,
                             features = gene_vec,
                             reduction = "esvdumap")
ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_layer23_esvd_umap_genes.png"),
                plot1, device = "png", width = 11, height = 8, units = "in")

