rm(list=ls())

library(Seurat)

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_large.RData")

head(covariates)
nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates[,-2], esvd_res$b_mat[,-2])
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

sns[["pred"]] <- Seurat::CreateAssayObject(counts = t(pred_mat))
Seurat::DefaultAssay(sns) <- "pred"
sns <- Seurat::NormalizeData(sns, normalization.method = "LogNormalize", scale.factor = 10000)
sns[["pred"]]@var.features <- rownames(sns)
sns <- Seurat::ScaleData(sns)
sns <- Seurat::RunPCA(sns, features = colnames(pred_mat), verbose = F,
                      reduction.name = "esvdpca")

set.seed(10)
sns <- Seurat::RunUMAP(sns, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(sns, reduction = "esvdpca", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS, with covariates)\neSVD, Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_esvd_covariates_large_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

rm(list = c("pred_mat")); gc()

######################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(sns@meta.data)[which(sns@meta.data$celltype == "L2/3")]

sns[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(sns, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human brain (SNS, with covariates)\neSVD, Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup7/sns_esvd_covariates_large_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

png("../../../../out/fig/writeup7/sns_esvd_covariates_large_asd_hist.png",
    width = 1800, height = 1500, units = "px", res = 300)
hist(esvd_res$b_mat[,3], main = "Human brain (SNS, with covariates)\neSVD, ASD coef. histogram",
     col = "gray", xlab = "ASD coefficient", breaks = 50)
graphics.off()
