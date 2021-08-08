rm(list=ls())
load("../../../../out/writeup6b/writeup6b_10x_mouseretinal_esvd.RData")

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
retinal[["pred"]] <- glm_assay

Seurat::DefaultAssay(retinal) <- "pred"
retinal <- Seurat::NormalizeData(retinal, normalization.method = "LogNormalize", scale.factor = 10000)
retinal[["pred"]]@var.features <- rownames(retinal)
retinal <-  Seurat::ScaleData(retinal)
retinal <- Seurat::RunPCA(retinal, features = Seurat::VariableFeatures(retinal), verbose = F,
                           reduction.name = "esvdpca")

set.seed(10)
retinal <- Seurat::RunUMAP(retinal, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(retinal, reduction = "esvdpca", group.by = "Cluster", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse retinal (10x)\neSVD, Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/10x_mouseretinal_esvd_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(retinal@meta.data)

retinal[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(retinal, reduction = "esvdfactorumap", group.by = "Cluster", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse retinal (10x)\neSVD, Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/10x_mouseretinal_esvd_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

#########################

png( "../../../../out/fig/writeup6b/10x_mouseretinal_esvd_poisson_libraryhist.png",
     height = 1200, width = 1200, res = 300, units = "px")
graphics::hist(esvd_res$b_mat[,2], col = "gray", main = "Coefficient of log-depth",
               xlab = "Log depth", breaks = 25)
graphics::lines(x = rep(1,2), y = c(0,1e6), col = "red", lwd = 2, lty = 2)
graphics.off()
