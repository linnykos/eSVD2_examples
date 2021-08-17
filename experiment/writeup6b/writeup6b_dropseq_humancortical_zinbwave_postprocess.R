rm(list=ls())
load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_zinbwave.RData")

# I'm not sure how to get rid of the normalization
pred_mat <- (1 - t(zinbwave::getPi(zinb_res))) * t(zinbwave::getMu(zinb_res))
colnames(pred_mat) <- colnames(mat)
rownames(pred_mat) <- rownames(mat)

zinb_assay <- Seurat::CreateAssayObject(counts = pred_mat)
cortical[["zinb"]] <- zinb_assay

Seurat::DefaultAssay(cortical) <- "zinb"
cortical <- Seurat::NormalizeData(cortical, normalization.method = "LogNormalize", scale.factor = 10000)
cortical[["zinb"]]@var.features <- rownames(cortical)
cortical <-  Seurat::ScaleData(cortical)
cortical <- Seurat::RunPCA(cortical, features = Seurat::VariableFeatures(cortical), verbose = F,
                     reduction.name = "zinbpca")

set.seed(10)
cortical <- Seurat::RunUMAP(cortical, reduction = "zinbpca", dims = 1:30, reduction.name = "zinbumap")

plot1 <- Seurat::DimPlot(cortical, reduction = "zinbumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\nZINB-WaVE, Full")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_zinbwave_full_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(zinbwave::getW(zinb_res))@cell.embeddings
rownames(tmp) <- rownames(cortical@meta.data)

cortical[["zinbfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(cortical, reduction = "zinbfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\nZINB-WaVE, Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_zinbwave_factor_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
