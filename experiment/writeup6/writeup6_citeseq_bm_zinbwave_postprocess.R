rm(list=ls())
load("../../../../out/writeup6/writeup6_citeseq_bm_zinbwave.RData")

# I'm not sure how to get rid of the normalization
pred_mat <- (1 - t(zinbwave::getPi(zinb_res))) * t(zinbwave::getMu(zinb_res))
colnames(pred_mat) <- colnames(mat)
rownames(pred_mat) <- rownames(mat)

zinb_assay <- Seurat::CreateAssayObject(counts = pred_mat)
bm[["zinb"]] <- zinb_assay

Seurat::DefaultAssay(bm) <- "zinb"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm[["zinb"]]@var.features <- rownames(bm)
bm <-  Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, features = Seurat::VariableFeatures(bm), verbose = F,
                     reduction.name = "zinbpca")

set.seed(10)
bm <- Seurat::RunUMAP(bm, reduction = "zinbpca", dims = 1:30, reduction.name = "zinbumap")

plot1 <- Seurat::DimPlot(bm, reduction = "zinbumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nZINB-WaVE, Full")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_zinbwave_full_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(zinbwave::getW(zinb_res))@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["zinbfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "zinbfactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nZINB-WaVE, Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_zinbwave_factor_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
