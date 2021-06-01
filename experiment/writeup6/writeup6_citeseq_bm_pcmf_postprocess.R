rm(list=ls())
load("../../../../out/writeup6/writeup6_citeseq_bm_pcmf.RData")

# Not entirely sure if this is correct...
pred_mat <- pcmf_res$factor$U %*% t(pcmf_res$factor$V)
colnames(pred_mat) <- colnames(mat)
rownames(pred_mat) <- rownames(mat)

pcmf_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
bm[["pcmf"]] <- pcmf_assay

Seurat::DefaultAssay(bm) <- "pcmf"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm[["pcmf"]]@var.features <- rownames(bm)
bm <-  Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, features = Seurat::VariableFeatures(bm), verbose = F,
                     reduction.name = "pcmfpca")

set.seed(10)
bm <- Seurat::RunUMAP(bm, reduction = "pcmfpca", dims = 1:30, reduction.name = "pcmfumap")

plot1 <- Seurat::DimPlot(bm, reduction = "pcmfumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\npCMF, Full")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_pcmf_full_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(pcmf_res$factor$U)@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["pcmffactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "pcmffactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\npCMF, Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_pcmf_factor_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
