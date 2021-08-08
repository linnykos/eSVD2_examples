rm(list=ls())
load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_pcmf.RData")

# Not entirely sure if this is correct...
pred_mat <- pcmf_res$factor$U %*% t(pcmf_res$factor$V)
colnames(pred_mat) <- colnames(mat)
rownames(pred_mat) <- rownames(mat)

pcmf_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
cortical[["pcmf"]] <- pcmf_assay

Seurat::DefaultAssay(cortical) <- "pcmf"
cortical <- Seurat::NormalizeData(cortical, normalization.method = "LogNormalize", scale.factor = 10000)
cortical[["pcmf"]]@var.features <- rownames(cortical)
cortical <-  Seurat::ScaleData(cortical)
cortical <- Seurat::RunPCA(cortical, features = Seurat::VariableFeatures(cortical), verbose = F,
                     reduction.name = "pcmfpca")

set.seed(10)
cortical <- Seurat::RunUMAP(cortical, reduction = "pcmfpca", dims = 1:30, reduction.name = "pcmfumap")

plot1 <- Seurat::DimPlot(cortical, reduction = "pcmfumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\npCMF, Full")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_pcmf_full_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(pcmf_res$factor$U)@cell.embeddings
rownames(tmp) <- rownames(cortical@meta.data)

cortical[["pcmffactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(cortical, reduction = "pcmffactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Droq-seq)\npCMF, Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_humancortical_pcmf_factor_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
