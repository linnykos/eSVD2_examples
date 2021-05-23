rm(list=ls())

load("../../../../out/writeup6/writeup6_citeseq_bm_glmpca_poisson.RData")

glm_assay <- Seurat::CreateAssayObject(counts = pred_mat)
bm[["glm"]] <- glm_assay

Seurat::DefaultAssay(bm) <- "glm"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm[["glm"]]@var.features <- rownames(bm)
bm <-  Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, features = Seurat::VariableFeatures(bm), verbose = F)

set.seed(10)
bm <- Seurat::RunUMAP(bm, dims = 1:30, reduction.name = "glmumap")

plot1 <- Seurat::DimPlot(bm, reduction = "glmumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nGLM-PCA, Full, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_glmpca_full_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "glmfactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nGLM-PCA, Factor, Poisson")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6/citeseq_bm_glmpca_factor_poisson_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
