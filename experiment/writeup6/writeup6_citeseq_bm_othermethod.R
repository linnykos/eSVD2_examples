# from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
bm <- SeuratData::LoadData(ds = "bmcite")

Seurat::DefaultAssay(bm) <- "RNA"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm <-  Seurat::FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)

mat <- bm[["RNA"]]@counts[Seurat::VariableFeatures(bm),]
set.seed(10)
glmpca_res <- glmpca::glmpca(mat, L = 30, fam = "poi",
                             ctl = list(verbose = T), minibatch = "stochastic")
pred_mat <- glmpca:::predict.glmpca(glmpca_res)

###########

glm_assay <- Seurat::CreateAssayObject(counts = pred_mat)
bm[["glm"]] <- glm_assay

Seurat::DefaultAssay(bm) <- "glm"
bm <- Seurat::NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm <-  Seurat::FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)
bm <-  Seurat::ScaleData(bm)
bm <- Seurat::RunPCA(bm, features = Seurat::VariableFeatures(bm),
                     verbose = F)

set.seed(10)
bm <- Seurat::RunUMAP(bm, dims = 1:30, reduction.name = "glmumap")

plot1 <- Seurat::DimPlot(bm, reduction = "glmumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nGLM-PCA, Full")
ggplot2::ggsave(filename = "../../../out/fig/writeup6/citeseq_bm_glmpca_full_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(bm@meta.data)

bm[["glmfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "glmfactorumap", group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Human bone marrow (CITE-seq)\nGLM-PCA, Factor")
ggplot2::ggsave(filename = "../../../out/fig/writeup6/citeseq_bm_glmpca_factor_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")
