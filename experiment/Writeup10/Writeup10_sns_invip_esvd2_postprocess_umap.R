rm(list=ls())
library(Seurat)
load("../../../../out/Writeup10/Writeup10_sns_invip_esvd2.RData")

apply(esvd_res_full$b_mat, 2, quantile)

sns <- Seurat::NormalizeData(sns)
sns <- Seurat::ScaleData(sns, features = sns[["RNA"]]@var.features)
sns <- Seurat::RunPCA(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns, dims = 1:50,
                       reduction.name = 'umap.rna',
                       reduction.key = 'rnaUMAP_')

plot1 <-  Seurat::DimPlot(sns, reduction = "umap.rna",
                          group.by = c("diagnosis", "Seqbatch", "sex", "individual"),
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup10/sns_invip_lognorm_umap.png"),
                plot1, device = "png", width = 12, height = 10, units = "in")

###########################

sing_val <- sqrt(eSVD2:::.l2norm(esvd_res_full$covariate[,"diagnosis_ASD"]) * eSVD2:::.l2norm(esvd_res_full$b_mat[,"diagnosis_ASD"]))
tmp <- esvd_res_full$covariate[,"diagnosis_ASD"]/eSVD2:::.l2norm(esvd_res_full$covariate[,"diagnosis_ASD"])*sing_val
x_mat <- cbind(esvd_res_full$x_mat, tmp)
set.seed(10)
tmp <- Seurat::RunUMAP(x_mat)@cell.embeddings
rownames(tmp) <- rownames(sns@meta.data)

sns[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                        key = "esvdfactorumap_",
                                                        assay = "RNA")
plot1 <-  Seurat::DimPlot(sns, reduction = "esvdfactorumap",
                          group.by = c("diagnosis", "Seqbatch", "sex", "individual"),
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup10/sns_invip_esvd2_umap.png"),
                plot1, device = "png", width = 12, height = 10, units = "in")
