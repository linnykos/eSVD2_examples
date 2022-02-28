rm(list=ls())
library(Seurat)
load("../../../../out/Writeup10/Writeup10_sns_invip_esvd2.RData")

apply(esvd_res_full$b_mat, 2, quantile)

tab <- table(sns@meta.data[,"individual"], sns@meta.data[,"diagnosis"])
indiv2 <- sapply(1:nrow(tab), function(i){
  if(tab[i,1] != 0) paste0("A.", rownames(tab)[i]) else paste0("C.", rownames(tab)[i])
})
indiv_vec <- sapply(1:nrow(sns@meta.data), function(i){
  idv <- sns@meta.data[i,"individual"]
  indiv2[which(rownames(tab) == idv)]
})
indiv_vec <- factor(indiv_vec)
sns$individual2 <- indiv_vec

sns <- Seurat::NormalizeData(sns)
sns <- Seurat::ScaleData(sns, features = sns[["RNA"]]@var.features)
sns <- Seurat::RunPCA(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns, dims = 1:50,
                       reduction.name = 'umap.rna',
                       reduction.key = 'rnaUMAP_')

plot1 <-  Seurat::DimPlot(sns, reduction = "umap.rna",
                          group.by = c("diagnosis", "Seqbatch", "sex", "individual2"),
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
                          group.by = c("diagnosis", "Seqbatch", "sex", "individual2"),
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup10/sns_invip_esvd2_umap.png"),
                plot1, device = "png", width = 12, height = 10, units = "in")
