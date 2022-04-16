rm(list=ls())
load("../../../../out/Writeup11b/Writeup11b_adams_T_esvd.RData")
library(Seurat)
source("multiple_testing.R")
source("plotting.R")

mat <- as.matrix(Matrix::t(adams[["RNA"]]@counts[adams[["RNA"]]@var.features,]))
case_control_variable <- "Disease_Identity_IPF"
nuisance_vec[is.na(nuisance_vec)] <- 0
res <- compute_posterior(alpha_max = 50,
                         case_control_variable = case_control_variable,
                         esvd_res = esvd_res_full,
                         mat = mat,
                         nuisance_lower_quantile = 0.01,
                         nuisance_vec = nuisance_vec)

#################

metadata <- adams@meta.data
case_individuals <- unique(metadata[which(metadata$Disease_Identity == "IPF"),"Subject_Identity"])
control_individuals <- unique(metadata[which(metadata$Disease_Identity == "Control"),"Subject_Identity"])
teststat_vec <- compute_test_statistics(case_individuals = case_individuals,
                                        control_individuals = control_individuals,
                                        covariate_individual = "Subject_Identity",
                                        metadata = metadata,
                                        posterior_mean_mat = res$posterior_mean_mat,
                                        posterior_var_mat = res$posterior_var_mat,
                                        verbose = T)

adams_teststat_vec <- teststat_vec
save(adams_teststat_vec,
     file = "../../../../out/Writeup11b/Writeup11b_adams_T_esvd_teststat.RData")

##########

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("B", "Macrophage", "Macrophage Alveolar", "NK"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("Macrophages_disease_vs_control_.csv", "Monocytes_disease_vs_control_.csv",
              "B_Cells_disease_vs_control_.csv", "NK_Cells_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
other_genes <- unique(c(adams_df_genes_others, habermann_df_genes_others))

hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_list <- list(de_genes,
                  setdiff(other_genes, de_genes),
                  setdiff(unique(c(hk_genes, cycling_genes)), c(other_genes, de_genes)))
names(gene_list) <- c("Published DE gene", "Other interest genes", "Housekeeping gene")

png("../../../../out/fig/Writeup11b/adams_T_esvd2_teststat_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
histogram_plot(col_template_vec = c(2,4,3),
               gene_list = gene_list,
               teststat_vec = teststat_vec,
               hist_spacing = 0.5,
               main = "Adams: Histogram of test statistic")
graphics.off()

png("../../../../out/fig/Writeup11b/adams_T_esvd2_teststat_histogram_separate.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
histogram_plot(col_template_vec = c(2,4,3),
               gene_list = gene_list,
               teststat_vec = teststat_vec,
               bool_separate = T,
               hist_spacing = 0.5)
graphics.off()

##########################

vec1 <- esvd_res_full$b_mat[,case_control_variable]
names(vec1) <- colnames(mat)
vec2 <- esvd_init$log_pval_vec
names(vec2) <- colnames(mat)
png("../../../../out/fig/Writeup11b/adams_T_esvd2_additional_histogram.png",
    height = 1000, width = 2000,
    units = "px", res = 300)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))
histogram_plot(col_template_vec = c(2,4,3),
               gene_list = gene_list,
               teststat_vec = vec2,
               hist_spacing = 0.5,
               main = "Adams: Initial p-value")

histogram_plot(col_template_vec = c(2,4,3),
               gene_list = gene_list,
               teststat_vec = vec1,
               hist_spacing = 0.5,
               main = "Adams: Diagnosis coefficient")
graphics.off()

#########################

sing_val <- sqrt(eSVD2:::.l2norm(esvd_res_full$covariate[,case_control_variable]) * eSVD2:::.l2norm(esvd_res_full$b_mat[,case_control_variable]))
tmp <- esvd_res_full$covariate[,case_control_variable]/eSVD2:::.l2norm(esvd_res_full$covariate[,case_control_variable])*sing_val
x_mat <- cbind(esvd_res_full$x_mat, tmp)
set.seed(10)
tmp <- Seurat::RunUMAP(x_mat)@cell.embeddings
rownames(tmp) <- rownames(adams@meta.data)

adams[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                              key = "esvdfactorumap_",
                                                              assay = "RNA")
plot1 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Disease_Identity",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot1 <- plot1 + Seurat::NoLegend()
plot2 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Subject_Identity",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot2 <- plot2 + Seurat::NoLegend()
plot3 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Gender",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot3 <- plot3 + Seurat::NoLegend()
plot4 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Tobacco",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot4 <- plot4 + Seurat::NoLegend()
plot5 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/adams_T_esvd2_umap.png"),
                plot5, device = "png", width = 10, height = 10, units = "in")

####

set.seed(10)
tmp <- Seurat::RunUMAP(esvd_res_full$x_mat)@cell.embeddings
rownames(tmp) <- rownames(adams@meta.data)

adams[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                              key = "esvdfactorumap_",
                                                              assay = "RNA")
plot1 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Disease_Identity",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot1 <- plot1 + Seurat::NoLegend()
plot2 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Subject_Identity",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot2 <- plot2 + Seurat::NoLegend()
plot3 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Gender",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot3 <- plot3 + Seurat::NoLegend()
plot4 <-  Seurat::DimPlot(adams, reduction = "esvdfactorumap",
                          group.by = "Tobacco",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot4 <- plot4 + Seurat::NoLegend()
plot5 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/adams_T_esvd2_umap-noDisease.png"),
                plot5, device = "png", width = 10, height = 10, units = "in")
