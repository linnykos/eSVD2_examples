rm(list=ls())
load("../../../../out/Writeup11b/Writeup11b_habermann_T_esvd.RData")
library(Seurat)
source("multiple_testing.R")
source("plotting.R")
source("reparameterization.R")

# reparameterize the results
nat_mat <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat) +
  tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
esvd_res_full2 <- final_reparameterization(esvd_res = esvd_res_full,
                                           exclude_vars = colnames(esvd_res_full$covariates)[grep("^Sample_Name", colnames(esvd_res_full$covariates))])
nat_mat2 <- tcrossprod(esvd_res_full2$x_mat, esvd_res_full2$y_mat) +
  tcrossprod(esvd_res_full2$covariates, esvd_res_full2$b_mat)
sum(abs(nat_mat - nat_mat2))
tmp <- crossprod(esvd_res_full2$x_mat); diag(tmp) <- 0
sum(abs(tmp))
tmp <- stats::cor(esvd_res_full2$x_mat); diag(tmp) <- 0
sum(abs(tmp))
tmp <- crossprod(esvd_res_full2$x_mat, esvd_res_full2$covariates[,-grep("^Sample_Name", colnames(esvd_res_full2$covariates))])
sum(abs(tmp))
esvd_res_full <- esvd_res_full2

mat <- as.matrix(Matrix::t(habermann[["RNA"]]@counts[habermann[["RNA"]]@var.features,]))
case_control_variable <- "Diagnosis_IPF"
nuisance_vec[is.na(nuisance_vec)] <- 0
res <- compute_posterior(alpha_max = 50,
                         case_control_variable = case_control_variable,
                         esvd_res = esvd_res_full,
                         mat = mat,
                         nuisance_lower_quantile = 0.01,
                         nuisance_vec = nuisance_vec)

#################

metadata <- habermann@meta.data
case_individuals <- unique(metadata[which(metadata$Diagnosis == "IPF"),"Sample_Name"])
control_individuals <- unique(metadata[which(metadata$Diagnosis == "Control"),"Sample_Name"])
teststat_vec <- compute_test_statistics(case_individuals = case_individuals,
                                        control_individuals = control_individuals,
                                        covariate_individual = "Sample_Name",
                                        metadata = metadata,
                                        posterior_mean_mat = res$posterior_mean_mat,
                                        posterior_var_mat = res$posterior_var_mat,
                                        verbose = T)

habermann_teststat_vec <- teststat_vec
save(habermann_teststat_vec,
     file = "../../../../out/Writeup11b/Writeup11b_habermann_T_esvd_teststat.RData")


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

png("../../../../out/fig/Writeup11b/habermann_T_esvd2_teststat_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
histogram_plot(col_template_vec = c(2,4,3),
               gene_list = gene_list,
               teststat_vec = teststat_vec,
               hist_spacing = 0.5,
               main = "Habermann: Histogram of test statistic")
graphics.off()

png("../../../../out/fig/Writeup11b/habermann_T_esvd2_teststat_histogram_separate.png",
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
png("../../../../out/fig/Writeup11b/habermann_T_esvd2_additional_histogram.png",
    height = 1000, width = 2000,
    units = "px", res = 300)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))
histogram_plot(col_template_vec = c(2,4,3),
               gene_list = gene_list,
               teststat_vec = vec2,
               hist_spacing = 0.5,
               main = "Habermann: Initial p-value")

histogram_plot(col_template_vec = c(2,4,3),
               gene_list = gene_list,
               teststat_vec = vec1,
               hist_spacing = 0.5,
               main = "Habermann: Diagnosis coefficient")
graphics.off()

#########################

sing_val <- sqrt(eSVD2:::.l2norm(esvd_res_full$covariate[,case_control_variable]) * eSVD2:::.l2norm(esvd_res_full$b_mat[,case_control_variable]))
tmp <- esvd_res_full$covariate[,case_control_variable]/eSVD2:::.l2norm(esvd_res_full$covariate[,case_control_variable])*sing_val
x_mat <- cbind(esvd_res_full$x_mat, tmp)
set.seed(10)
tmp <- Seurat::RunUMAP(x_mat)@cell.embeddings
rownames(tmp) <- rownames(habermann@meta.data)

habermann[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                              key = "esvdfactorumap_",
                                                              assay = "RNA")
plot1 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Diagnosis",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot1 <- plot1 + Seurat::NoLegend()
plot2 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Sample_Name",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot2 <- plot2 + Seurat::NoLegend()
plot3 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Gender",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot3 <- plot3 + Seurat::NoLegend()
plot4 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Tobacco",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot4 <- plot4 + Seurat::NoLegend()
plot5 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/habermann_T_esvd2_umap.png"),
                plot5, device = "png", width = 10, height = 10, units = "in")

####

set.seed(10)
tmp <- Seurat::RunUMAP(esvd_res_full$x_mat)@cell.embeddings
rownames(tmp) <- rownames(habermann@meta.data)

habermann[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                              key = "esvdfactorumap_",
                                                              assay = "RNA")
plot1 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Diagnosis",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot1 <- plot1 + Seurat::NoLegend()
plot2 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Sample_Name",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot2 <- plot2 + Seurat::NoLegend()
plot3 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Gender",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot3 <- plot3 + Seurat::NoLegend()
plot4 <-  Seurat::DimPlot(habermann, reduction = "esvdfactorumap",
                          group.by = "Tobacco",
                          label = TRUE,
                          repel = TRUE, label.size = 2.5)
plot4 <- plot4 + Seurat::NoLegend()
plot5 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup11b/habermann_T_esvd2_umap-noDisease.png"),
                plot5, device = "png", width = 10, height = 10, units = "in")

###############

load("../../../../out/Writeup11b/Writeup11b_habermann_T_esvd_teststat.RData")
load("../../../../out/Writeup11b/Writeup11b_adams_T_esvd_teststat.RData")

col_template_vec <- c(2,4,3)
idx_list <- lapply(gene_list, function(gene_vec){
  which(names(habermann_teststat_vec) %in% gene_vec)
})
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(habermann_teststat_vec))
for(i in 1:length(idx_list)){
  col_vec[idx_list[[i]]] <- col_template_vec[i]
}

xlim <- quantile(habermann_teststat_vec, probs = c(0.01, 0.99))*1.1
ylim <- quantile(adams_teststat_vec, probs = c(0.01, 0.99))*1.1

png("../../../../out/fig/Writeup11b/habermann_adams_teststatistic.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(habermann_teststat_vec, adams_teststat_vec,
     xlab = "Habermann", ylab = "Adams", pch = 16, col = col_vec,
     xlim = xlim, ylim = ylim)
graphics.off()

png("../../../../out/fig/Writeup11b/habermann_adams_teststatistic_separate.png",
    height = 1200, width = 3600, res = 300, units = "px")
par(mfrow = c(1,3))
for(i in 1:length(idx_list)){
  plot(habermann_teststat_vec[idx_list[[i]]], adams_teststat_vec[idx_list[[i]]],
       xlab = "Habermann", ylab = "Adams", pch = 16, col = col_template_vec[i],
       xlim = xlim, ylim = ylim)
  lines(c(-15,15), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-15,15), lwd = 2, lty = 2)
}
graphics.off()
